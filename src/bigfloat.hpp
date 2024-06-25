#ifndef BIG_FLOAT
#define BIG_FLOAT
#include <algorithm>
#include <array>
#include <iostream>
#include <ostream>
typedef unsigned long size_t;

typedef union DoubleAndLong {
  double val;
  long binary;
} DoubleAndLong;
template <size_t bytes> struct FixedFloat {
  // it should always have more bytes than a double
  FixedFloat(double init)
      : size_mantissa((size_t)(bytes / 1.3)),
        size_exponent(bytes - size_mantissa) {
    for (size_t i = 0; i < bytes; i++) {
      data[i] = 0;
      if (i < size_mantissa)
        working_mantissa[i] = 0;
    }
    // init data
    DoubleAndLong conv;
    conv.val = init;
    const long init_data = conv.binary;
    if (!is_zero(init)) {
      // first the fraction (6 bytes + 4 bits)
      // exponent (1 byte + 3 bits)
      const long exponent = (init_data & (0x7FFl << 52)) >> 52;
      data[size_mantissa] = exponent & 0xFFl;
      data[size_mantissa + 1] = (exponent & 0x700l) >> 8;
      // now adapt the bias by adding as many binary 1s as the exponential has
      for (int i = 0; i < 52; i++) {
        const size_t d = init_data & (1l << i);
        const size_t bits = size_mantissa * 8 - 52 + i;
        const size_t byte = bits / 8;
        const char bit = bits % 8;
        if (d)
          data[byte] |= 1 << bit;
        else
          data[byte] &= ~(1 << bit);
      }
      // additional digits
      {
        const size_t double_expo = 11;
        char carry = 0;
        for (size_t i = size_mantissa * 8 + double_expo - 1;
             i < (size_mantissa + size_exponent) * 8 - 1; i++) {
          const size_t byte = i / 8;
          const size_t bit = i % 8;
          const unsigned char old_val = (data[byte] & (1 << bit)) == 0 ? 0 : 1;
          // addition
          const char new_val = old_val + carry + 1;
          if (new_val == 1) {
            data[byte] |= (1 << bit); // set to 1 (carry and old_val were 0)
          } else if (new_val == 2) {
            // set to 0 and set carry (either carry or old_val was set)
            data[byte] &= ~(1 << bit); // the digit was 1, we got no carry
                                       // -> it becomes 10, carry is set
            // else the digit was 0, we got a carry and a 1 -> keep carry
            carry = 1;
          } else { // 3
                   // carry is set, value is set
            data[byte] |= (1 << bit);
            carry = 1;
          }
        }
        if (carry == 1) {
          // highest exponent bit becomes 1
          data[size_mantissa + size_exponent - 1] |= 0x80;
        }
      }
    }
    sign = (init_data & (1l << 63)) ? 1 : 0;
  }
  FixedFloat() {
    size_mantissa = (size_t)(bytes / 1.3);
    size_exponent = bytes - size_mantissa;
    for (size_t i = 0; i < bytes; i++) {
      data[i] = 0;
      if (i < size_mantissa)
        working_mantissa[i] = 0;
    }
  }
  FixedFloat<bytes> operator+(FixedFloat<bytes> b) const {
    FixedFloat<bytes> c = FixedFloat<bytes>(*this);
    c += b;
    return c;
  }
  FixedFloat<bytes> operator-(FixedFloat<bytes> b) const {
    FixedFloat<bytes> c = FixedFloat<bytes>(*this);
    c -= b;
    return c;
  }
  void operator+=(FixedFloat<bytes> b) {
    // check if abs of a is > abs b
    if (b.is_zero())
      return; // do nothing
    if (is_zero()) {
      // set to b
      for (size_t i = 0; i < data.size(); i++)
        data[i] = b.data[i];
      sign = b.sign;
      return;
    }
    const long shift = calculate_mantissa_shift(b);
    const long shift_a = shift < 0 ? -shift : 0;
    const long shift_b = shift > 0 ? shift : 0;
    const char sign_a_tmp = sign;
    const char sign_b_tmp = b.sign;
    sign = 0;
    b.sign = 0;
    const bool a_greater = *this > b || *this == b;
    sign = sign_a_tmp;
    b.sign = sign_b_tmp;
    if (shift_a != 0) {
      // copy new exponent of b
      for (size_t i = 0; i < size_exponent; i++)
        data[size_mantissa + i] = b.data[size_mantissa + i];
    }
    // now we can carry out integer addition of mantissa
    // test on sign
    if (sign == b.sign) {
      perform_mantissa_addition(*this, b, shift_a, shift_b, false);
    } else {
      // if b negative -> subtraction
      if (a_greater) // normal a - b
        perform_mantissa_addition(*this, b, shift_a, shift_b, true);
      else { // we calculate b - a and negate
        perform_mantissa_addition(b, *this, shift_b, shift_a, true);
        sign = sign ? 0 : 1;
      }
    }
  }
  void operator-=(FixedFloat<bytes> b) {
    if (b.is_zero())
      return; // do nothing
    if (is_zero()) {
      // set to b
      for (size_t i = 0; i < data.size(); i++)
        data[i] = b.data[i];
      sign = b.sign ? 0 : 1; // flip bs sign
      return;
    }
    const long shift = calculate_mantissa_shift(b);
    const long shift_a = shift < 0 ? -shift : 0;
    const long shift_b = shift > 0 ? shift : 0;
    // check if abs of a is > abs b
    const char sign_a_tmp = sign;
    const char sign_b_tmp = b.sign;
    sign = 0;
    b.sign = 0;
    const bool a_greater = this->operator>(b) || this->operator==(b);
    sign = sign_a_tmp;
    b.sign = sign_b_tmp;
    if (shift_a != 0) {
      // copy new exponent of b
      for (size_t i = 0; i < size_exponent; i++)
        data[size_mantissa + i] = b.data[size_mantissa + i];
    }
    // now we can carry out integer addition of mantissa
    // test on sign
    if (sign == b.sign) {
      if (a_greater) // normal a - b
        perform_mantissa_addition(*this, b, shift_a, shift_b, true);
      else { // we calculate b - a and negate
        perform_mantissa_addition(b, *this, shift_b, shift_a, true);
        sign = sign_a_tmp ? 0 : 1;
      }
    } else {
      // (-a) - b is the same as -(a + b), a - (-b) is the same as a + b (add
      // the abs and keep the sign of a)
      perform_mantissa_addition(*this, b, shift_a, shift_b, false);
    }
  }
  FixedFloat<bytes> operator-() const {
    FixedFloat<bytes> c(*this);
    c.sign = sign ? 0 : 1;
    return c;
  }
  double operator*() const {
    if (is_zero())
      return 0.0;
    // copy the exponent
    DoubleAndLong conv;
    long final = 0;
    // final exponent only has 11 bits, just copy 3 bytes from the original
    for (int i = 0; i < 2; i++)
      final |= ((long)data[size_mantissa + i] << (i * 8 + 52));
    // the old bias was 2^(size_exponent - 1) - 1, the new one is 2^10 -1, i.e.
    // we subtract the bits added in the constructor
    // since only one remaining bit (the 11th) is affected by this no loop is
    // necessary
    if ((data[size_mantissa + 1] & 0x4) !=
        0) {                // the first one that is subtracted
      final &= ~(1l << 62); // set it to 0
    } else
      final |= (1l << 62); // set it to 1
    // just copy as much from the mantissa as possible
    for (long i = size_mantissa * 8 - 52; i < size_mantissa * 8; i++) {
      const size_t byte = i / 8;
      const char bit = i % 8;
      const char d = (data[byte] & (1 << bit)) >> bit;
      final |= (long)d << (i - (size_mantissa * 8 - 52));
    }
    // look at the first missing 3 bits for rounding
    char round = 0;
    for (int i = 0; i < 3; i++) {
      const size_t byte = (size_mantissa * 8 - 53 - i) / 8;
      const char bit = (size_mantissa * 8 - 53 - i) % 8;
      const char d = data[byte] & (1 << bit);
      if (d)
        round |= (1 << (2 - i));
    }
    if (round >= 4)
      final += 1; // TODO handle overflow
    if (sign)
      final |= 1l << 63;
    else
      final &= ~(1l << 63);
    conv.binary = final;
    return conv.val;
  }
  bool operator==(const FixedFloat<bytes> b) const {
    // all bits just have to be the same
    if (b.sign != sign)
      return false;
    // test upper bits of exponent
    for (size_t i = 0; i < std::min(size_exponent, b.size_exponent); i++) {
      if (data[size_mantissa + size_exponent - 1 - i] !=
          b.data[b.size_mantissa + b.size_exponent - 1 - i])
        return false;
    }
    // test upper bits of mantissa
    for (size_t i = 0; i < std::min(size_mantissa, b.size_mantissa); i++) {
      if (data[size_mantissa - 1 - i] != b.data[b.size_mantissa - 1 - i])
        return false;
    }
    return true;
  }
  bool operator!=(const FixedFloat<bytes> b) const { return !(*this == b); }
  bool operator>(FixedFloat<bytes> b) const {
    if (sign && !b.sign)
      return true;
    if (!sign && b.sign)
      return false;
    const bool true_value =
        sign && b.sign ? false
                       : true; // because we have to invert for negative values
    const long shift = calculate_mantissa_shift(b);
    // because value is 1.mantissa, if the exponent is smaller -> the value is
    // smaller
    if (shift != 0)
      return shift > 0 ? true_value : !true_value;
    // a > b if a has a i s.t. i = max{i | a[i] = 1 and b[i] = 0} and b a index
    // j s.t. j = max{j | b[j] = 1 and a[i] = 0} and i > j
    for (long i = size_mantissa * 8 - 1; i >= 0; i--) {
      const size_t byte = i / 8;
      const char bit = i % 8;
      const char data_a = data[byte] & (1 << bit);
      const char data_b = b.data[byte] & (1 << bit);
      if (data_a && !data_b)
        return true_value;
      else if (!data_a && data_b)
        return !true_value;
    }
    // they are the same (always false)
    return false;
  }
  bool operator<(const FixedFloat<bytes> b) const {
    return !(*this > b) && (*this != b);
  }
  // mul
  FixedFloat<bytes> operator*(FixedFloat<bytes> b) const {
    FixedFloat<bytes> res(*this);
    res *= b;
    return res;
  }
  void operator*=(FixedFloat<bytes> b) {
    using namespace std;
    // if exact one is negative, the sign is negative
	int oldsign = sign;
	double vala = **this;
    sign = (sign == b.sign) ? 0 : 1;
    if (is_zero())
      return; // do nothing
    if (b.is_zero()) {
      // set to zero
      for (size_t i = 0; i < data.size(); i++)
        data[i] = 0;
      return;
    }
    // we put the mantissa of a in its working mantissa
    for (size_t i = 0; i < size_mantissa; i++) {
      working_mantissa[i] = data[i];
      data[i] = 0;
    }
    // sum exponents
    char carry = 0;
    for (size_t i = 0; i < size_exponent * 8; i++) {
      const size_t byte = i / 8;
      const char bit = i % 8;
      const char data_a = (data[size_mantissa + byte] & (1l << bit)) >> bit;
      const char data_b = (b.data[size_mantissa + byte] & (1l << bit)) >> bit;
      const char sum = data_a + data_b + carry;
      if (sum % 2 == 0)
        data[size_mantissa + byte] &= ~(1 << bit);
      else
        data[size_mantissa + byte] |= (1 << bit);
      if (sum > 1)
        carry = 1;
      else
        carry = 0;
    }
    // and we have to subtract the bias once
    {
      bool was_carry = carry;
      carry = 0;
      for (size_t i = 0; i < size_exponent * 8; i++) {
        const size_t byte = i / 8;
        const char bit = i % 8;
        const char data_a = (data[size_mantissa + byte] & (1l << bit)) >> bit;
        const char data_b = i < size_exponent * 8 - 1 ? 1 : 0;
        const char sum = data_a - data_b - carry;
        if (sum % 2 == 0 || (data_b == 0 && !was_carry &&
                             sum < 0)) // set to 0 if last bit can't borrow
          data[size_mantissa + byte] &= ~(1 << bit);
        else
          data[size_mantissa + byte] |= (1 << bit);
        if (sum < 0)
          carry = 1;
        else
          carry = 0;
      }
    }
    bool has_one = false; // if the implicit one is present in the result
    // multiply mantissa
    for (size_t i = 0; i < size_mantissa * 8; i++) {
      const size_t byte = i / 8;
      const char bit = i % 8;
      if (working_mantissa[byte] & (1 << bit)) {
        // add complete b, shifted by (i + 1) in reverse order
        const size_t shift_b = size_mantissa * 8 - i;
        char carry = 0;
        // we set carry to 1 if we want to round
        {
          int round = 0;
          for (int k = 0; k < 3; k++) {
            if (shift_b > k) {
              const size_t byte_bj = (shift_b - 1 - k) / 8;
              const char bit_bj = (shift_b - 1 - k) % 8;
              const char data = (b.data[byte_bj] & (1 << bit_bj)) >> bit_bj;
              if (data)
                round |= (1 << (2 - k));
            }
          }
          if (round >= 4)
            carry = 1;
        }
        // add b shifted to the result
        for (size_t j = 0; j < size_mantissa * 8; j++) {
          const size_t byte_aj = j / 8;
          const size_t bit_aj = j % 8;
          const size_t byte_bj = (j + shift_b) / 8;
          const size_t bit_bj = (j + shift_b) % 8;
          char data_bj;
          if (j + shift_b == size_mantissa * 8)
            data_bj = 1;
          else if (j + shift_b < size_mantissa * 8)
            data_bj = ((b.data[byte_bj] & (1 << bit_bj)) >> bit_bj);
          else
            data_bj = 0;
          const char data_aj = (data[byte_aj] & (1 << bit_aj)) >> bit_aj;
          const char sum = data_aj + data_bj + carry;
          if (sum > 1)
            carry = 1;
          else
            carry = 0;
          if (sum % 2 == 0)
            data[byte_aj] &= ~(1 << bit_aj);
          else
            data[byte_aj] |= (1 << bit_aj);
        }
        if (carry && !has_one) {
          has_one = true;
        } else
          // if carry is still set -> we carry to 1., which becomes 10., we
          // shift mantissa one right and add one to exponent
          if (carry && has_one) {
            shift_and_add_exponent(0);
          }
      }
    }
    // add complete mantissa again because 1.
    {
      char carry = 0;
      for (size_t i = 0; i < size_mantissa * 8; i++) {
        const size_t byte = i / 8;
        const size_t bit = i % 8;
        const char data_b = (b.data[byte] & (1 << bit)) >> bit;
        const char data_a = (data[byte] & (1 << bit)) >> bit;
        const char sum = data_a + data_b + carry;
        if (sum > 1)
          carry = 1;
        else
          carry = 0;
        if (sum % 2 == 0)
          data[byte] &= ~(1 << bit);
        else
          data[byte] |= (1 << bit);
      }
      // we have has_one.data + 1.b + carry.0
      const char sum = (has_one ? 1 : 0) + 1 + carry;
      if (sum > 1)
        shift_and_add_exponent((sum % 2 == 0 ? 0 : 1));
      // i have no idea why i would need the following line
      // add_one_to_exponent();
    }
  }

  // protected:
  size_t size_mantissa; // in bytes
  size_t size_exponent; // in bytes
  char sign = 1;
  std::array<char, bytes> data; // exponent, mantissa, i.e. 0 starts at the
                                // mantissa, it indexes as for byte operations
  std::array<char, bytes>
      working_mantissa; // needed as interim memory for multiplication
  /** adds one to the exponent and shifts mantissa one to the right, the
   * inserted bit is configures by `first_digit` */
  void shift_and_add_exponent(char first_digit) {
    // shift mantissa
    for (long j = 1; j < size_mantissa * 8; j++) {
      const size_t byte_aj = j / 8;
      const size_t bit_aj = j % 8;
      const size_t byte_dj = (j - 1) / 8;
      const size_t bit_dj = (j - 1) % 8;
      const char data_a = data[byte_aj] & (1 << bit_aj);
      if (data_a)
        data[byte_dj] |= 1 << bit_dj;
      else
        data[byte_dj] &= ~(1 << bit_dj);
    }
    if (first_digit)
      // set first bit to 1
      data[size_mantissa - 1] |= (1 << 7);
    else
      data[size_mantissa - 1] &= ~(1 << 7);
    // add one to exponent
    add_one_to_exponent();
  }
  void add_one_to_exponent() {
    char carry = 1;
    for (size_t j = 0; j < size_exponent * 8 && carry; j++) {
      const size_t byte_aj = j / 8;
      const size_t bit_aj = j % 8;
      const char data_a =
          (data[size_mantissa + byte_aj] & (1 << bit_aj)) >> bit_aj;
      const char sum = data_a + carry;
      if (sum > 1)
        carry = 1;
      else
        carry = 0;
      if (sum % 2 == 0)
        data[size_mantissa + byte_aj] &= ~(1 << bit_aj);
      else
        data[size_mantissa + byte_aj] |= (1 << bit_aj);
    }
  }
  void sub_one_from_exponent() {
    char carry = 1;
    for (size_t j = 0; j < size_exponent * 8 && carry; j++) {
      const size_t byte_aj = j / 8;
      const size_t bit_aj = j % 8;
      const char data_a =
          (data[size_mantissa + byte_aj] & (1 << bit_aj)) >> bit_aj;
      const char sum = data_a - carry;
      if (sum < 0)
        carry = 1;
      else
        carry = 0;
      if (sum % 2 == 0)
        data[size_mantissa + byte_aj] &= ~(1 << bit_aj);
      else
        data[size_mantissa + byte_aj] |= (1 << bit_aj);
    }
  }
  long calculate_mantissa_shift(const FixedFloat<bytes> &b) const {
    // exponent may differ -> the one with the lower exponent has to be right
    // shifted to match the higher
    // count difference between this and b in exponent. If > 0 -> b has to be
    // right shifted, if < 0 this has to be right shifted
    // we manually subtract a - b in shift_a and b - a in shift_b s.t. we don't
    // have to care about negative numbers
    long shift_a = 0;
    long shift_b = 0;
    char carry_b = 0;
    char carry_a = 0;
    for (size_t i = 0; i < size_exponent * 8; i++) {
      const size_t byte = i / 8;
      const char bit = i % 8;
      const char data_a = (data[size_mantissa + byte] & (1l << bit)) >> bit;
      const char data_b = (b.data[size_mantissa + byte] & (1l << bit)) >> bit;
      const char sum_a = data_a - data_b - carry_a;
      if (sum_a == 0)
        carry_a = 0;
      else if (sum_a == 1)
        shift_a |= (1l << i);
      else if (sum_a == -1) {
        // borrow 2
        shift_a |= (1l << i);
        carry_a = 1;
      } else if (sum_a == -2) {
        carry_a = 1;
      }
      const char sum_b = data_b - data_a - carry_b;
      if (sum_b == 0)
        carry_b = 0;
      else if (sum_b == 1)
        shift_b |= (1l << i);
      else if (sum_b == -1) {
        // borrow 2
        shift_b |= (1l << i);
        carry_b = 1;
      } else if (sum_b == -2) {
        carry_b = 1;
      }
    }
    return carry_a == 0 ? shift_a : -shift_b;
  }
  /**
   * Performs binary addition of mantissas and stores the result in the data of
   * this.
   */
  void perform_mantissa_addition(FixedFloat<bytes> &a, FixedFloat<bytes> &b,
                                 size_t shift_a, size_t shift_b,
                                 bool negate_b) {
    const size_t byte_shift_a = shift_a / 8;
    const size_t byte_shift_b = shift_b / 8;
    const size_t bit_shift_a = shift_a % 8;
    const size_t bit_shift_b = shift_b % 8;
    char carry = 0;
    // we set carry to 1 if we want to round
    {
      int round = 0;
      for (int k = 0; k < 3; k++) {
        if (shift_b > k) {
          const size_t byte_bj = (shift_b - 1 - k) / 8;
          const char bit_bj = (shift_b - 1 - k) % 8;
          const char data = (b.data[byte_bj] & (1 << bit_bj)) >> bit_bj;
          if (data)
            round |= (1 << (2 - k));
        }
      }
      if (round >= 4)
        carry = 1;
    }
    // addition
    for (size_t i = 0; i < size_mantissa * 8; i++) {
      const size_t byte = i / 8;
      const char bit = i % 8;
      char data_a = 0;
      // since it is 1.mantissa, the 1 appears
      if ((byte + byte_shift_a == size_mantissa - 1 &&
           bit + bit_shift_a == 8) ||
          (bit_shift_a == 0 && bit == 0 &&
           byte_shift_a + byte == size_mantissa)) {
        data_a = 1;
      } else if (byte + byte_shift_a < size_mantissa &&
                 (byte + byte_shift_a != size_mantissa - 1 ||
                  bit + bit_shift_a < 8)) {
        const size_t byte_a = byte + byte_shift_a + ((bit + bit_shift_a) / 8);
        const char bit_a = (bit + bit_shift_a) % 8;
        data_a = (a.data[byte_a] & (1 << (bit_a))) >> bit_a;
      }
      char data_b = 0;
      if ((byte + byte_shift_b == size_mantissa - 1 &&
           bit + bit_shift_b == 8) ||
          (bit_shift_b == 0 && bit == 0 &&
           byte_shift_b + byte == size_mantissa)) {
        data_b = 1;
      } else if (byte + byte_shift_b < size_mantissa &&
                 (byte + byte_shift_b != size_mantissa - 1 ||
                  bit + bit_shift_b < 8)) {
        const size_t byte_b = byte + byte_shift_b + ((bit + bit_shift_b) / 8);
        const char bit_b = (bit + bit_shift_b) % 8;
        data_b = (b.data[byte_b] & (1 << (bit_b))) >> bit_b;
      }
      const char val =
          negate_b ? data_a - data_b - carry : data_a + data_b + carry;
      if (val % 2 == 0)
        data[byte] &= ~(1 << bit);
      else
        data[byte] |= (1 << bit);
      if (!negate_b && val > 1 || negate_b && val < 0)
        carry = 1;
      else
        carry = 0;
    }
    if (!negate_b) {
      bool was_carry = carry != 0;
      // both had a 1 before the comma -> increment exponent and shift mantissa
      // to right
      if (shift_a == 0 && shift_b == 0)
        carry = 1;
      if (carry) {
        // shift right
        for (size_t i = 1; i < size_mantissa * 8; i++) {
          const size_t byte = i / 8;
          const char bit = i % 8;
          const char d = (data[byte] & (1 << bit));
          const size_t d_byte = bit == 0 ? byte - 1 : byte;
          const char d_bit = bit == 0 ? 7 : bit - 1;
          if (d)
            data[d_byte] |= (1 << d_bit);
          else
            data[d_byte] &= ~(1 << d_bit);
        }
        // add one to exponent
        for (size_t i = 0; i < size_exponent * 8 && carry != 0; i++) {
          const size_t byte = i / 8;
          const char bit = i % 8;
          const char d = (data[size_mantissa + byte] & (1 << bit));
          if (d) {
            data[size_mantissa + byte] &= ~(1 << bit);
          } else {
            carry = 0;
            data[size_mantissa + byte] |= (1 << bit);
          }
        }
        if (was_carry && shift_a == 0 && shift_b == 0)
          data[size_mantissa - 1] |= (1 << 7);
        else
          data[size_mantissa - 1] &= ~(1 << 7);
      }
    } else {
      if ((shift_a == shift_b && carry == 0) ||
          (shift_a == 0 && shift_b != 0 && carry != 0)) {
        normalize();
      }
    }
  }
  void normalize() {
    // left shift mantissa (and subtract 1 to
    // exponent each) to correct this
    size_t shifting = 1;
    for (long i = size_mantissa * 8 - 1; i >= 0; i--) {
      const int byte = i / 8;
      const char bit = i % 8;
      if (data[byte] & (1 << bit))
        break;
      else
        shifting++;
    }
    // left shift mantissa
    for (long i = size_mantissa * 8 - 1; i >= shifting; i--) {
      const size_t byte = i / 8;
      const char bit = i % 8;
      const size_t byte_sh = (i - shifting) / 8;
      const char bit_sh = (i - shifting) % 8;
      const char d = data[byte_sh] & (1 << bit_sh);
      if (d)
        data[byte] |= 1 << bit;
      else
        data[byte] &= ~(1 << bit);
    }
    // subtract from exponent (the 13092481982367th time i write manual
    // addition here)
    char carry = 0;
    for (int i = 0; i < 64; i++) {
      const int byte = i / 8;
      const char bit = i % 8;
      const char data_a = (data[size_mantissa + byte] & (1 << bit)) >> bit;
      const char data_b = (shifting & (1 << i)) >> i;
      const char sum = data_a - data_b - carry;
      if (sum % 2 == 0)
        data[size_mantissa + byte] &= ~(1 << bit);
      else
        data[size_mantissa + byte] |= (1 << bit);
      if (sum < 0)
        carry = 1;
      else
        carry = 0;
    }
  }
  bool is_zero() const {
    for (char d : data)
      if (d)
        return false;
    return true;
  }
  bool is_zero(double a) const {
    DoubleAndLong conv;
    conv.val = a;
    // set sign to zero
    conv.binary &= ~(1l << 63);
    return conv.binary == 0;
  }
};

#endif

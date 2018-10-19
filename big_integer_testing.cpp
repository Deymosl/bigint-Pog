#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <vector>
#include <utility>
#include <gtest/gtest.h>

#include "big_integer.h"

TEST(correctness, two_plus_two)
{
    EXPECT_EQ(BigInteger(2) + BigInteger(2), BigInteger(4));
    EXPECT_EQ(BigInteger(2) + 2             , 4); // implicit converion from int must work
    EXPECT_EQ(2              + BigInteger(2), 4);
}

TEST(correctness, default_ctor)
{
    BigInteger a;
    BigInteger b = 0;
    EXPECT_EQ(a, 0);
    EXPECT_EQ(a, b);
}

TEST(correctness, ctor_limits)
{
    BigInteger a = std::numeric_limits<int>::min();
    BigInteger b = std::numeric_limits<int>::max();
    EXPECT_EQ(a + b, -1);
}

TEST(correctness, copy_ctor)
{
    BigInteger a = 3;
    BigInteger b = a;

    EXPECT_EQ(a, b);
    EXPECT_EQ(b, 3);
}

TEST(correctness, copy_ctor_real_copy)
{
    BigInteger a = 3;
    BigInteger b = a;
    a = 5;

    EXPECT_EQ(b, 3);
}

TEST(correctness, copy_ctor_real_copy2)
{
    BigInteger a = 3;
    BigInteger b = a;
    b = 5;

    EXPECT_EQ(a, 3);
}

TEST(correctness, assignment_operator)
{
    BigInteger a = 4;
    BigInteger b = 7;
    b = a;

    EXPECT_TRUE(a == b);
}

TEST(correctness, self_assignment)
{
    BigInteger a = 5;
    a = a;

    EXPECT_TRUE(a == 5);
}

TEST(correctness, assignment_return_value)
{
    BigInteger a = 4;
    BigInteger b = 7;
    (a = b) = a;

    EXPECT_TRUE(a == 7);
    EXPECT_TRUE(b == 7);
}

TEST(correctness, comparisons)
{
    BigInteger a = 100;
    BigInteger b = 100;
    BigInteger c = 200;

    EXPECT_TRUE(a == b);
    EXPECT_TRUE(a != c);
    EXPECT_TRUE(a < c);
    EXPECT_TRUE(c > a);
    EXPECT_TRUE(a <= a);
    EXPECT_TRUE(a <= b);
    EXPECT_TRUE(a <= c);
    EXPECT_TRUE(c >= a);
}

TEST(correctness, compare_zero_and_minus_zero)
{
    BigInteger a;
    BigInteger b = -a;

    EXPECT_TRUE(a == b);
}

TEST(correctness, add)
{
    BigInteger a = 5;
    BigInteger b = 20;

    EXPECT_TRUE(a + b == 25);

    a += b;
    EXPECT_TRUE(a == 25);
}

TEST(correctness, add_signed)
{
    BigInteger a = 5;
    BigInteger b = -20;

    EXPECT_TRUE(a + b == -15);

    a += b;
    EXPECT_TRUE(a == -15);
}

TEST(correctness, add_return_value)
{
    BigInteger a = 5;
    BigInteger b = 1;

    (a += b) += b;
    EXPECT_EQ(a, 7);
}

TEST(correctness, sub)
{
    BigInteger a = 20;
    BigInteger b = 5;

    EXPECT_TRUE(a - b == 15);

    a -= b;
    EXPECT_TRUE(a == 15);
}

TEST(correctness, sub_signed)
{
    BigInteger a = 5;
    BigInteger b = 20;

    EXPECT_TRUE(a - b == -15);

    a -= b;
    EXPECT_TRUE(a == -15);

    a -= -100;
    EXPECT_TRUE(a == 85);
}

TEST(correctness, sub_return_value)
{
    BigInteger a = 5;
    BigInteger b = 1;

    (a -= b) -= b;
    EXPECT_EQ(a, 3);
}

TEST(correctness, mul)
{
    BigInteger a = 5;
    BigInteger b = 20;

    EXPECT_TRUE(a * b == 100);

    a *= b;
    EXPECT_TRUE(a == 100);
}

TEST(correctness, mul_signed)
{
    BigInteger a = -5;
    BigInteger b = 20;

    EXPECT_TRUE(a * b == -100);

    a *= b;
    EXPECT_TRUE(a == -100);
}

TEST(correctness, mul_return_value)
{
    BigInteger a = 5;
    BigInteger b = 2;

    (a *= b) *= b;
    EXPECT_EQ(a, 20);
}

TEST(correctness, div_)
{
    BigInteger a = 20;
    BigInteger b = 5;
    BigInteger c = 20;

    EXPECT_TRUE(a / b == 4);
    EXPECT_TRUE(a % b == 0);

    a /= b;
    EXPECT_TRUE(a == 4);

    c %= b;
    EXPECT_TRUE(c == 0);
}

TEST(correctness, div_int_min)
{
    BigInteger a = std::numeric_limits<int>::min();
    EXPECT_TRUE((a / a) == (a / std::numeric_limits<int>::min()));
}

TEST(correctness, div_int_min_2)
{
    BigInteger a = std::numeric_limits<int>::min();
    BigInteger b = -1;
    BigInteger c = a / b;
    EXPECT_TRUE(c == (a / -1));
    EXPECT_TRUE((c - std::numeric_limits<int>::max()) == 1);
}

TEST(correctness, div_signed)
{
    BigInteger a = -20;
    BigInteger b = 5;

    EXPECT_TRUE(a / b == -4);
    EXPECT_TRUE(a % b == 0);
}

TEST(correctness, div_rounding)
{
    BigInteger a = 23;
    BigInteger b = 5;

    EXPECT_TRUE(a / b == 4);
    EXPECT_TRUE(a % b == 3);
}

TEST(correctness, div_rounding_negative)
{
    BigInteger a = 23;
    BigInteger b = -5;
    BigInteger c = -23;
    BigInteger d = 5;

    EXPECT_TRUE(a / b == -4);
    EXPECT_TRUE(c / d == -4);
    EXPECT_TRUE(a % b == 3);
    EXPECT_TRUE(c % d == -3);
}

TEST(correctness, div_return_value)
{
    BigInteger a = 100;
    BigInteger b = 2;

    (a /= b) /= b;
    EXPECT_EQ(a, 25);
}

TEST(correctness, unary_plus)
{
    BigInteger a = 123;
    BigInteger b = +a;

    EXPECT_TRUE(a == b);

    // this code should not compile:
    // &+a;
}

TEST(correctness, negation)
{
    BigInteger a = 666;
    BigInteger b = -a;

    EXPECT_TRUE(b == -666);
}

TEST(correctness, negation_int_min)
{
    BigInteger a = std::numeric_limits<int>::min();
    BigInteger b = -a;

    EXPECT_EQ(std::numeric_limits<int>::max(), b - 1);
}

TEST(correctness, and_)
{
    BigInteger a = 0x55;
    BigInteger b = 0xaa;

    EXPECT_TRUE((a & b) == 0);
    EXPECT_TRUE((a & 0xcc) == 0x44);
    a &= b;
    EXPECT_TRUE(a == 0);
}

TEST(correctness, and_signed)
{
    BigInteger a = 0x55;
    BigInteger b = 0xaa;

    EXPECT_TRUE((b & -1) == 0xaa);
    EXPECT_TRUE((a & (0xaa - 256)) == 0);
    EXPECT_TRUE((a & (0xcc - 256)) == 0x44);

    BigInteger c = 0x55;
    BigInteger d = 0xcc;
    EXPECT_EQ(BigInteger(0x44), c & d);
}

TEST(correctness, and_return_value)
{
    BigInteger a = 7;

    (a &= 3) &= 6;
    EXPECT_EQ(a, 2);
}

TEST(correctness, or_)
{
    BigInteger a = 0x55;
    BigInteger b = 0xaa;

    EXPECT_TRUE((a | b) == 0xff);
    a |= b;
    EXPECT_TRUE(a == 0xff);

    BigInteger c = 0x55;
    BigInteger d = 0xcc;
    EXPECT_EQ(BigInteger(0xdd), c | d);
}

TEST(correctness, or_signed)
{
    BigInteger a = 0x55;
    BigInteger b = 0xaa;

    EXPECT_TRUE((a | (b - 256)) == -1);
}

TEST(correctness, or_return_value)
{
    BigInteger a = 1;

    (a |= 2) |= 4;
    EXPECT_EQ(a, 7);
}

TEST(correctness, xor_)
{
    BigInteger a = 0xaa;
    BigInteger b = 0xcc;

    EXPECT_TRUE((a ^ b) == 0x66);

    BigInteger c = 0x55;
    BigInteger d = 0xcc;
    EXPECT_EQ(BigInteger(0x99), c ^ d);
}

TEST(correctness, xor_signed)
{
    BigInteger a = 0xaa;
    BigInteger b = 0xcc;

    EXPECT_TRUE((a ^ (b - 256)) == (0x66 - 256));
}

TEST(correctness, xor_return_value)
{
    BigInteger a = 1;

    (a ^= 2) ^= 1;
    EXPECT_EQ(a, 2);
}

TEST(correctness, not_)
{
    BigInteger a = 0xaa;

    EXPECT_TRUE(~a == (-a - 1));
}

TEST(correctness, shl_)
{
    BigInteger a = 23;

    a << 5;

    std::cout << "ok";

    EXPECT_TRUE((a << 5) == 23 * 32);

    a <<= 5;
    EXPECT_TRUE(a == 23 * 32);
}

TEST(correctness, shl_return_value)
{
    BigInteger a = 1;

    (a <<= 2) <<= 1;
    EXPECT_EQ(a, 8);
}

TEST(correctness, shr_)
{
    BigInteger a = 23;

    EXPECT_EQ(a >> 2, 5);

    a >>= 2;
    EXPECT_EQ(a, 5);
}

TEST(correctness, shr_31)
{
    BigInteger a = 65536;

    EXPECT_EQ((a*a) >> 31, 2);
}

TEST(correctness, shr_signed)
{
    BigInteger a = -1234;

    EXPECT_EQ(a >> 3, -155);

    a >>= 3;
    EXPECT_EQ(a, -155);
}

TEST(correctness, shr_return_value)
{
    BigInteger a = 64;

    (a >>= 2) >>= 1;
    EXPECT_EQ(a, 8);
}

TEST(correctness, add_long)
{
    BigInteger a("10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
    BigInteger b(                                                     "100000000000000000000000000000000000000");
    BigInteger c("10000000000000000000000000000000000000000000000000000100000000000000000000000000000000000000");

    EXPECT_EQ(a + b, c);
}

TEST(correctness, add_long_signed)
{
    BigInteger a("-1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
    BigInteger b( "1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");

    EXPECT_EQ(a + b, 0);
}

TEST(correctness, add_long_signed2)
{
    BigInteger a("-1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
    BigInteger b(                                                     "100000000000000000000000000000000000000");
    BigInteger c( "-999999999999999999999999999999999999999999999999999900000000000000000000000000000000000000");

    EXPECT_EQ(a + b, c);
}

TEST(correctness, add_long_pow2)
{
    BigInteger a( "18446744073709551616");
    BigInteger b("-18446744073709551616");
    BigInteger c( "36893488147419103232");

    EXPECT_EQ(a + a, c);
    EXPECT_EQ(b + c, a);
    EXPECT_EQ(c + b, a);
}

TEST(correctness, sub_long)
{
    BigInteger a("10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
    BigInteger b(                                                     "100000000000000000000000000000000000000");
    BigInteger c( "9999999999999999999999999999999999999999999999999999900000000000000000000000000000000000000");

    EXPECT_EQ(a - b, c);
}

TEST(correctness, mul_long)
{
    BigInteger a("10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
    BigInteger b(                                                     "100000000000000000000000000000000000000");
    BigInteger c("10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
                  "00000000000000000000000000000000000000");

    EXPECT_EQ(a * b, c);
}

TEST(correctness, mul_long_signed)
{
    BigInteger a("-1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
    BigInteger b(                                                     "100000000000000000000000000000000000000");
    BigInteger c("-1000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000"
                  "00000000000000000000000000000000000000");

    EXPECT_EQ(a * b, c);
}

TEST(correctness, mul_long_signed2)
{
    BigInteger a("-100000000000000000000000000");
    BigInteger c("100000000000000000000000000"
                  "00000000000000000000000000");

    EXPECT_EQ(a * a, c);
}

TEST(correctness, mul_long_pow2)
{
    BigInteger a("18446744073709551616");
    BigInteger b("340282366920938463463374607431768211456");
    BigInteger c("115792089237316195423570985008687907853269984665640564039457584007913129639936");

    EXPECT_EQ(a * a, b);
    EXPECT_EQ(b * b, c);
}

TEST(correctness, div_long)
{
    BigInteger a("10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
    BigInteger b(                                                     "100000000000000000000000000000000000000");
    BigInteger c("100000000000000000000000000000000000000000000000000000");

    EXPECT_EQ(a / b, c);
}

TEST(correctness, div_long_signed)
{
    BigInteger a("-10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
    BigInteger b(                                                      "100000000000000000000000000000000000000");
    BigInteger c("-100000000000000000000000000000000000000000000000000000");

    EXPECT_EQ(a / b, c);
}

TEST(correctness, div_long_signed2)
{
    BigInteger a("-10000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000");
    BigInteger b(                                                     "-100000000000000000000000000000000000000");
    BigInteger c( "100000000000000000000000000000000000000000000000000000");

    EXPECT_EQ(a / b, c);
}

TEST(correctness, negation_long)
{
    BigInteger a( "10000000000000000000000000000000000000000000000000000");
    BigInteger c("-10000000000000000000000000000000000000000000000000000");

    EXPECT_EQ(-a, c);
    EXPECT_EQ(a, -c);
}

TEST(correctness, string_conv)
{
    EXPECT_EQ(to_string(BigInteger("100")), "100");
    EXPECT_EQ(to_string(BigInteger("0100")), "100");
    EXPECT_EQ(to_string(BigInteger("0")), "0");
    EXPECT_EQ(to_string(BigInteger("-0")), "0");
    EXPECT_EQ(to_string(BigInteger("-1000000000000000")), "-1000000000000000");
}


namespace
{
    unsigned const number_of_iterations = 1;
    size_t const number_of_multipliers = 10;

    int myrand()
    {
        int val = rand() - RAND_MAX / 2;
        if (val != 0)
            return val;
        else
            return 1;
    }
}

TEST(correctness, mul_div_randomized)
{
    for (unsigned itn = 0; itn != number_of_iterations; ++itn)
    {
        std::vector<int> multipliers;

        for (size_t i = 0; i != number_of_multipliers; ++i)
            multipliers.push_back(myrand());

        BigInteger accumulator = 1;

        for (size_t i = 0; i != number_of_multipliers; ++i)
            accumulator *= multipliers[i];

        std::random_shuffle(multipliers.begin(), multipliers.end());

        for (size_t i = 1; i != number_of_multipliers; ++i)
            accumulator /= multipliers[i];

        EXPECT_TRUE(accumulator == multipliers[0]);
    }
}

namespace
{
    template <typename T>
    void erase_unordered(std::vector<T>& v, typename std::vector<T>::iterator pos)
    {
        std::swap(v.back(), *pos);
        v.pop_back();
    }

    template <typename T>
    T extract_random_element(std::vector<T>& v)
    {
        size_t index = rand() % v.size();
        T copy = v[index];
        erase_unordered(v, v.begin() + index);
        return copy;
    }

    template <typename T>
    void merge_two(std::vector<T>& v)
    {
        assert(v.size() >= 2);

        T a = extract_random_element(v);
        T b = extract_random_element(v);

        T ab = a * b;
        ASSERT_TRUE(ab / a == b);
        ASSERT_TRUE(ab / b == a);

        v.push_back(ab);
    }

    template <typename T>
    T merge_all(std::vector<T> v)
    {
        assert(!v.empty());

        while (v.size() >= 2)
            merge_two(v);

        return v[0];
    }
}

TEST(correctness, mul_merge_randomized)
{
    for (unsigned itn = 0; itn != number_of_iterations; ++itn)
    {
        std::vector<BigInteger> x;
        for (size_t i = 0; i != number_of_multipliers; ++i)
            x.push_back(myrand());

        BigInteger a = merge_all(x);
        BigInteger b = merge_all(x);

        EXPECT_TRUE(a == b);
    }
}

namespace
{
    BigInteger rand_big(size_t size)
    {
        BigInteger result = rand();

        for (size_t i = 0; i != size; ++i)
        {
            result *= RAND_MAX;
            result += rand();
        }

        return result;
    }
}

TEST(correctness, div_randomized)
{
    for (size_t itn = 0; itn != number_of_iterations * number_of_multipliers; ++itn)
    {
        BigInteger divident = rand_big(10);
        BigInteger divisor = rand_big(6);
        BigInteger quotient = divident / divisor;
        BigInteger residue = divident % divisor;
        ASSERT_EQ(divident - quotient * divisor, residue);
        EXPECT_GE(residue, 0);
        EXPECT_LT(residue, divisor);
    }
}
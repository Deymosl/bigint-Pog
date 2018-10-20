//
//

#ifndef BIG_INT_BIGINT_H
#define BIG_INT_BIGINT_H

#include <algorithm>
#include <string>
#include <vector>
#include <climits>
#include <cmath>
#include <cstdint>

class BigInteger {
public:
    BigInteger();

    BigInteger(const BigInteger &source);

    BigInteger(int source);

    explicit BigInteger(std::string source);

    ~BigInteger();

    BigInteger operator+() const;

    BigInteger operator-() const;

    BigInteger operator~() const;

    BigInteger &operator=(const BigInteger &source);

    BigInteger &operator+=(const BigInteger &source);

    BigInteger &operator-=(const BigInteger &source);

    BigInteger &operator/=(const BigInteger &source);

    BigInteger &operator*=(const BigInteger &source);

    BigInteger &operator%=(const BigInteger &source);

    BigInteger &operator|=(const BigInteger &source);

    BigInteger &operator&=(const BigInteger &source);

    BigInteger &operator^=(const BigInteger &source);

    BigInteger &operator<<=(unsigned int rhs);

    BigInteger &operator>>=(unsigned int lhs);

    friend BigInteger operator+(BigInteger q, const BigInteger &w);

    friend BigInteger operator-(BigInteger q, const BigInteger &w);

    friend BigInteger operator*(BigInteger q, const BigInteger &w);

    friend BigInteger operator/(BigInteger q, const BigInteger &w);

    friend BigInteger operator%(BigInteger q, const BigInteger &w);


    friend BigInteger operator&(BigInteger q, const BigInteger &w);

    friend BigInteger operator|(BigInteger q, const BigInteger &w);

    friend BigInteger operator^(BigInteger q, const BigInteger &w);

    friend BigInteger operator<<(BigInteger q, unsigned int w);

    friend BigInteger operator>>(BigInteger q, unsigned int w);

    friend bool operator==(const BigInteger &q, const BigInteger &w);

    friend bool operator!=(const BigInteger &q, const BigInteger &w);

    friend bool operator>=(const BigInteger &q, const BigInteger &w);

    friend bool operator<=(const BigInteger &q, const BigInteger &w);

    friend bool operator>(const BigInteger &q, const BigInteger &w);

    friend bool operator<(const BigInteger &q, const BigInteger &w);

    friend std::string to_string(const BigInteger &q);

    friend std::ostream &operator<<(std::ostream &s, const BigInteger &q);

    size_t length() const;

    uint32_t operator[](size_t q) const;

private:
    std::vector<uint32_t> data;
    bool sign = false;

    std::vector<uint32_t> quotient(uint32_t q) const;

    BigInteger inverse() const;

    template<class Func>
    BigInteger &applyBitwiseOperation(const BigInteger &q, Func function);

    void shrink();
};

std::vector<uint32_t> &correct(std::vector<uint32_t> &q);

void divide(const std::vector<uint32_t> &q, const std::vector<uint32_t> &w, std::vector<uint32_t> &res);

bool isLess(const std::vector<uint32_t> &q, const std::vector<uint32_t> &w);

void add(const std::vector<uint32_t> &q, const std::vector<uint32_t> &w, std::vector<uint32_t> &res);

void sub(const std::vector<uint32_t> &q, const std::vector<uint32_t> &w, std::vector<uint32_t> &res);

void mul(const std::vector<uint32_t> &q, const std::vector<uint32_t> &w, std::vector<uint32_t> &res);

#endif //BIG_INT_BIGINT_H

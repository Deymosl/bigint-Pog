#include <iostream>
#include <cstdint>
#include <stdint-gcc.h>
#include "big_integer.h"

uint32_t get(const std::vector<uint32_t> &q, int w) {
    return (w < (int) q.size() ? q[w] : 0);
}

size_t BigInteger::length() const {
    return data.size();
}

void BigInteger::shrink() {
    correct(data);
    if (data.back() == 0)
        sign = false;
}

std::vector<uint32_t> &correct(std::vector<uint32_t> &q) {
    while (q.size() > 1 && q.back() == 0)
        q.pop_back();
    return q;
}

BigInteger::BigInteger() {
    sign = false;
    data.push_back(0);
}

BigInteger::BigInteger(int source) {
    sign = source < 0;
    if (sign) {
        source = -(source + 1);
        data.push_back(static_cast<uint32_t>(source) + 1);
    } else data.push_back(source);
}

BigInteger::BigInteger(const BigInteger &source) {
    sign = source.sign;
    data = source.data;
    shrink();
}

BigInteger::BigInteger(std::string source) {
    sign = false;
    BigInteger result(0);
    reverse(source.begin(), source.end());
    if (source.back() == '-') {
        source.pop_back();
        sign = true;
    }
    while (!source.empty()) {
        result *= 10;
        result += source.back() - '0';
        source.pop_back();
    }
    data = result.data;
    shrink();
}

BigInteger::~BigInteger() {
    sign = false;
    data.clear();
}

BigInteger BigInteger::operator-() const {
    if (*this == 0)
        return *this;
    BigInteger result = BigInteger(*this);
    result.sign = !sign;
    return result;
}

BigInteger BigInteger::operator+() const {
    return *this;
}

BigInteger BigInteger::operator~() const {
    BigInteger res(*this);
    return -(res + 1);
}

BigInteger &BigInteger::operator=(const BigInteger &q) {
    sign = q.sign;
    data = q.data;
    return *this;
}


BigInteger &BigInteger::operator+=(const BigInteger &q) {
    if (sign == q.sign) {
        add(data, q.data, data);
    } else {
        bool s = sign;
        const BigInteger &a1 = s ? q : *this;
        const BigInteger &a2 = !s ? q : *this;
        sign = isLess(a1.data, a2.data);
        sub(a1.data, a2.data, data);
    }
    return *this;
}

std::vector<uint32_t> BigInteger::quotient(uint32_t q) const {
    BigInteger res;
    res.data.resize(data.size());
    uint64_t carry = 0;
    for (int i = data.size() - 1; i >= 0; i--) {
        uint64_t temp = data[i] + (carry << 32);
        res.data[i] = temp / q;
        carry = temp % q;
    }
    correct(res.data);
    return res.data;
}

bool isLess(const std::vector<uint32_t> &q, const std::vector<uint32_t> &w) {
    if (q.size() < w.size())
        return true;
    if (q.size() > w.size())
        return false;
    return lexicographical_compare(q.rbegin(), q.rend(), w.rbegin(), w.rend());
}

void mul(const std::vector<uint32_t> &q, const std::vector<uint32_t> &w, std::vector<uint32_t> &res) {
    int n = q.size();
    res.resize(q.size() + w.size());
    res.back() = 0;
    uint32_t carry = 0;
    int j = 0;
    uint64_t cur = 0;
    for (int i = n - 1; i >= 0; i--) {
        uint32_t t = q[i];
        for (j = 0, carry = 0; j < (int) w.size() || carry > 0; j++) {
            cur = (uint64_t) res[i + j] * (j != 0) + (uint64_t) t * (uint64_t) get(w, j);
            res[i + j] = cur + carry;
            carry = (cur + carry) >> 32;
        }
    }
    correct(res);
}

void add(const std::vector<uint32_t> &q, const std::vector<uint32_t> &w, std::vector<uint32_t> &res) {
    res.resize(std::max(q.size(), w.size()) + 1);
    uint32_t carry = 0;
    for (unsigned int i = 0; i < std::max(q.size(), w.size()); i++) {
        uint32_t qi = get(q, i);
        uint32_t wi = get(w, i);
        res[i] = qi + wi + carry;
        carry = (UINT32_MAX - qi < wi + carry); //hm
    }
    if (carry)
        res.back() = 1;
    correct(res);
}

void sub(const std::vector<uint32_t> &q, const std::vector<uint32_t> &w, std::vector<uint32_t> &res) {
    bool less = isLess(q, w);
    const std::vector<uint32_t> &a = (less ? w : q);
    const std::vector<uint32_t> &b = (less ? q : w);
    res.resize(std::max(a.size(), b.size()));
    for (int i = a.size() - 1; i >= 0; i--) {
        uint32_t ai = get(a, i);
        uint32_t bi = get(b, i);
        res[i] = ai - bi + (ai < bi) * UINT32_MAX + (ai < bi);
        if (ai < bi) {
            int j = i + 1;
            while (res[j] == 0) {
                res[j]--;
                j++;
            }
            res[j]--;
        }
    }
    correct(res);
}

void divide(const std::vector<uint32_t> &q, const std::vector<uint32_t> &w, std::vector<uint32_t> &res) {
    const unsigned int fac = (unsigned int) ((1ll << 32) / ((unsigned long long) (w.back() + 1)));
    res.resize(0);
    std::vector<uint32_t> qFactorized;
    mul(q, std::vector<uint32_t>(1, fac), qFactorized);
    std::vector<uint32_t> wFactorized;
    mul(w, std::vector<uint32_t>(1, fac), wFactorized);
    int n = qFactorized.size();
    std::vector<uint32_t> temp;
    for (int i = n - 1; i >= 0; i--) {
        temp.insert(temp.begin(), qFactorized[i]);
        correct(temp);
        if (isLess(temp, wFactorized)) {
            if (temp.back() == 0)
                res.push_back(0);
            continue;
        }
        uint64_t dividend = uint64_t(temp.back());
        if (temp.size() != wFactorized.size()) {
            dividend <<= 32;
            dividend += temp[temp.size() - 2];
        }
        uint64_t divider = uint64_t(wFactorized.back());
        uint64_t quot = dividend / divider;
        quot = std::min(quot, ((uint64_t) 1 << 32) - 1);
        std::vector<uint32_t> buffer;
        mul(wFactorized, std::vector<uint32_t>(1, quot), buffer);
        if (isLess(temp, buffer)) {
            sub(buffer, wFactorized, buffer);
            quot--;
        }
        if (isLess(temp, buffer)) {
            uint64_t l = 0, r = quot;
            while (r - l > 1) {
                uint64_t mid = (r + l) / 2;
                std::vector<uint32_t> tt;
                mul(wFactorized, std::vector<uint32_t>(1, mid), tt);
                if (isLess(temp, tt)) {
                    r = mid;
                } else {
                    l = mid;
                }
            }
            quot = l;
            mul(wFactorized, std::vector<uint32_t>(1, quot), buffer);
        }
        sub(temp, buffer, temp);
        res.push_back(quot);
    }
    reverse(res.begin(), res.end());
}

BigInteger &BigInteger::operator/=(const BigInteger &q) {
    if (q.data.size() == 1) {
        data = quotient(q.data[0]);
        sign ^= q.sign;
    } else if (q != 0) {
        if (isLess(data, q.data)) {
            data = std::vector<uint32_t>(1, 0);
            sign = false;
        } else {
            sign ^= q.sign;
            std::vector<uint32_t> res;
            divide(data, q.data, res);
            data = res;
        }
    }
    correct(data);
    if (data.back() == 0)
        sign = false;
    return *this;
}

BigInteger &BigInteger::operator*=(const BigInteger &q) {
    sign ^= q.sign;
    mul(data, q.data, data);
    if (data.back() == 0)
        sign = false;
    return *this;
}

BigInteger &BigInteger::operator%=(const BigInteger &q) {
    BigInteger temp(*this);
    temp /= q;
    temp *= q;
    *this -= temp;
    return *this;
}

bool operator==(const BigInteger &q, const BigInteger &w) {
    if (q.sign != w.sign || q.length() != w.length())
        return false;
    for (unsigned int i = q.length(); i > 0; i--) {
        if (q[i - 1] != w[i - 1])
            return false;
    }
    return true;
}

bool operator!=(const BigInteger &q, const BigInteger &w) {
    return !(q == w);
}

bool operator<(const BigInteger &q, const BigInteger &w) {
    if (q.sign != w.sign)
        return (q.sign > w.sign);
    if (q.length() != w.length()) {
        return (q.length() < w.length()) ^ q.sign;
    }
    return isLess(q.sign ? w.data : q.data, q.sign ? q.data : w.data);
}

bool operator>(const BigInteger &q, const BigInteger &w) {
    return w < q;
}

bool operator>=(const BigInteger &q, const BigInteger &w) {
    return !(q < w);
}

bool operator<=(const BigInteger &q, const BigInteger &w) {
    return !(q > w);
}

template<class Func>
BigInteger &BigInteger::applyBitwiseOperation(const BigInteger &q, Func function) {
    data.resize(std::max(data.size(), q.data.size()));
    if (sign)
        data = inverse().data;
    const BigInteger &right = (q.sign ? q.inverse() : q);
    for (size_t i = 0; i < length(); i++)
        data[i] = function(data[i], right[i]);
    correct(data);
    sign = function(sign, q.sign);
    if (sign)
        data = inverse().data;
    return *this;

}

BigInteger BigInteger::inverse() const {
    BigInteger res(*this);
    for (auto &i : res.data) i = UINT32_MAX - i;
    res.sign ^= true;
    return res + 1;
}

BigInteger &BigInteger::operator&=(const BigInteger &q) {
    return applyBitwiseOperation(q, std::bit_and<uint32_t>());
}

BigInteger operator&(BigInteger q, const BigInteger &w) {
    return q &= w;
}

BigInteger &BigInteger::operator|=(const BigInteger &q) {
    return applyBitwiseOperation(q, std::bit_or<uint32_t>());
}

BigInteger operator|(BigInteger q, const BigInteger &w) {
    return q |= w;
}

BigInteger &BigInteger::operator^=(const BigInteger &q) {
    return applyBitwiseOperation(q, std::bit_xor<uint32_t>());
}

BigInteger operator^(BigInteger q, const BigInteger &w) {
    return q ^= w;
}

int getDeg(uint32_t q) {
    int res = 0;
    while (q) {
        res++;
        q /= 2;
    }
    return res;
}

BigInteger &BigInteger::operator<<=(unsigned int q) {
    int maxDeg = getDeg(data.back());
    int rem = q % 32;
    int carry = (maxDeg + rem > 32);
    std::vector<uint32_t> res;
    res.resize(data.size() + q / 32 + carry);
    copy(data.begin(), data.end(), res.begin());
    for (int i = (int) res.size() - 1; i >= (int) res.size() - (int) data.size(); i--) {
        res[i] <<= rem;
        if (i != (int) res.size() - (int) data.size())
            res[i] += res[i - 1] >> (32 - maxDeg);
    }
    data = res;
    return *this;
}

BigInteger operator<<(BigInteger q, unsigned int w) {
    return (q <<= w);
}

BigInteger &BigInteger::operator>>=(unsigned int q) {
    std::vector<uint32_t> res;
    res.resize(std::max(0, (int) ((int) data.size() - q / 32)));
    copy(data.begin() + q / 32, data.end(), res.begin());
    if (res.empty())
        res.push_back(0);
    else {
        uint32_t next = 0;
        for (int i = 0; i < q; i++) next *= 2, next += 1;
        for (int i = 0; i < res.size(); i++) {
            res[i] >>= q;
            if (i != res.size() - 1)
                res[i] += (res[i + 1] & next) << (sizeof(uint32_t) * 8 - q);
        }
    }
    correct(res);
    data = res;
    if (sign)
        add(data, std::vector<uint32_t>(1, 1), data);
    shrink();
    return *this;
}

BigInteger operator>>(BigInteger q, unsigned int w) {
    BigInteger f = (q >>= w);
    return f;
}

uint32_t BigInteger::operator[](size_t pos) const {
    if (pos < length())
        return data[pos];
    return 0;
}

BigInteger &BigInteger::operator-=(const BigInteger &q) {
    return *this += -q;
}

BigInteger operator+(BigInteger q, const BigInteger &w) {
    return q += w;
}

BigInteger operator-(BigInteger q, const BigInteger &w) {
    return q -= w;
}

BigInteger operator*(BigInteger q, const BigInteger &w) {
    return q *= w;
}

BigInteger operator/(BigInteger q, const BigInteger &w) {
    return q /= w;
}

BigInteger operator%(BigInteger q, const BigInteger &w) {
    return q %= w;
}

std::string to_string(BigInteger const &q) {
    if (q == 0)
        return "0";
    std::string res = "";
    BigInteger w = q;
    while (w != 0) {
        res += '0' + (w % 10).data.back();
        w /= 10;
    }
    if (q.sign)
        res.push_back('-');
    reverse(res.begin(), res.end());
    return res;
}

std::ostream &operator<<(std::ostream &os, const BigInteger &q) {
    std::string res = to_string(q);
    os << res;
    return os;
}

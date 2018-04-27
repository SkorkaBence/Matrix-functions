#pragma once

#include <iostream>
#include <string>

namespace sbl {

    template<typename T>
    class basic_rational {
    public:
        class divisionByZeroException {};
    private:

        T numerator;
        T denominator;

        T abs(T a) {
            if (a >= 0) {
                return a;
            } else {
                return -1 * a;
            }
        }

        T GreatestCommonDivisor(T a, T b) {
            a = abs(a);
            b = abs(b);
            if (a == 0 || b == 0 || a == 1 || b == 1) {
                return 1;
            }
            while (a != b) {
                if (a > b) {
                    a -= b;
                } else {
                    b -= a;
                }
            }
            return a;
        }

        void Simplification() {
            if (denominator < 0) {
                numerator *= -1;
                denominator *= -1;
            }

            T div = GreatestCommonDivisor(numerator, denominator);
            numerator /= div;
            denominator /= div;
            if (denominator == 0) {
                throw divisionByZeroException();
            }
            if (numerator == 0) {
                denominator = 1;
            }
        }

    public:

        basic_rational() : numerator(0), denominator(1) {}
        basic_rational(T n) : numerator(n), denominator(1) {}
        basic_rational(T n, T d) : numerator(n), denominator(d) {
            Simplification();
        }
        basic_rational(const basic_rational& r) : numerator(r.numerator), denominator(r.denominator) {}

        ~basic_rational() {}

        void prettyPrint(std::ostream& out) {
            if (denominator == 1) {
                out << numerator;
            } else {
                out << numerator << "/" << denominator;
            }
        }

        friend std::ostream& operator << (std::ostream& out, const basic_rational& dt) {
            out << dt.numerator << " " << dt.denominator;
            return out;
        }

        friend std::istream& operator >> (std::istream& in, basic_rational& dt) {
            in >> dt.numerator;
            in >> dt.denominator;
            return in;
        }

        basic_rational& operator = (const basic_rational& ot) {
            numerator = ot.numerator;
            denominator = ot.denominator;
            return *this;
        }

        basic_rational& operator += (const basic_rational& ot) {
            numerator = numerator * ot.denominator + ot.numerator * denominator;
            denominator *= ot.denominator;
            Simplification();
            return *this;
        }

        basic_rational operator + (const basic_rational& ot) const {
            basic_rational r(numerator, denominator);
            r += ot;
            return r;
        }

        basic_rational& operator -= (const basic_rational& ot) {
            numerator = numerator * ot.denominator - ot.numerator * denominator;
            denominator *= ot.denominator;
            Simplification();
            return *this;
        }

        basic_rational operator - (const basic_rational& ot) const {
            basic_rational r(numerator, denominator);
            r -= ot;
            return r;
        }

        basic_rational& operator *= (const basic_rational& ot) {
            numerator *= ot.numerator;
            denominator *= ot.denominator;
            Simplification();
            return *this;
        }

        basic_rational operator * (const basic_rational& ot) const {
            basic_rational r(numerator, denominator);
            r *= ot;
            return r;
        }

        basic_rational& operator /= (const basic_rational& ot) {
            numerator *= ot.denominator;
            denominator *= ot.numerator;
            Simplification();
            return *this;
        }

        basic_rational operator / (const basic_rational& ot) const {
            basic_rational r(numerator, denominator);
            r /= ot;
            return r;
        }

        basic_rational& operator = (T ot) {
            numerator = ot;
            denominator = 1;
            return *this;
        }

        basic_rational& operator += (T ot) {
            numerator += ot * denominator;
            Simplification();
            return *this;
        }

        basic_rational operator + (T ot) const {
            basic_rational r(numerator, denominator);
            r += ot;
            return r;
        }

        basic_rational& operator -= (T ot) {
            numerator -= ot * denominator;
            denominator *= ot.denominator;
            Simplification();
            return *this;
        }

        basic_rational operator - (T ot) const {
            basic_rational r(numerator, denominator);
            r -= ot;
            return r;
        }

        basic_rational& operator *= (T ot) {
            numerator *= ot;
            Simplification();
            return *this;
        }

        basic_rational operator * (T ot) const {
            basic_rational r(numerator, denominator);
            r *= ot;
            return r;
        }

        basic_rational& operator /= (T ot) {
            denominator *= ot;
            Simplification();
            return *this;
        }

        basic_rational operator / (T ot) const {
            basic_rational r(numerator, denominator);
            r /= ot;
            return r;
        }

        // rac <-> rac ----------------------------------------------------------------

        bool operator < (const basic_rational& ot) const {
            return (numerator * ot.denominator) < (ot.numerator * denominator);
        }

        bool operator > (const basic_rational& ot) const {
            return ot < *this;
        }

        bool operator <= (const basic_rational& ot) const {
            return !(*this > ot);
        }

        bool operator >= (const basic_rational& ot) const {
            return !(*this < ot);
        }

        bool operator == (const basic_rational& ot) const {
            return (*this <= ot) && (*this >= ot);
        }

        bool operator != (const basic_rational& ot) const {
            return (*this < ot) || (*this > ot);
        }

        // rac <-> T ----------------------------------------------------------------

        bool operator < (T ot) const {
            basic_rational r(ot);
            return *this < r;
        }

        bool operator > (T ot) const {
            basic_rational r(ot);
            return r < *this;
        }

        bool operator <= (T ot) const {
            return !(*this > ot);
        }

        bool operator >= (T ot) const {
            return !(*this < ot);
        }

        bool operator == (T ot) const {
            return (*this <= ot) && (*this >= ot);
        }

        bool operator != (T ot) const {
            return (*this < ot) || (*this > ot);
        }

        // T <-> rac ----------------------------------------------------------------

        friend bool operator < (T ot, const basic_rational& th) {
            return th > ot;
        }

        friend bool operator > (T ot, const basic_rational& th) {
            return ot < th;
        }

        friend bool operator <= (T ot, const basic_rational& th) {
            return !(th > ot);
        }

        friend bool operator >= (T ot, const basic_rational& th) {
            return !(th < ot);
        }

        friend bool operator == (T ot, const basic_rational& th) {
            return (th <= ot) && (th >= ot);
        }

        friend bool operator != (T ot, const basic_rational& th) {
            return (th < ot) || (th > ot);
        }

        T val() const {
            return numerator / denominator;
        }

        double double_val() const {
            return numerator / denominator;
        }

        operator double() const {
            return double_val();
        }

    };

    typedef basic_rational<int> rational;

}

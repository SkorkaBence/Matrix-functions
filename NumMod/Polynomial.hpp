#pragma once

#include <iostream>
#include <vector>

#ifndef NEGATIVE_INFINITY
#define NEGATIVE_INFINITY -2147483647
#endif // !NEGATIVE_INFINITY

namespace sbl {

    template<typename T>
    class basic_polynomial {
    public:
        class divisionByZeroException {};
    private:
        T* numarr = nullptr;
        unsigned int size = 0;

        void deleteArray() {
            if (numarr != nullptr) {
                delete[] numarr;
                numarr = nullptr;
            }
        }

        void resize(unsigned int newsize) {
            if (size == newsize) {

            } else if (newsize == 0) {
                deleteArray();
            } else {
                int keep = newsize > size ? size : newsize;
                T* newarr = new T[newsize];
                for (int i = 0; i < keep; i++) {
                    newarr[i] = numarr[i];
                }
                delete[] numarr;
                numarr = newarr;
            }
            size = newsize;
        }

        void push_back(T val) {
            resize(size + 1);
            numarr[size - 1] = val;
        }

        void removeLeadingZeroes() {
            int newsize = size;
            while (newsize > 0 && numarr[newsize - 1] == 0) {
                newsize--;
            }
            resize(newsize);
        }

    public:
        basic_polynomial() {}
        basic_polynomial(const std::vector<T>& nums) {
            resize(nums.size());
            for (int i = 0; i < nums.size(); i++) {
                numarr[i] = nums[i];
            }
            removeLeadingZeroes();
        }
        basic_polynomial(const basic_polynomial& copy) {
            resize(copy.size);
            for (unsigned i = 0; i < copy.size; i++) {
                numarr[i] = copy.numarr[i];
            }
            removeLeadingZeroes();
        }

        friend std::ostream& operator << (std::ostream& out, const basic_polynomial& dt) {
            out << "P(x) =";
            if (dt.size == 0) {
                out << " 0";
            } else {
                for (int i = dt.size - 1; i >= 0; i--) {
                    out << " (" << dt.numarr[i] << ")";
                    if (i > 0) {
                        out << "x^" << i << " +";
                    }
                }
            }
            return out;
        }

        T operator () (T inp) const {
            T r = 0;
            for (unsigned int i = 0; i < size; i++) {
                T x = numarr[i];
                for (int j = 0; j < i; j++) {
                    x *= inp;
                }
                r += x;
            }
            return r;
        }

        int deg() const {
            if (size == 0) {
                return NEGATIVE_INFINITY;
            }
            int i = size - 1;
            while (i >= 0 && numarr[i] == 0) {
                i--;
            }
            if (i <= 0) {
                return NEGATIVE_INFINITY;
            }
            return i;
        }

        T& operator [] (unsigned int pos) {
            return numarr[pos];
        }

        const T& operator [] (unsigned int pos) const {
            return numarr[pos];
        }

        void setValue(unsigned int pos, T val) {
            if (pos >= size) {
                int oldsize = size;
                resize(pos + 1);
                for (unsigned int i = oldsize; i < size; i++) {
                    numarr[i] = 0;
                }
            }
            numarr[pos] = val;
        }

        T getValue(unsigned int pos) const {
            if (pos >= 0 && size > pos) {
                return numarr[pos];
            } else {
                return 0;
            }
        }

        bool isNull() const {
            if (size == 0) {
                return true;
            }
            unsigned int i = size - 1;
            while (i >= 0 && numarr[i] == 0) {
                i--;
            }
            return (i < 0);
        }

        // Pol <-> T ----------------------------------------------------------------

        basic_polynomial& operator += (T val) {
            setValue(0, getValue(0) + val);
            removeLeadingZeroes();
            return *this;
        }

        basic_polynomial operator + (T val) const {
            basic_polynomial p(*this);
            p += val;
            return p;
        }

        basic_polynomial& operator -= (T val) {
            setValue(0, getValue(0) - val);
            removeLeadingZeroes();
            return *this;
        }

        basic_polynomial operator - (T val) const {
            basic_polynomial p(*this);
            p -= val;
            return p;
        }

        basic_polynomial& operator *= (T val) {
            for (unsigned int i = 0; i < size; i++) {
                numarr[i] *= val;
            }
            removeLeadingZeroes();
            return *this;
        }

        basic_polynomial operator * (T val) const {
            basic_polynomial p(*this);
            p *= val;
            return p;
        }

        basic_polynomial& operator /= (T val) {
            for (unsigned int i = 0; i < size; i++) {
                numarr[i] /= val;
            }
            removeLeadingZeroes();
            return *this;
        }

        basic_polynomial operator / (T val) const {
            basic_polynomial p(*this);
            p /= val;
            return p;
        }

        bool operator == (T val) {
            if (size > 0) {
                bool ok = (numarr[0] == val);
                for (unsigned int i = 1; ok && i < size; i++) {
                    ok = ok & (numarr[i] == 0);
                }
                return ok;
            } else {
                return false;
            }
        }

        bool operator != (T val) {
            return !(*this == val);
        }

        // T <-> Pol ----------------------------------------------------------------

        friend basic_polynomial operator + (T val, const basic_polynomial& pol) {
            basic_polynomial p(pol);
            p += val;
            p.removeLeadingZeroes();
            return p;
        }

        friend basic_polynomial operator - (T val, const basic_polynomial& pol) {
            basic_polynomial p(pol);
            for (unsigned int i = 0; i < p.size; i++) {
                p.numarr[i] *= -1;
            }
            p += val;
            p.removeLeadingZeroes();
            return p;
        }

        friend basic_polynomial operator * (T val, const basic_polynomial& pol) {
            basic_polynomial p(pol);
            p *= val;
            p.removeLeadingZeroes();
            return p;
        }

        friend basic_polynomial operator / (T val, const basic_polynomial& pol) {
            basic_polynomial p(pol);
            for (unsigned int i = 0; i < p.size; i++) {
                if (p.numarr[i] != 0) {
                    p.numarr[i] = 1 / p.numarr[i];
                }
            }
            p *= val;
            p.removeLeadingZeroes();
            return p;
        }

        friend bool operator == (T val, const basic_polynomial& pol) {
            return pol == val;
        }

        friend bool operator != (T val, const basic_polynomial& pol) {
            return pol != val;
        }

        // Pol <-> Pol ----------------------------------------------------------------

        basic_polynomial& operator = (const basic_polynomial& pol) {
            resize(pol.size);
            for (unsigned int i = 0; i < size; i++) {
                numarr[i] = pol.numarr[i];
            }
            removeLeadingZeroes();
            return *this;
        }

        basic_polynomial& operator += (const basic_polynomial& pol) {
            int d1 = deg();
            int d2 = pol.deg();
            int mx = (d1 > d2 ? d1 : d2);
            for (int i = 0; i <= mx; i++) {
                setValue(i, getValue(i) + pol.getValue(i));
            }
            return *this;
        }

        basic_polynomial operator + (const basic_polynomial& pol) const {
            basic_polynomial p(*this);
            p += pol;
            return p;
        }

        basic_polynomial& operator -= (const basic_polynomial& pol) {
            int min = size > pol.size ? pol.size : size;
            int max = size < pol.size ? pol.size : size;
            resize(max);
            for (unsigned int i = 0; i < min; i++) {
                setValue(i, getValue(i) - pol.getValue(i));
            }
            removeLeadingZeroes();
            return *this;
        }

        basic_polynomial operator - (const basic_polynomial& pol) const {
            basic_polynomial p(*this);
            p -= pol;
            return p;
        }

        basic_polynomial& operator *= (const basic_polynomial& pol) {
            *this = *this * pol;
            return *this;
        }

        basic_polynomial operator * (const basic_polynomial& pol) const {
            basic_polynomial p;
            int s1 = deg();
            int s2 = pol.deg();
            for (int i = 0; i < ((s1 + 1) + (s2 + 1)); i++) {
                T g = 0;
                for (int k = 0; k <= i; k++) {
                    g += getValue(k) * pol.getValue(i - k);
                }
                p.setValue(i, g);
            }
            return p;
        }

        basic_polynomial& operator /= (const basic_polynomial& pol) {
            *this = *this / pol;
            return *this;
        }

        basic_polynomial operator / (const basic_polynomial& pol) const {
            basic_polynomial q;
            basic_polynomial r;
            divideWithResidue(pol, q, r);
            return q;
        }

        void divideWithResidue(const basic_polynomial& pol, basic_polynomial& sol, basic_polynomial& res) const {
            if (pol.isNull()) {
                throw divisionByZeroException();
            }

            basic_polynomial q;
            basic_polynomial r(*this);

            while (!r.isNull() && r.deg() >= pol.deg()) {
                T leading_v = r.getValue(r.deg()) / pol.getValue(pol.deg());
                int leadign_p = r.deg() - pol.deg();
                basic_polynomial t;
                t.setValue(leadign_p, leading_v);
                q += t;
                r -= t * pol;
            }

            sol = q;
            res = r;
        }

        basic_polynomial derivate() {
            basic_polynomial p;
            for (int i = size - 1; i > 0; i--) {
                p.setValue(i - 1, getValue(i) * i);
            }
            return p;
        }

        basic_polynomial derivate(unsigned int n) {
            basic_polynomial p(*this);
            for (unsigned int i = 0; i < n; i++) {
                p = p.derivate();
            }
            return p;
        }

    };

    typedef basic_polynomial<int> polynomial;

}

#pragma once

#include <iostream>

namespace sbl {

    template<typename T>
    class basic_complex {
    private:

        T re;
        T im;

    public:

        basic_complex() : re(0), im(0) {}
        basic_complex(T r) : re(r), im(0) {}
        basic_complex(T r, T i) : re(r), im(i) {}
        basic_complex(const basic_complex& copy) : re(copy.re), im(copy.im) {}

        ~basic_complex() {}

        friend std::ostream& operator << (std::ostream& ou, const basic_complex& cn) {
            ou << "(" << cn.re << "|" << cn.im << ")";
            return ou;
        }

        basic_complex& operator = (const basic_complex& oc) {
            re = oc.re;
            im = oc.im;
            return *this;
        }

        basic_complex& operator += (const basic_complex& oc) {
            re += oc.re;
            im += oc.im;
            return *this;
        }

        basic_complex operator + (const basic_complex& oc) const {
            basic_complex r(re, im);
            r += oc;
            return r;
        }

        basic_complex& operator -= (const basic_complex& oc) {
            re -= oc.re;
            im -= oc.im;
            return *this;
        }

        basic_complex operator - (const basic_complex& oc) const {
            basic_complex r(re, im);
            r -= oc;
            return r;
        }

        basic_complex& operator *= (const basic_complex& oc) {
            T newre = re * oc.re - im * oc.im;
            T newim = re * oc.im + im * oc.re;
            re = newre;
            im = newim;
            return *this;
        }

        basic_complex operator * (const basic_complex& oc) const {
            basic_complex r(re, im);
            r *= oc;
            return r;
        }

        basic_complex& operator /= (const basic_complex& oc) {
            T d = oc.re * oc.re + oc.im * oc.im;
            T newre = re * oc.re + im * oc.im;
            newre /= d;
            T newim = im * oc.re - re * oc.im;
            newim /= d;
            re = newre;
            im = newim;
            return *this;
        }

        basic_complex operator / (const basic_complex& oc) const {
            basic_complex r(re, im);
            r /= oc;
            return r;
        }

        basic_complex& operator = (T oc) {
            re = oc;
            im = 0;
            return *this;
        }

        basic_complex& operator += (T oc) {
            re += oc;
            return *this;
        }

        basic_complex operator + (T oc) const {
            basic_complex r(re, im);
            r += oc;
            return r;
        }

        basic_complex& operator -= (T oc) {
            re -= oc;
            return *this;
        }

        basic_complex operator - (T oc) const {
            basic_complex r(re, im);
            r -= oc;
            return r;
        }

        basic_complex& operator *= (T oc) {
            re *= oc;
            im *= oc;
            return *this;
        }

        basic_complex operator * (T oc) const {
            basic_complex r(re, im);
            r *= oc;
            return r;
        }

        basic_complex& operator /= (T oc) {
            re /= oc;
            im /= oc;
            return *this;
        }

        basic_complex operator / (T oc) const {
            basic_complex r(re, im);
            r /= oc;
            return r;
        }

        // Comp <-> Comp ----------------------------------------------------------------

        bool operator == (const basic_complex& oc) const {
            return (re == oc.re) && (im == oc.im);
        }

        bool operator != (const basic_complex& oc) const {
            return !(*this == oc);
        }

        // Comp <-> T ----------------------------------------------------------------

        bool operator == (T oc) const {
            return (re == oc) && (im == 0);
        }

        bool operator != (T oc) const {
            return !(*this == oc);
        }

        // T <-> Comp ----------------------------------------------------------------

        friend bool operator == (T oc, const basic_complex& th) {
            return (th.re == oc) && (th.im == 0);
        }

        friend bool operator != (T oc, const basic_complex& th) {
            return !(th == oc);
        }

    };

    typedef basic_complex<int> complex;

}
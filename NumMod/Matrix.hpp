#pragma once

#include <iostream>
#include <random>

namespace sbl {

    template<typename T = int>
    class Matrix {
    public:
        class SquareMatrixRequired {};
        class InvalidSizeException {};
        class ProcessStuckException {};

        struct LU {
            Matrix L;
            Matrix U;
        };

        struct LDU {
            Matrix L;
            Matrix D;
            Matrix U;
        };

        Matrix(const Matrix& copy) {
            resize(copy.Width, copy.Height, 0);
            for (unsigned i = 0; i < real_size; i++) {
                arr[i] = copy.arr[i];
            }
        }

        Matrix(unsigned int w, unsigned int h, T def) {
            resize(w, h, def);
        }

        Matrix(unsigned int w, unsigned int h) : Matrix(w, h, 0) {}

        Matrix() : Matrix(0, 0) {}

        ~Matrix() {
            delete[] arr;
        }

        const T& operator () (unsigned i, unsigned j) const {
            return arr[i * Width + j];
        }

        T& operator () (unsigned i, unsigned j) {
            return arr[i * Width + j];
        }

        unsigned getWidth() const {
            return Width;
        }

        unsigned getHeight() const {
            return Height;
        }

        friend std::ostream& operator << (std::ostream& out, const Matrix& m) {
            out << m.Width << " " << m.Height;
            for (unsigned i = 0; i < m.Height; i++) {
                out << std::endl;
                for (unsigned j = 0; j < m.Width; j++) {
                    out << m(i, j) << '\t';
                }
            }
            return out;
        }

        friend std::istream& operator >> (std::istream& in, Matrix& m) {
            unsigned w;
            unsigned h;
            in >> w;
            in >> h;
            m.resize(w, h, 0);
            for (unsigned i = 0; i < h; i++) {
                for (unsigned j = 0; j < w; j++) {
                    in >> m(i, j);
                }
            }
            return in;
        }

        Matrix& operator = (Matrix o) {
            resize(o.Width, o.Height, 0);
            for (int i = 0; i < real_size; i++) {
                arr[i] = o.arr[i];
            }
            return *this;
        }

        Matrix& operator += (const Matrix& o) {
            if (o.Width != Width || o.Height != Height) {
                throw InvalidSizeException();
            }
            for (int i = 0; i < real_size; i++) {
                arr[i] += o.arr[i];
            }
            return *this;
        }

        Matrix operator + (const Matrix& o) const {
            Matrix newm(*this);
            newm += o;
            return newm;
        }

        Matrix& operator -= (const Matrix& o) {
            if (o.Width != Width || o.Height != Height) {
                throw InvalidSizeException();
            }
            for (int i = 0; i < real_size; i++) {
                arr[i] -= o.arr[i];
            }
            return *this;
        }

        Matrix operator - (const Matrix& o) const {
            Matrix newm(*this);
            newm -= o;
            return newm;
        }

        Matrix& operator *= (const Matrix& o) {
            *this = *this * o;
            return *this;
        }

        Matrix operator * (const Matrix& o) const {
            if (Width != o.Height) {
                throw InvalidSizeException();
            }

            Matrix a(o.Width, Height, 0);

            for (unsigned i = 0; i < a.Height; i++) {
                for (unsigned j = 0; j < a.Width; j++) {
                    T val = 0;
                    for (unsigned k = 0; k < Width; k++) {
                        val += (*this)(i, k) * o(k, j);
                    }
                    a(i, j) = val;
                }
            }

            return a;
        }

        Matrix operator ~ () {
            Matrix a(Height, Width);
            for (unsigned i = 0; i < Height; i++) {
                for (unsigned j = 0; j < Width; j++) {
                    a(j, i) = (*this)(i, j);
                }
            }
            return a;
        }

        void cleanResize(unsigned w, unsigned h) {
            resize(w, h, 0);
        }

        void randomize() {
            std::default_random_engine generator;
            std::uniform_int_distribution<int> distribution(1, 10);

            for (unsigned i = 0; i < real_size; i++) {
                arr[i] = distribution(generator);
            }
        }

    private:

        void resize(unsigned w, unsigned h, T def) {
            if (arr != nullptr) {
                delete[] arr;
            }
            Width = w;
            Height = h;
            real_size = Width * Height;
            arr = new T[real_size];
            for (unsigned i = 0; i < real_size; i++) {
                arr[i] = def;
            }
        }

        unsigned int Width;
        unsigned int Height;
        unsigned int real_size;

        T * arr = nullptr;
    };

}
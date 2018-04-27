#pragma once

#include <cmath>
#include "Matrix.hpp"
#include "Polynomial.hpp"

namespace sbl {

    template<typename T>
    void GaussianElimination(Matrix<T>& m) {
        GaussianElimination(m, 0);
    }

    template<typename T>
    void GaussianElimination(Matrix<T>& m, unsigned step) {
        unsigned w = m.getWidth();
        unsigned h = m.getHeight();
        if (w < h) {
            throw Matrix<T>::InvalidSizeException();
        }
        if (step >= h) {
            GaussianSubstitution(m);
            return;
        }

        T mainItem = m(step, step);

        if (mainItem == 0) {
            throw Matrix<T>::ProcessStuckException();
        }

        for (unsigned i = step + 1; i < h; ++i) {
            T referenceItem = m(i, step);
            T multiplier = referenceItem / mainItem;
            for (unsigned j = step; j < w; j++) {
                m(i, j) = m(i, j) - (m(step, j) * multiplier);
            }
        }

        GaussianElimination(m, step + 1);
    }

    template<typename T>
    void GaussianSubstitution(Matrix<T>& m) {
        GaussianSubstitution(m, 0);
    }

    template<typename T>
    void GaussianSubstitution(Matrix<T>& m, unsigned step) {
        unsigned w = m.getWidth();
        unsigned h = m.getHeight();

        int step_level = h - step - 1;
        if (step_level < 0) {
            return;
        }

        T mainItem = m(step_level, step_level);

        for (unsigned j = step_level; j < w; ++j) {
            m(step_level, j) /= mainItem;
        }

        for (unsigned i = 0; i < step_level; ++i) {
            T referenceItem = m(i, step_level);
            for (unsigned j = step_level; j < w; ++j) {
                m(i, j) = m(i, j) - (referenceItem * m(step_level, j));
            }
        }

        GaussianSubstitution(m, step + 1);
    }

    template<typename T>
    void Inverse(Matrix<T>& m) {
        unsigned Width = m.getWidth();
        unsigned Height = m.getHeight();

        if (Width != Height) {
            throw Matrix<T>::SquareMatrixRequired();
        }

        Matrix<T> inv(Width * 2, Height);

        for (unsigned i = 0; i < Height; ++i) {
            for (unsigned j = 0; j < Width; ++j) {
                inv(i, j) = m(i, j);
                inv(i, j + Width) = (i == j ? 1 : 0);
            }
        }

        GaussianElimination(inv);

        for (unsigned i = 0; i < Height; ++i) {
            for (unsigned j = 0; j < Width; ++j) {
                m(i, j) = inv(i, j + Width);
            }
        }
    }

    template<typename T>
    typename Matrix<T>::LDU LDUDecomposition(const Matrix<T>& original) {
        unsigned Width = original.getWidth();
        unsigned Height = original.getHeight();

        if (Width != Height) {
            throw Matrix<T>::SquareMatrixRequired();
        }

        Matrix<T> m(original);

        for (unsigned step = 0; step < Width - 1; ++step) {
            T divider = m(step, step);
            for (unsigned i = step + 1; i < Width; ++i) {
                m(i, step) /= divider;
            }
            for (unsigned i = step + 1; i < Width; ++i) {
                for (unsigned j = step + 1; j < Width; ++j) {
                    m(i, j) -= m(step, j) * m(i, step);
                }
            }
        }

        typename Matrix<T>::LDU res;
        res.L.cleanResize(Width, Height);
        res.D.cleanResize(Width, Height);
        res.U.cleanResize(Width, Height);

        for (unsigned i = 0; i < Height; ++i) {
            for (unsigned j = 0; j < Width; ++j) {
                if (i == j) {
                    res.L(i, j) = 1;
                    res.D(i, j) = m(i, j);
                    res.U(i, j) = 1;
                } else if (i > j) {
                    res.L(i, j) = m(i, j);
                } else if (i < j) {
                    res.U(i, j) = (m(i, j) / m(i, i));
                }
            }
        }

        return res;
    }

    template<typename T>
    T det(const Matrix<T> m) {
        typename Matrix<T>::LDU ldu = LDUDecomposition(m);

        T det = 1;
        for (unsigned i = 0; i < m.getHeight(); ++i) {
            det *= ldu.D(i, i);
        }

        return det;
    }

    template<typename T>
    Matrix<T> HouseholderMatrix(const Matrix<T>& m) {
        if (m.getWidth() != 1) {
            throw Matrix<T>::InvalidSizeException();
        }

        unsigned h = m.getHeight();

        Matrix<T> identity(h, h);
        for (unsigned i = 0; i < h; i++) {
            identity(i, i) = 1;
        }

        Matrix<T> mT = ~m;

        return identity - 2 * m * mT;
    }

    template<typename T>
    T norm1(const Matrix<T>& m) {
        T norm = 0;

        for (unsigned j = 0; j < m.getWidth(); ++j) {
            T sum = 0;
            for (unsigned i = 0; i < m.getHeight(); ++i) {
                T v = m(i, j);
                sum += (v >= 0 ? v : -v);
            }
            if (sum > norm) {
                norm = sum;
            }
        }

        return norm;
    }

    template<typename T>
    T normInf(const Matrix<T>& m) {
        T norm = 0;

        for (unsigned i = 0; i < m.getHeight(); ++i) {
            T sum = 0;
            for (unsigned j = 0; j < m.getWidth(); ++j) {
                T v = m(i, j);
                sum += (v >= 0 ? v : -v);
            }
            if (sum > norm) {
                norm = sum;
            }
        }

        return norm;
    }

    template<typename T>
    T normFrob(const Matrix<T>& m) {
        T sum = 0;

        for (unsigned i = 0; i < m.getHeight(); ++i) {
            for (unsigned j = 0; j < m.getWidth(); ++j) {
                T v = m(i, j);
                sum += v * v;
            }
        }

        return sqrt(sum);
    }

    /*template<typename T>
    Matrix<T> EigenValues(const Matrix<T>& m) {
        typename Matrix<T>::LDU ldu = LDUDecomposition(m);
    }*/

    template<typename T>
    Matrix<T> CompanionMatrix(basic_polynomial<T> polinom) {
        int deg = polinom.deg();
        if (deg < 0) {
            throw Matrix<T>::InvalidSizeException();
        }

        if (polinom.getValue(deg) != 1) {
            throw Matrix<T>::InvalidSizeException();
        }

        Matrix<T> m(deg, deg);
        for (int i = 0; i < deg; ++i) {
            m(i, deg - 1) = -1 * polinom.getValue(i);
            if (i < deg - 1) {
                m(i + 1, i) = 1;
            }
        }
        return m;
    }

    template<typename T>
    T trace(Matrix<T> matrix) {
        unsigned Width = matrix.getWidth();
        unsigned Height = matrix.getHeight();

        if (Width != Height) {
            throw Matrix<T>::SquareMatrixRequired();
        }

        T tr = 0;
        for (unsigned int i = 0; i < Width; ++i) {
            tr += matrix(i, i);
        }

        return tr;
    }

}
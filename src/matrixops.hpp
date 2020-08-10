#pragma once

#include <Eigen/Dense>

namespace khover {

struct rowops;
struct colops;

struct rowops {
    using dual_t = colops;

    template<class Derived>
    static inline
    std::size_t size(const Eigen::MatrixBase<Derived> &mat)
        noexcept
    {
        return static_cast<std::size_t>(mat.cols());
    }

    template<class Derived>
    static inline
    decltype(auto) at(Eigen::MatrixBase<Derived> &mat, std::size_t i)
        noexcept
    {
        return mat.row(i);
    }

    template<class Derived>
    static inline
    void swap(Eigen::MatrixBase<Derived> &mat, std::size_t i1, std::size_t i2)
        noexcept
    {
        mat.row(i1).swap(mat.row(i2));
    }

    template<class Derived, class T>
    static inline
    void scalar(Eigen::MatrixBase<Derived> &mat, std::size_t i, T a)
        noexcept
    {
        mat.row(i) *= a;
    }

    template <class Derived,class T>
    static inline
    void axpy(Eigen::MatrixBase<Derived> &mat, T a, std::size_t isrc, std::size_t itgt)
        noexcept
    {
        mat.row(itgt) += a * mat.row(isrc);
    }

    template <class Derived, class F>
    static inline
    std::size_t find_nonzero(
        const Eigen::MatrixBase<Derived> &mat,
        std::size_t i, F f)
        noexcept
    {
        std::size_t c = 0;

        for(; c < static_cast<std::size_t>(mat.cols()); ++c) {
            if (mat(i,c) != 0) {
                f(c, mat(i,c));
                break;
            }
        }

        // Since Eigen >= 3.4.0, STL iterator will come.
        // for(auto x : mat.row(i)) {
        //     if(x != 0) {
        //         f(c,x);
        //         break;
        //     }
        //     ++c;
        // }

        return c;
    }

    // The same as find_nonzero but range check is not performed.
    template <class Derived, class F>
    static inline
    std::size_t find_nonzero_unsafe(
        const Eigen::MatrixBase<Derived> &mat,
        std::size_t i, F f)
        noexcept
    {
        std::size_t c = 0;

        for(; c < static_cast<std::size_t>(mat.cols()); ++c) {
            if (mat.coeff(i,c) != 0) {
                f(c, mat.coeff(i,c));
                break;
            }
        }

        // Since Eigen >= 3.4.0, STL iterator will come.
        // for(auto x : mat.row(i)) {
        //     if(x != 0) {
        //         f(c,x);
        //         break;
        //     }
        //     ++c;
        // }

        return c;
    }
};

struct colops {
    using dual_t = rowops;

    template<class Derived>
    static inline
    std::size_t size(const Eigen::MatrixBase<Derived> &mat)
        noexcept
    {
        return static_cast<std::size_t>(mat.rows());
    }

    template<class Derived>
    static inline
    decltype(auto) at(Eigen::MatrixBase<Derived> &mat, std::size_t i)
        noexcept
    {
        return mat.col(i);
    }

    template<class Derived>
    static inline
    void swap(Eigen::MatrixBase<Derived> &mat, std::size_t j1, std::size_t j2)
        noexcept
    {
        mat.col(j1).swap(mat.col(j2));
    }

    template<class Derived, class T>
    static inline
    void scalar(Eigen::MatrixBase<Derived> &mat, std::size_t j, T a)
        noexcept
    {
        mat.col(j) = a * mat.col(j);
        //mat.col(j) *= a;
    }

    template <class Derived,class T>
    static inline
    void axpy(Eigen::MatrixBase<Derived> &mat, T a, std::size_t jsrc, std::size_t jtgt)
        noexcept
    {
        mat.col(jtgt) = mat.col(jtgt) + a * mat.col(jsrc);
    }

    template <class Derived, class F>
    static inline
    std::size_t find_nonzero(
        const Eigen::MatrixBase<Derived> &mat,
        std::size_t j, F f)
        noexcept
    {
        std::size_t r = 0;

        for(; r < static_cast<std::size_t>(mat.rows()); ++r) {
            if (mat(r,j) != 0) {
                f(r, mat(r,j));
                break;
            }
        }

        // for(auto x : mat.col(j)) {
        //     if(x != 0) {
        //         f(r,x);
        //         break;
        //     }
        //     ++r;
        // }
        return r;
    }

    template <class Derived, class F>
    static inline
    std::size_t find_nonzero_unsafe(
        const Eigen::MatrixBase<Derived> &mat,
        std::size_t j, F f)
        noexcept
    {
        std::size_t r = 0;

        for(; r < static_cast<std::size_t>(mat.rows()); ++r) {
            if (mat.coeff(r,j) != 0) {
                f(r, mat.coeff(r,j));
                break;
            }
        }
        return r;
    }
};

}

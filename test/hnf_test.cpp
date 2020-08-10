#include <iostream>
#include <cstdint>
#include <Eigen/Dense>

#include "hnf.hpp"

template<class Derived>
bool is_rowHNF(const Eigen::MatrixBase<Derived> &mat) noexcept
{
    for(std::size_t i = 0, j = 0; j < mat.cols() && i < mat.rows(); ++j) {
        for(std::size_t k = i+1; k < mat.rows(); ++k)
            if (mat(k,j) != 0)
                return false;
        if (mat(i,j) != 0) {
            for(std::size_t k = 0; k < i; ++k)
                if (mat(k,j) < 0 || mat(i,j) < mat(k,j))
                    return false;
            ++i;
        }
    }

    return true;
}

template<class Derived>
bool is_colHNF(const Eigen::MatrixBase<Derived> &mat) noexcept
{
    return is_rowHNF(mat.transpose().eval());
}

template<class Derived>
bool is_unimodular(const Eigen::MatrixBase<Derived> &mat) noexcept
{
    return std::abs(mat.determinant()) == 1;
}

int main(int argc, char*argv[])
{
    std::cout << "==========" << std::endl;
    std::cout << "Invalid size matrices should be rejected." << std::endl;
    std::cout << "==========" << std::endl;
    {
        using namespace khover;

        Eigen::Matrix<int,6,4> mat, a, b;
        Eigen::Matrix<int,4,6> c, d;
        Eigen::Matrix<int,3,3> e;

        auto fst = hnf_LLL<rowops>(mat,std::tie(b,c,d),std::tie(a,b));
        auto snd = hnf_LLL<rowops>(mat,std::tie(c,d),std::tie(a,b,c));
        auto trd = hnf_LLL<colops>(mat,std::tie(c,d,e),std::tie(a,b));
        auto fth = hnf_LLL<colops>(mat,std::tie(c,d),std::tie(e,a,b));

        if(fst || snd || trd || fth)
        {
            std::cout << "Failed to detect wrong size matrices" << std::endl;
            return -1;
        }
    }

    std::cout << "==========" << std::endl;
    std::cout << "Row HNF" << std::endl;
    std::cout << "==========" << std::endl;
    {
        Eigen::Matrix<std::int64_t, 3,4> mat0;
        mat0 <<
            1, 2, 3, 4,
            5, 6, 7, 8,
            9, 10, 11, 12;
        auto mat = mat0;
        auto u = Eigen::Matrix<std::int64_t,3,3>::Identity().eval();
        auto rk = khover::hnf_LLL<typename khover::rowops>(mat,std::tie(u),std::tuple<>{});

        if (!is_rowHNF(mat) || u*mat != mat0 || !is_unimodular(u) || !rk || *rk!=2)
        {
            std::cout << "mat=" << std::endl;
            std::cout << mat << std::endl;
            std::cout << std::endl;
            std::cout << "u=" << std::endl;
            std::cout << u << std::endl;
            std::cout << "rk=" << (rk ? *rk : -1) << std::endl;
            return -1;
        }
    }

    std::cout << "==========" << std::endl;
    std::cout << "Column HNF" << std::endl;
    std::cout << "==========" << std::endl;
    {
        Eigen::Matrix<std::int64_t, 5,4> mat0;
        mat0 <<
            1, 2, 4, 7,
            2, 3, 4, 5,
            6, 7, 8, 9,
            1, 1, 2, 3,
            5, 8, 13, 21;
        auto mat = mat0;
        Eigen::Matrix<std::int64_t,4,4> u = decltype(u)::Identity();
        Eigen::Matrix<std::int64_t,4,4> v = decltype(v)::Identity();
        auto rk = khover::hnf_LLL<typename khover::colops>(mat,std::tie(u),std::tie(v));

        if (!is_colHNF(mat) || mat*u != mat0
            || v*u != Eigen::Matrix<std::int64_t,4,4>::Identity()
            || !rk || *rk != 4)
        {
            std::cout << "mat=" << std::endl;
            std::cout << mat << std::endl;
            std::cout << std::endl;
            std::cout << "u=" << std::endl;
            std::cout << u << std::endl;
            std::cout << "rk=" << (rk ? *rk : -1) << std::endl;
            return -1;
        }
    }

    return 0;
}

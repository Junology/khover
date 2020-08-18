#include <iostream>
#include <cstdint>
#include <Eigen/Dense>

#include "hnf.hpp"

template<class Derived>
bool is_rowHNF(const Eigen::MatrixBase<Derived> &mat) noexcept
{
    for(int i = 0, j = 0; j < mat.cols() && i < mat.rows(); ++j) {
        for(int k = i+1; k < mat.rows(); ++k)
            if (mat(k,j) != 0)
                return false;
        if (mat(i,j) != 0) {
            for(int k = 0; k < i; ++k)
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
            std::cout << "u*mat=" << std::endl;
            std::cout << u*mat << std::endl;
            return -1;
        }
    }
    {
        Eigen::Matrix<std::int64_t, 3,4> mat0;
        mat0 <<
            0, 1, 2, 3,
            4, 5, 6, 7,
            8, 9,10,11;
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
            std::cout << "u*mat=" << std::endl;
            std::cout << u*mat << std::endl;
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
    {
        Eigen::Matrix<std::int64_t,3,3> mat0;
        mat0 <<
            1, 1, 0,
           -1, 0, 1,
            0,-1,-1;

        Eigen::Matrix<std::int64_t,3,3> u0;
        u0 <<
            1, 0, 0,
           -1, 1, 0,
            1, 0, 1;
        decltype(mat0) mat = mat0;
        decltype(u0) u = u0;
        Eigen::Matrix<std::int64_t,3,3> v = decltype(v)::Identity();
        auto rk = khover::hnf_LLL<typename khover::colops>(
            mat,std::tie(u),std::tie(v));

        if (!is_colHNF(mat) || mat*u != mat0*u0
            || v*u != u0
            || !rk || *rk != 2)
        {
            ERR_MSG(
                "mat=\n" << mat << std::endl
                << "u=\n" << u << std::endl
                << "v=\n" << v << std::endl
                << "rk=" << (rk ? *rk : -1) << std::endl
                << "mat0*u0=\n" << mat0*u0 << std::endl
                << "mat*u=\n" << mat*u << std::endl
                << "v*u=\n" << v*u
                );
            return -1;
        }
    }
    {
        Eigen::Matrix<std::int64_t, 3,4> mat0;
        mat0 <<
            0, 1, 2, 3,
            4, 5, 6, 7,
            8, 9,10,11;
        auto mat = mat0;
        auto u = Eigen::Matrix<std::int64_t,4,4>::Identity().eval();
        auto rk = khover::hnf_LLL<typename khover::colops>(mat,std::tie(u),std::tuple<>{});

        if (!is_colHNF(mat) || mat*u != mat0 || !is_unimodular(u) || !rk || *rk!=2)
        {
            ERR_MSG(
                "mat=\n" << mat << std::endl
                << "u=\n" << u << std::endl
                << "rk=" << (rk ? *rk : -1) << std::endl
                << "mat*u=\n" << mat*u
                );
            return -1;
        }
    }
    {
        Eigen::Matrix<std::int64_t, 9,27> mat0;
        mat0 <<
            1, 0, 0, 0, 0, 1, 0, 1, 3, 0, 0, 1, 0, 0, 0,-1, 0, 0, 0, 1, 3,-1, 0, 0,-3, 0, 0,
            0, 1, 0, 1, 0,-3, 0,-3,-8, 1, 0,-3, 0, 0, 0, 0,-1, 0, 0,-3,-8, 0,-1, 0, 0,-3, 0,
            0, 0, 1, 0, 1, 3, 1, 3, 6, 0, 1, 3, 0, 0, 0, 0, 0,-1, 1, 3, 6, 0, 0,-1, 0, 0,-3,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 3, 1, 3, 0, 0, 0, 3, 1, 3, 8, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0,-6, 0, 0,-8, 0, 0, 0, 0, 0,-8, 0, 8, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 6, 1, 3, 9, 0, 0, 0, 1, 3, 9, 0, 0, 8,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 0,-3, 0, 1, 0, 0, 0,-3, 0, 1,-6, 2, 6,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 1,-3,-3, 0, 0, 0, 1,-3,-3, 0,-12,-16,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,-1, 0, 1, 0, 0, 0, 0, 0, 1, 0, 2, 6, 6;
        auto mat = mat0;
        auto u = Eigen::Matrix<std::int64_t,27,27>::Identity().eval();
        auto v = Eigen::Matrix<std::int64_t,27,27>::Identity().eval();
        auto rk = khover::hnf_LLL<khover::colops>(mat,std::tie(u),std::tie(v));

        if (!is_colHNF(mat)
            || mat*u != mat0
            || !(u*v).isIdentity() ||
            !rk || *rk!=7)
        {
            ERR_MSG(
                "mat=\n" << mat << std::endl
                << "u=\n" << u << std::endl
                << "rk=" << (rk ? *rk : -1) << std::endl
                << "mat*u=\n" << mat*u
                );
            return -1;
        }
        std::cout << mat.leftCols(1+*rk) << std::endl;
    }

    std::cout << "==========" << std::endl;
    std::cout << "HNF of zero matrices" << std::endl;
    std::cout << "==========" << std::endl;
    {
        using matrix_t = Eigen::Matrix<std::int64_t,Eigen::Dynamic,Eigen::Dynamic>;
        matrix_t zmat;

        for(std::size_t r = 0; r < 5; ++r) {
            for(std::size_t c = 0; c < 5; ++c) {
                // test row HNF
                zmat = matrix_t::Zero(r,c);
                auto rowrk = khover::hnf_LLL<khover::rowops>(zmat, {}, {});
                if (!rowrk || *rowrk != 0 || !zmat.isZero()) {
                    std::cout << "row HNF of zero is not zero!" << std::endl;
                    std::cout << "rk=" << (rowrk ? *rowrk : -1) << std::endl;
                    std::cout << zmat << std::endl;
                    return -1;
                }
                // test column HNF
                zmat.setZero();
                auto colrk = khover::hnf_LLL<khover::colops>(zmat, {}, {});
                if (!colrk || *colrk!=0 || !zmat.isZero()) {
                    std::cout << "column HNF of zero is not zero!" << std::endl;
                    std::cout << "rk=" << (colrk ? *colrk : -1) << std::endl;
                    std::cout << zmat << std::endl;
                    return -1;
                }
            }
        }
    }

    return 0;
}

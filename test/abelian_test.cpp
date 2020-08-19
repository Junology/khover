#include <iostream>
#include <functional>

#include "abelian.hpp"

#include "debug/debug.hpp"

using namespace khover;

template <class T>
bool has_same_span(
    Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> m1,
    Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> m2 )
{
    hnf_LLL<colops>(m1, {}, {});
    hnf_LLL<colops>(m2, {}, {});

    return m1.rows() == m2.rows()
        && m1.cols() == m2.cols()
        && m1 == m2;
}

int main(int argc, char* argv[])
{
    {
        AbelianGroup a(
            (AbelianGroup::matrix_t(3,4) << 0,1,2,3,4,5,6,7,8,9,10,11).finished(),
            std::false_type{},
            3);
        auto shouldbe = (Eigen::Matrix<int64_t,5,1>() << 4,0,0,0,0).finished();

        a.reduce(std::tuple<>{}, std::tuple<>{});
        if(a.ngens() != 5
           || a.nrels() != 1
           || a.get_repmatrix() != shouldbe)
        {
            ERR_MSG(a.ngens() << "x" << a.nrels() << "\n" << a.get_repmatrix());
            return -1;
        }
    }
    {
        auto mat = (Eigen::Matrix<int64_t,4,5>() <<
                    6, 0,-4,-4, 2,
                    0, 1, 1, 2, 1,
                    2, 1, 5, 6, 7,
                    8, 2, 2, 4,10 ).finished();
        auto shouldbe = (Eigen::Matrix<int64_t,5,2>() <<
                         2, 0, 0, 16, 0, 0, 0, 0, 0, 0).finished();

        AbelianGroup a(mat, std::false_type{}, 2);

        a.reduce(std::tuple<>{}, std::tuple<>{});

        if(a.ngens() != 5
           || a.nrels() != 2
           || a.get_repmatrix() != shouldbe)
        {
            ERR_MSG(a.get_repmatrix());
            return -1;
        }
    }

    {
        auto mat = (Eigen::Matrix<int64_t,5,4>() <<
                    4, 0, 0, 0,
                    0, 1, 0, 0,
                    1, 0, 6, 0,
                    3, 0, 3, 12,
                    2, 0, 2, 8).finished();
        AbelianGroup a(mat, std::true_type{}, 2);
        // -> a = Z^2 + Z/96

        if(auto res = a.compute(); res.freerank != 3
           || std::accumulate(
               res.torsions.begin(), res.torsions.end(), 1, std::multiplies<int64_t>{}) != 96) {
            ERR_MSG("Wrong abelian group:" << std::endl << res.pretty());
            return -1;
        }
    }

    // Reduction with morphisms
    std::cout << "---" << std::endl;
    std::cout << "Reduction with morphisms" << std::endl;
    {
        auto mat = (Eigen::Matrix<int64_t,5,4>() <<
                    4, 0, 0, 0,
                    0, 1, 0, 0,
                    1, 0, 6, 0,
                    3, 0, 3, 12,
                    2, 0, 2, 8).finished();
        AbelianGroup a(mat, std::true_type{}, 2);
        // -> a = Z^2 + Z/96 by the previous test.
        // The torsion generator: [1 0 5 3 2]^t

        // Notice the this doesn't compile.
        // auto u = Eigen::Matrix<int64_t,7,7>::Identity().eval();
        AbelianGroup::matrix_t u = Eigen::Matrix<int64_t,7,7>::Identity().eval();
        AbelianGroup::matrix_t v = Eigen::Matrix<int64_t,7,7>::Identity().eval();
        a.reduce(std::tie(u), std::tie(v));

        auto u_shouldbe
            = (decltype(u)(u.rows(), 3) <<
               4, 0,-1,
               0, 0, 0,
               1, 6,-5,
               3, 3,-3,
               2, 2,-2,
               0, 0, 0,
               0, 0, 0).finished();

        if(u.leftCols(3) != u_shouldbe
           || v*u != AbelianGroup::matrix_t::Identity(6,6)) {
            ERR_MSG(
                "Wrong post-morphism transformation." << std::endl
                << "Matrix u:" << std::endl
                << u << std::endl
                << "Shouldbe" << std::endl
                << u_shouldbe << std::endl);
            return -1;
        }

        auto enh_mat = Eigen::Matrix<int64_t,7,4>::Zero().eval();
        enh_mat.topLeftCorner<5,4>().noalias() = mat;
        auto check_v
            = (decltype(v)(a.ngens(),a.nrels()+4) <<
               a.get_repmatrix(), v*enh_mat
                ).finished();
        if(auto rk = hnf_LLL<colops>(check_v, {}, {});
           !rk || *rk != 3
            )
        {
            std::cerr << "Wrong pre-morphism transformation." << std::endl;
            return -1;
        }
    }

    // Image computation test
    std::cout << "---" << std::endl;
    std::cout << "Image computation test" << std::endl;
    {
        auto mat = (Eigen::Matrix<int64_t,5,4>() <<
                    4, 0, 0, 0,
                    0, 1, 0, 0,
                    1, 0, 6, 0,
                    3, 0, 3, 12,
                    2, 0, 2, 8).finished();
        AbelianGroup a(mat, std::true_type{}, 2);
        // -> a = Z^2 + Z/96 by the previous test.
        // The torsion generator: [1 0 5 3 2]^t

        auto u
            = (AbelianGroup::matrix_t(7,6) <<
               4, 0,-1, 0, 0, 0,
               0, 0, 0, 0, 0, 0,
               1, 6,-5,-2, 0, 0,
               3, 3,-3,-1, 0, 0,
               2, 2,-2,-1, 0, 0,
               0, 0, 0, 0, 1, 0,
               0, 0, 0, 0, 0, 1).finished();
        // It is verified that u induces the identity.

        auto coeffmat = Eigen::Matrix<int64_t,6,6>::Zero().eval();
        coeffmat.diagonal() << 1, 1, 16, 3, 0, 5;
        coeffmat.topLeftCorner(2, 2) << -9, 8, 5, 11;

        auto enh_mat = Eigen::Matrix<int64_t,7,4>::Zero().eval();
        enh_mat.topLeftCorner<5,4>().noalias() = mat;
        auto [imgrp,b] = *a.image(u*coeffmat,std::true_type{});

        if(!has_same_span<int64_t>(b*imgrp.get_repmatrix(), enh_mat)) {
            std::cerr << "Wrong image inclusion:" << std::endl;
            std::cerr << "b=" << std::endl << b << std::endl;
            std::cerr << "b*rep=" << std::endl
                      <<  b*imgrp.get_repmatrix() << std::endl;
            std::cerr << "Shouldbe" << std::endl << enh_mat << std::endl;
            return -1;
        }

        if(auto im = imgrp.compute();
           im.freerank != 2
           || std::accumulate(
               im.torsions.begin(), im.torsions.end(), 1,
               std::multiplies<int64_t>{}) != 6) // 6 = 96/16
        {
            std::cout << u*coeffmat << std::endl;
            ERR_MSG("Wrong image." << std::endl << im.pretty());
            return -1;
        }
    }

    return 0;
}

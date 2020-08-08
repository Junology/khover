#include <iostream>
#include <cstdint>
#include <utility>

#include <Eigen/Dense>

#include "utils.hpp"
#include "matrixops.hpp"

template <class T, class Derived>
void manipulate(Eigen::MatrixBase<Derived> &mat)
{
    T::scalar(mat,0,3);
    T::swap(mat,2,3);
    T::axpy(mat,5,1,0);

    T::dual_t::axpy(mat,7,1,2);
    T::dual_t::axpy(mat,2,1,3);
}

int main(int argc, char*argv[])
{
    Eigen::Matrix<int64_t,4,4> mat = Eigen::Matrix<int64_t,4,4>::Identity();

    manipulate<typename khover::rowops>(mat);
    {
        Eigen::Matrix<int64_t,4,4> shouldbe;
        shouldbe <<
            3, 5, 35, 10,
            0, 1,  7,  2,
            0, 0,  0,  1,
            0, 0,  1,  0;

        if(mat != shouldbe) {
            std::cout << "First test failed" << std::endl;
            std::cout << mat << std::endl;
            return -1;
        }
    }

    manipulate<typename khover::colops>(mat);
    {
        Eigen::Matrix<int64_t,4,4> shouldbe;
        shouldbe <<
            34, 5, 10, 35,
             5, 1,  2,  7,
            35, 7, 15, 49,
            10, 2,  4, 15;

        if(mat != shouldbe) {
            std::cout << "Second test failed" << std::endl;
            std::cout << mat << std::endl;
            return -1;
        }
    }

    Eigen::Matrix<int64_t,2,5> mat2;
    mat2 <<
        1, 2, 3, 4, 5,
        0, 1, 1, 2, 3;

    khover::for_each_tuple(
        std::tie(mat,mat2),
        [](auto& m){ khover::rowops::axpy(m,-2,1,0); } );
    {
        Eigen::Matrix<int64_t,4,4> shouldbe1;
        shouldbe1 <<
            24, 3,  6, 21,
             5, 1,  2,  7,
            35, 7, 15, 49,
            10, 2,  4, 15;

        Eigen::Matrix<int64_t,2,5> shouldbe2;
        shouldbe2 <<
            1, 0, 1, 0, -1,
            0, 1, 1, 2, 3;
        if(mat != shouldbe1 || mat2 != shouldbe2) {
            std::cout << "Third test failed" << std::endl;
            std::cout << mat << std::endl;
            return -1;
        }
    }

    return 0;
}


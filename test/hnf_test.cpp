#include <iostream>
#include <cstdint>
#include <Eigen/Dense>

#include "hnf.hpp"

int main(int argc, char*argv[])
{
    {
        Eigen::Matrix<std::int64_t, 3,4> mat;
        mat <<
            1, 2, 3, 4,
            5, 6, 7, 8,
            9, 10, 11, 12;
        auto u = Eigen::Matrix<std::int64_t,3,3>::Identity();
        khover::hnf_row(u,mat);

        Eigen::Matrix<std::int64_t, 3,4> shouldbeM;
        shouldbeM <<
            1, 2, 3, 4,
            0, 4, 8, 12,
            0, 0, 0, 0;
        Eigen::Matrix<std::int64_t,3,3> shouldbeU;
        shouldbeU <<
            1, 0, 0,
            5, -1, 0,
            9, -2, 1;

        if (mat != shouldbeM || u != shouldbeU) {
            std::cerr << "Failed!" << std::endl;
            return -1;
        }
    }

    return 0;
}

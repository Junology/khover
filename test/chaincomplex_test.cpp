#include <iostream>

#include "chaincomplex.hpp"
#include "debug/debug.hpp"

using namespace khover;

using matrix_t = ChainIntegral::matrix_t;

int main(int argc, char* argv[])
{
    ///*
    {
        matrix_t mata =
            (matrix_t(3,6) <<
             1, 1, 1, 1, 0, 0,
             -1,-1,0, 0, 1, 1,
             0, 0,-1,-1,-1,-1).finished();
        matrix_t matb =
            (matrix_t(6,3) <<
             1, 0, 1,
             0, 1, 0,
             -1,0, 0,
             0,-1,-1,
             1, 1, 0,
             0, 0, 1).finished();

        if(!(mata*matb).isZero()) {
            ERR_MSG("Wrong matrix multiplication.");
            return -1;
        }
    }
    //*/
    {
        ChainIntegral ch(
            -1,
            Eigen::Matrix<int64_t,0,1>::Zero());
        ch.prepend(
            (matrix_t(1,3) << 1, 1, 1).finished());
        ch.prepend(
            (matrix_t(3,6) <<
             1, 1, 1, 1, 0, 0,
             -1,-1,0, 0, 1, 1,
             0, 0,-1,-1,-1,-1).finished());
        ch.prepend(
            (matrix_t(6,3) <<
             1, 0, 1,
             0, 1, 0,
             -1,0, 0,
             0,-1,-1,
             1, 1, 0,
             0, 0, 1).finished());

        auto h = ch.compute();

        if(h.size() != 5) {
            ERR_MSG("Wrong length of homology groups.");
            return -1;
        }

        if(auto h0 = h[0].compute();
           h0.first != 0 || !h0.second.empty()) {
            ERR_MSG("Wrong -1st homology: " << h0);
            return -1;
        }
        if(auto h1 = h[1].compute();
           h1.first != 0 || !h1.second.empty()) {
            ERR_MSG("Wrong 0th homology: " << h1);
            return -1;
        }
        if(auto h2 = h[2].compute();
           h2.first != 0 || !h2.second.empty()) {
            ERR_MSG("Wrong 1st homology: " << h2);
            return -1;
        }
        if(auto h3 = h[3].compute();
           h3.first != 1 || !h3.second.empty()) {
            ERR_MSG("Wrong 2nd homology: " << h3);
            return -1;
        }
        if(auto h4 = h[4].compute();
           h4.first != 0 || !h4.second.empty()) {
            ERR_MSG("Wrong 3rd homology: " << h4);
            return -1;
        }
    }

    return 0;
}

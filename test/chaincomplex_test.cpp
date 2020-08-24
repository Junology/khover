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

    // Rank computation
    {
        ChainIntegral ch(
            -1,
            (matrix_t(1,3) << 1, 1, 1).finished(),
            (matrix_t(3,6) <<
             1, 1, 1, 1, 0, 0,
             -1,-1,0, 0, 1, 1,
             0, 0,-1,-1,-1,-1).finished(),
            (matrix_t(6,3) <<
             1, 0, 1,
             0, 1, 0,
             -1,0, 0,
             0,-1,-1,
             1, 1, 0,
             0, 0, 1).finished()
            );
        if(ch.rank(-2) != 0
           || ch.rank(-1) != 1
           || ch.rank(0) != 3
           || ch.rank(1) != 6
           || ch.rank(2) != 3
           || ch.rank(3) != 0)
        {
            ERR_MSG("Wrong ranks:{"
                    << ch.rank(-1) << ", "
                    << ch.rank(0) << ", "
                    << ch.rank(1) << ", "
                    << ch.rank(2) << ", "
                    << ch.rank(3) << "}");
            return EXIT_FAILURE;
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
           h0.freerank != 0 || !h0.torsions.empty()) {
            ERR_MSG("Wrong -1st homology: " << h0.pretty());
            return -1;
        }
        if(auto h1 = h[1].compute();
           h1.freerank != 0 || !h1.torsions.empty()) {
            ERR_MSG("Wrong 0th homology: " << h1.pretty());
            return -1;
        }
        if(auto h2 = h[2].compute();
           h2.freerank != 0 || !h2.torsions.empty()) {
            ERR_MSG("Wrong 1st homology: " << h2.pretty());
            return -1;
        }
        if(auto h3 = h[3].compute();
           h3.freerank != 1 || !h3.torsions.empty()) {
            ERR_MSG("Wrong 2nd homology: " << h3.pretty());
            return -1;
        }
        if(auto h4 = h[4].compute();
           h4.freerank != 0 || !h4.torsions.empty()) {
            ERR_MSG("Wrong 3rd homology: " << h4.pretty());
            return -1;
        }
    }

    return 0;
}

#include <iostream>

#include "chaincomplex.hpp"

using namespace khover;

using matrix_t = ChainIntegral::matrix_t;

template<class T>
std::ostream& operator<<(std::ostream& os, std::vector<T> const& vec) {
    if (vec.empty()) {
        os << "{}";
        return os;
    }

    os << "{" << vec.front();
    for(auto itr = std::next(std::begin(vec)); itr != std::end(vec); ++itr)
        os << ", " << *itr;
    os << "}";
    return os;
}

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
            std::cerr << "Wrong matrix multiplication." << std::endl;
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
            std::cerr << "Wrong length of homology groups." << std::endl;
            return -1;
        }

        if(auto h0 = h[0].compute();
           h0.first != 0 || !h0.second.empty()) {
            std::cerr << "Wrong -1st homology: " << h0 << std::endl;
            return -1;
        }
        if(auto h1 = h[1].compute();
           h1.first != 0 || !h1.second.empty()) {
            std::cerr << "Wrong 0th homology: " << h1 << std::endl;
            return -1;
        }
        if(auto h2 = h[2].compute();
           h2.first != 0 || !h2.second.empty()) {
            std::cerr << "Wrong 1st homology: " << h2 << std::endl;
            return -1;
        }
        if(auto h3 = h[3].compute();
           h3.first != 1 || !h3.second.empty()) {
            std::cerr << "Wrong 2nd homology: " << h3 << std::endl;
            return -1;
        }
        if(auto h4 = h[4].compute();
           h4.first != 0 || !h4.second.empty()) {
            std::cerr << "Wrong 3rd homology: " << h4 << std::endl;
            return -1;
        }
    }

    return 0;
}

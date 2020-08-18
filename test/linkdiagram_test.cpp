#include <iostream>
#include <algorithm>

#include "linkdiagram.hpp"

using namespace khover;

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
    auto six_two = read_gauss_code(
        {1,-2, 5,-6, 3,-1, 2,-4, 6,-5, 4,-3},
        {std::make_pair(1,false)});

    if(!six_two) {
        std::cerr << "Failed to load a knot." << std::endl;
        return -1;
    }

    if (six_two->nnegative() != 4) {
        std::cerr << "Wrong number of negative crossings: "
                  << six_two->nnegative() << " should be " << 4 << std::endl;
        return -1;
    }

    if (six_two->npositive() != 2) {
        std::cerr << "Wrong number of positive crossings: "
                  << six_two->npositive() << " should be " << 2 << std::endl;
        return -1;
    }

    for(unsigned long s = 0; s < 0b1000000 ; ++s) {
        if(six_two->cohDegree(s) != static_cast<int>(std::bitset<6>(s).count()) - 4)
        {
            std::cerr << "Wrong cohomological degree for " << std::bitset<6>(s) << std::endl;
            return -1;
        }
    }

    // Smoothing test: all 0
    {
        auto [ncomp, comptbl] = six_two->smoothing(0b000000u);
        std::vector<unsigned char> shouldbe = {0, 1, 2, 3, 4, 0, 1, 3, 2, 1, 3, 4};

        if (ncomp != 5
            || !std::equal(comptbl.begin(), comptbl.end(), shouldbe.begin())) {
            std::cerr << "Wrong components:" << std::endl;
            for(auto c : comptbl) {
                std::cout << static_cast<int>(c) << "; ";
            }
            std::cout << std::endl;
            return -1;
        }
    }

    // Smoothing test: 0b000001
    {
        auto [ncomp, comptbl] = six_two->smoothing(0b000001u);
        std::vector<unsigned char> shouldbe = {0, 1, 2, 3, 0, 0, 1, 3, 2, 1, 3, 0};

        if (ncomp != 4
            || !std::equal(comptbl.begin(), comptbl.end(), shouldbe.begin())) {
            std::cerr << "Wrong components:" << std::endl;
            for(auto c : comptbl) {
                std::cout << static_cast<int>(c) << "; ";
            }
            std::cout << std::endl;
            return -1;
        }
    }

    // Smoothing test: 0b000010
    {
        auto [ncomp, comptbl] = six_two->smoothing(0b000010u);
        std::vector<unsigned char> shouldbe = {0, 0, 1, 2, 3, 0, 0, 2, 1, 0, 2, 3};

        if (ncomp != 4
            || !std::equal(comptbl.begin(), comptbl.end(), shouldbe.begin())) {
            std::cerr << "Wrong components:" << std::endl;
            for(auto c : comptbl) {
                std::cerr << static_cast<int>(c) << "; ";
            }
            std::cerr << std::endl;
            std::cerr << "Shouldbe:" << std::endl;
            for(auto c : shouldbe) {
                std::cerr << static_cast<int>(c) << "; ";
            }
            std::cerr << std::endl;
            return -1;
        }
    }

    // Smoothing test: 0b000011
    {
        auto [ncomp,comptbl] = six_two->smoothing(0b000011u);
        std::vector<unsigned char> shouldbe = {0, 0, 1, 2, 0, 0, 0, 2, 1, 0, 2, 0};

        if (ncomp != 3
            || !std::equal(comptbl.begin(), comptbl.end(), shouldbe.begin())) {
            std::cerr << "Wrong components:" << std::endl;
            for(auto c : comptbl) {
                std::cout << static_cast<int>(c) << "; ";
            }
            std::cout << std::endl;
            return -1;
        }
    }

    // Adjacent state test:
    {
        state_t st0 = 0b010100;
        state_t st1 = 0b010101;
        state_t st2 = 0b011101;
        state_t st3 = 0b111101;

        if(six_two->stateCoeff(st0, st1) != 1
           || six_two->stateCoeff(st1, st0) != 0)
        {
            std::cerr << "Wrong state adjacency: " << std::endl
                      << st0 << std::endl << st1 << std::endl
                      << "Coeffs = ("
                      << six_two->stateCoeff(st0, st1) << ","
                      << six_two->stateCoeff(st1, st0) << ")"
                      << std::endl;
            return -1;
        }

        if(six_two->stateCoeff(st1, st2) != -1
           || six_two->stateCoeff(st2, st1) != 0)
        {
            std::cerr << "Wrong state adjacency: " << std::endl
                      << st1 << std::endl << st2 << std::endl
                      << "Coeffs = ("
                      << six_two->stateCoeff(st1, st2) << ","
                      << six_two->stateCoeff(st2, st1) << ")"
                      << std::endl;
            return -1;
        }

        if(six_two->stateCoeff(st2, st3) != 1
           || six_two->stateCoeff(st3, st2) != 0)
        {
            std::cerr << "Wrong state adjacency: " << std::endl
                      << st2 << std::endl << st3 << std::endl
                      << "Coeffs = ("
                      << six_two->stateCoeff(st2, st3) << ","
                      << six_two->stateCoeff(st3, st2) << ")"
                      << std::endl;
            return -1;
        }
    }

    // Skew-commutativity of cube.
    for(std::size_t k = 0; k < 32; ++k) {
        state_t st00 = k;

        for(std::size_t i = 0; i < 5; ++i) {
            if(st00.test(i))
                continue;

            for(std::size_t j = i+1; j < 6; ++j) {
                if(st00.test(j))
                    continue;

                state_t st01 = st00, st10 = st00, st11 = st00;
                st01.set(i);
                st10.set(j);
                st11.set(i);
                st11.set(j);

                if(six_two->stateCoeff(st00, st01)*six_two->stateCoeff(st01, st11)
                   != -six_two->stateCoeff(st00, st10)*six_two->stateCoeff(st10, st11))
                {
                    std::cerr << "Cube is not skew-commutative:" << std::endl
                              << st00 << std::endl
                              << "(i,j)=(" << i << "," << j << ")" << std::endl;
                    return -1;
                }
            }
        }
    }
    {
        auto trefoil = read_gauss_code(
            {1, -3, 2, -1, 3, -2}, {std::make_pair(1,false)});
        if(auto comps = trefoil->smoothing(state_t{0b011});
           comps.first != 1
           || comps.second[0] != 0
           || comps.second[1] != 0
           || comps.second[2] != 0
           || comps.second[3] != 0
           || comps.second[4] != 0
           || comps.second[5] != 0
            )
        {
            std::cerr << "Wrong components:" << std::endl;
            for(auto c : comps.second) {
                std::cout << static_cast<int>(c) << "; ";
            }
            std::cout << std::endl;
            return -1;
        }
    }
    return 0;
}

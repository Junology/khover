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
        {1,-2, 5,-6, 3,-1, 2,-4, 6,-5, 4,-3}, {false});

    if(!six_two) {
        std::cerr << "Failed to load a knot." << std::endl;
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
        auto comptbl = six_two->smoothing(0b000000u);
        std::vector<unsigned char> shouldbe = {0, 1, 2, 3, 4, 0, 1, 3, 2, 1, 3, 4};

        if (!std::equal(comptbl.begin(), comptbl.end(), shouldbe.begin())) {
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
        auto comptbl = six_two->smoothing(0b000001u);
        std::vector<unsigned char> shouldbe = {0, 1, 2, 3, 0, 0, 1, 3, 2, 1, 3, 0};

        if (!std::equal(comptbl.begin(), comptbl.end(), shouldbe.begin())) {
            std::cerr << "Wrong components:" << std::endl;
            for(auto c : comptbl) {
                std::cout << static_cast<int>(c) << "; ";
            }
            std::cout << std::endl;
            return -1;
        }
    }

    // Smoothing test: 0b000011
    {
        auto comptbl = six_two->smoothing(0b000011u);
        std::vector<unsigned char> shouldbe = {0, 0, 2, 3, 0, 0, 0, 3, 2, 0, 3, 0};

        if (!std::equal(comptbl.begin(), comptbl.end(), shouldbe.begin())) {
            std::cerr << "Wrong components:" << std::endl;
            for(auto c : comptbl) {
                std::cout << static_cast<int>(c) << "; ";
            }
            std::cout << std::endl;
            return -1;
        }
    }

    return 0;
}

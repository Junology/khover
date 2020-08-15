#include <iostream>
#include <vector>
#include <algorithm>

#include "states.hpp"

using namespace khover;

int main(int argc, char* argv[])
{
    constexpr std::size_t nbits = 9;
    constexpr std::size_t popcnt = 4;
    constexpr std::size_t binom = (9*8*7*6)/(1*2*3*4);

    std::vector<std::bitset<nbits>> tbl;

    for(std::size_t i = 0; i < binom; ++i) {
        tbl.push_back(bitsWithPop<nbits>(popcnt, i));
    }

    for(std::size_t n = 0; n < binom; ++n) {
        auto num = std::count(
            std::begin(tbl), std::end(tbl),
            std::bitset<nbits>(n));
        if (std::bitset<nbits>(n).count() != popcnt && num != 0)
        {
            std::cerr << "Wrong value emitted: "
                      << std::bitset<nbits>(n) << std::endl;
            return -1;
        }
        if (std::bitset<nbits>(n).count() == popcnt && num != 1)
        {
            std::cerr << "Wrong number of elements: "
                      << std::bitset<nbits>(n) << std::endl;
            return -1;
        }
    }

    for(std::size_t i = 0; i < binom; ++i) {
        if (bitsIndex(tbl[i]) != i) {
            std::cerr << "Wrong index on the bit sequence:" << std::endl
                      << tbl[i] << std::endl;
            return -1;
        }
    }

    return 0;
}

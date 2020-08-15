/*!
 * \file states.hpp
 * \author Jun Yoshida
 * \copyright (c) 2020 Jun Yoshida.
 * The project is released under the 2-clause BSD License.
 * \date August, 2020: created
 */

#pragma once

#include <cstdint>
#include <bitset>

#include "utils.hpp"

namespace khover {

enum Limits : std::size_t {
    max_crosses = 64,
    max_arcs = 64,
    max_components = max_arcs,
};

/*************
 *** Types ***
 *************/
//! The type representing states.
using state_t = std::bitset<max_crosses>;

//! The type indexing connected components.
using component_t = std::uint8_t;

//! The type representing enhancements on states.
using enhancement_t = std::bitset<max_components>;


/*****************
 *** Functions ***
 *****************/
//! Enumerate all the bitsets with constant pop-counts.
template<std::size_t n>
std::bitset<n>
bitsWithPop(std::size_t popcnt, std::size_t index)
    noexcept
{
    if(popcnt == 0)
        return std::bitset<n>{0u};

    if(index == 0)
        return (~std::bitset<n>{0u}) >> (n-popcnt);

    if(popcnt == 1)
        return std::bitset<n>{1u} << index;

    std::size_t highest = popcnt-1;
    std::size_t binomdiff = 1;
    while(index >= binomdiff) {
        index -= binomdiff;
        ++highest;
        binomdiff *= highest;
        binomdiff /= (highest+1)-popcnt;
    }

    return (std::bitset<n>{1u} << highest) | bitsWithPop<n>(popcnt-1, index);
}

//! Get the index of a bitset used in bitsWithPop function.
template<std::size_t n>
std::size_t
bitsIndex(std::bitset<n> bits)
    noexcept
{
    std::size_t popcnt = bits.count();
    std::size_t result = 0;

    for(std::size_t i = n-1; i >= 1 && popcnt > 0; --i) {
        if(bits.test(i)) {
            result += binom(i,popcnt);
            --popcnt;
        }
    }

    return result;
}

} // end namespace khover

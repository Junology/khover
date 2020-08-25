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

//! The class to carry data associated to a smoothing
//! \tparam T The type associating components to arcs.
//! It must satisfy
//! \code
//!   static_assert(
//!     std::is_same_v<
//!       std::decay_t<decltype(std::declval<T>()[std::declval<std::size_t>()])>,
//!       component_t>,
//!       "...");
//! \endcode
template <class T>
struct Smoothing {
    static_assert(
        std::is_same_v<
        std::decay_t<decltype(std::declval<T>()[std::declval<std::size_t>()])>,
        component_t
        >,
        "The template parameter T must have operator[] admitting indices of type std::size_t and returning component_t.");

    //! The state associated to the smoothing.
    state_t state;

    //! The number of components.
    std::size_t ncomp;

    //! Assignment of components to each arc.
    T arccomp;
};

//! The class to carry data associated to a smoothing and twisting of arcs.
template <class T>
struct SmoothTw : Smoothing<T> {
    std::bitset<max_arcs> twist;
};

/************************************
 *** Constants for bit operations ***
 ************************************/
template<std::size_t n>
inline constexpr
std::bitset<n> mask{~0lu};

template<std::size_t n>
inline constexpr
std::size_t maskbits = std::min(
    mask<n>.size(),
    static_cast<std::size_t>(
        std::numeric_limits<unsigned long long int>::digits));


/*************************************************
 *** Utility functions for state manipulations ***
 *************************************************/

//! Generate low-cut window
template <std::size_t n>
static inline
std::bitset<n> low_window(std::size_t lbits)
    noexcept
{
    return mask<n> >> (maskbits<n> - lbits);
}

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

//! Get the predecessor with the same popcount.
//! \pre bitsIndex(bits) > 0.
template<std::size_t n>
inline std::bitset<n>
predWithPop(std::bitset<n> const& bits)
    noexcept
{
    return bitsWithPop<n>(
        bits.count(),
        bitsIndex(bits)-1
        );
}


//! Insert a bit in a bitset.
//! \param bs The bitset to be inserted to.
//! \param pos The position where a new bit is inserted.
//! \param flag The newly inserted bit.
template<std::size_t n>
inline
std::bitset<n>
insertBit(std::bitset<n> const& bs, std::size_t pos, bool flag)
{
    return (bs & low_window<n>(pos))
        | (std::bitset<n>{static_cast<unsigned long>(flag)} << pos)
        | ((bs & ~low_window<max_crosses>(pos)) << 1);
}

} // end namespace khover

/*!
 * \file debug.cpp
 * \author Jun Yoshida
 * \copyright (c) 2020 Jun Yoshida.
 * The project is released under the 2-clause BSD License.
 * \date August, 2020: created
 */

#pragma once

#if not defined NDEBUG

#include <iostream>
#include <utility>
#include <vector>
#include <array>

#define DBG_MSG(x) std::cout << __FILE__ << ":" << __LINE__ << std::endl << x << std::endl

#define ERR_MSG(x) std::cerr << __FILE__ << ":" << __LINE__ << std::endl << x << std::endl

template<class T, class U>
inline std::ostream& operator<<(std::ostream& os, std::pair<T,U> p) {
    os << "(" << p.first << "," << p.second << ")";
    return os;
}

template<class T>
inline std::ostream& operator<<(std::ostream& os, std::vector<T> const& vec) {
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

template <>
inline std::ostream& operator<<(std::ostream& os, std::vector<uint8_t> const& vec) {
    if (vec.empty()) {
        os << "{}";
        return os;
    }

    os << "{" << static_cast<int>(vec.front());
    for(auto itr = std::next(std::begin(vec)); itr != std::end(vec); ++itr)
        os << ", " << static_cast<int>(*itr);
    os << "}";
    return os;
}

template<class T, std::size_t n>
inline std::ostream& operator<<(std::ostream& os, std::array<T,n> const& arr) {
    os << "{" << arr[0];
    for(std::size_t i = 1; i < n; ++i)
        os << ", " << arr[i];
    os << "}";
    return os;
}

#else

#define DBG_MSG(x)
#define ERR_MSG(x)

#endif

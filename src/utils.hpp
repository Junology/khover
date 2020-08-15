/*!
 * \file utils.hpp
 * \author Jun Yoshida
 * \copyright (c) 2020 Jun Yoshida.
 * The project is released under the 2-clause BSD License.
 * \date August, 2020: created
 */

#pragma once

#include <type_traits>
#include <utility>
#include <tuple>
#include <numeric>
#include <cstdlib>

namespace khover{

/****************************!
 *** \section type traits ***
 ****************************/
/*!
 * \class is_pubbase_of_template
 * Check if a class is publically derived from a template class
 */
template <template <class...> class TBase, class Derived>
class is_pubbase_of_template_impl {
private:
    template <class... Ts>
    static std::true_type is_pubbase_of_template_implfunc(TBase<Ts...> const*);
    // C style variadic length argument.
    static std::false_type is_pubbase_of_template_implfunc(...);

public:
    using type = decltype(is_pubbase_of_template_implfunc(std::declval<Derived*>()));
};

template <template <class...> class TBase, class Derived>
struct is_pubbase_of_template : public is_pubbase_of_template_impl<TBase,Derived>::type {};
//using is_pubbase_of_template = typename is_pubbase_of_template_impl<TBase, Derived>::type;


/**********************************!
 *** \section Utility functions ***
 **********************************/

//! The identity function
template <class T>
constexpr T&& identity(T&& t) noexcept
{
    return std::forward<T>(t);
}

//! Traverse tuples
template <class Tuple, class FuncT, class C=std::integral_constant<std::size_t, 0>>
constexpr void for_each_tuple(Tuple &&t, FuncT f, C = {}) noexcept
{
    if constexpr(C::value < std::tuple_size<std::decay_t<Tuple>>::value) {
        f(std::get<C::value>(t));
        for_each_tuple(t, f, std::integral_constant<std::size_t,C::value+1>{});
    }
}

//! Fold a tuple with a binary operator.
template <
    class T, class Tuple, class BinOp,
    class C=std::integral_constant<std::size_t,0>
    >
constexpr T foldl_tuple(T head, const Tuple &t, BinOp b, C = {}) noexcept
{
    if constexpr(C::value < std::tuple_size<Tuple>::value) {
        return foldl_tuple(
            b(head, std::get<C::value>(t)), t, b,
            std::integral_constant<std::size_t,C::value+1>{});
    }
    else {
        return head;
    }
}


/*!
 * \function array_to_tuple_with
 * Convert array to tuple.
 * If the element type is std::reference_wrapper<T>, the associated type is T& (see std::make_tuple for more details).
 */
namespace _impl {

template <class T, std::size_t n, class F, std::size_t...is>
constexpr auto map_to_tuple_impl(
    T const (&arr)[n], F&& f,
    std::index_sequence<is...>) noexcept
{
    return std::make_tuple(f(arr[is])...);
}

} // end namespace _impl

template <class T, std::size_t n, class F>
constexpr auto map_to_tuple(T const (&arr)[n], F&& f) noexcept
{
    return _impl::map_to_tuple_impl(arr, std::forward<F>(f),
                                    std::make_index_sequence<n>());
}


/*******************************!
 *** \section Math functions ***
 *******************************/

//! Integral division with "truncated toward -infinity."
template <class T>
constexpr
T floor_div(T a, T b) noexcept // division-by-zero is not an exception.
{
    auto q = a/b;
    auto r = a%b;
    return (r<0) ? q-1 : q;
}

//! simple square
template <class T>
inline constexpr T sqpow(T x) noexcept { return x*x; }

//! Compile-time non-negative integer power
//! It is verified that this function is faster than std::pow.
template<class T>
inline constexpr T cipow(T x, unsigned int n, T init = 1) noexcept
{
    T result = init;

    while(n) {
        if(n&0x1)
            result *= x;
        x *= x;
        n >>= 1;
    }

    return result;

    // The following version is faster than the above under
    //   > g++ --std=c++14 -O2 -march=native
    /*
    if (n==0)
        return 1;
    else
        return ((n&0x1) ? x : 1)*cipow(x*x,n>>1);
    */
}

inline constexpr
std::size_t
binom(std::size_t n, std::size_t k)
    noexcept
{
    if (k < 0 || k > n)
        return 0;

    if (k > n-k)
        k = n-k;

    if (k == 0)
        return 1;

    std::size_t result = 1;

    for (std::size_t i = 0; i < k; ++i)
    {
        result *= n-i;
        result /= (i+1);
    }

    return result;
}

}


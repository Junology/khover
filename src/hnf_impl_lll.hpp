/*!
 * \file hnf_impl_lll.hpp
 * \author Jun Yoshida
 * \copyright (c) 2020 Jun Yoshida.
 * The project is released under the 2-clause BSD License.
 * \date August, 2020: created
 */

#pragma once

#include <valarray>
#include <tuple>
#include <cmath>
#include <Eigen/Dense>

#include "utils.hpp"
#include "matrixops.hpp"

//* For debug
#include <iostream>

#define DBG_MSG std::cerr << __FILE__ << ":" << __LINE__ << std::endl

template<class T, class U>
std::ostream& operator<<(std::ostream& os, std::pair<T,U> p) {
    os << "(" << p.first << "," << p.second << ")";
    return os;
}
// */


namespace khover {

namespace _impl_LLL {

class Lambda_t {
public:
    using underlying_t = long double;
    static constexpr double delta = 0.75;

private:
    std::size_t m_size;
    std::valarray<underlying_t> m_diag;

    //! strictly lower triangle matrix with row-major
    std::valarray<underlying_t> m_offDiag;

protected:
    static constexpr std::size_t rhead(std::size_t i) noexcept {
        return i*(i-1)/2;
    }

public:
    Lambda_t(std::size_t n) noexcept
        : m_size(n),
          m_diag(static_cast<underlying_t>(1),n),
          m_offDiag(static_cast<underlying_t>(0),n*(n-1)/2)
    {}

    //! Element access (R/W)
    //! \pre i >= j >= 0
    underlying_t& operator()(std::size_t i, std::size_t j) {
        return i==j ? m_diag[i] : m_offDiag[rhead(i)+j];
    }

    //! Element access (Read only)
    underlying_t const & operator()(std::size_t i, std::size_t j) const {
        return const_cast<Lambda_t*>(this)->operator()(i,j);
    }

    //! \pre isrc != itgt
    void axpy(underlying_t coeff, std::size_t isrc, std::size_t itgt) noexcept {
        auto imin = std::min(isrc,itgt);

        if (imin <= 0) return;

        std::slice_array<underlying_t> tgt
            =  m_offDiag[std::slice(rhead(itgt),imin,1)];

        std::valarray<underlying_t> qsrc(
            m_offDiag[std::slice(rhead(isrc),imin,1)]);

        tgt += std::valarray<underlying_t>(coeff,imin) * qsrc;
    }

    //! Swap sliced rows of the size of the smaller row.
    void slice_swap(std::size_t i1, std::size_t i2) noexcept
    {
        auto imin = std::min(i1,i2);
        for(size_t j = 0; j < imin; ++j)
            std::swap(m_offDiag[rhead(i1)+j], m_offDiag[rhead(i2)+j]);
    }

    //! Update lambda in SWAP step in the HMM algorithm.
    void swap(std::size_t k) noexcept {
        for (size_t i = k+2; i < m_size-1; ++i) {
            underlying_t aux1 = m_offDiag[rhead(i)+k] / m_diag[k];
            underlying_t aux2 = m_offDiag[rhead(i)+k+1] / m_diag[k];
            m_offDiag[rhead(i)+k]
                = std::fma(aux1, m_offDiag[rhead(k+1)+k], aux2 * m_diag[k-1]);
            m_offDiag[rhead(i)+k+1]
                = std::fma(aux1, m_diag[k+1], (-aux2) * m_offDiag[rhead(k+1)+k]);
        }

        //! Be careful on overflows.
        m_diag[k] =
            fma(m_diag[k-1] / m_diag[k], m_diag[k+1],
                (m_offDiag[rhead(k+1)+k] / m_diag[k]) * m_offDiag[rhead(k+1)+k] );
    }

    //! Update lambda in MINUS step in the HMM algorithm.
    void negate(std::size_t k) noexcept{
        if (k <= 1) return;

        std::slice_array<underlying_t> sl = m_offDiag[std::slice(rhead(k),k,1)];
        sl *= std::valarray<underlying_t>(-1,k);
        for(size_t i = k+1; i < m_size; ++i)
            m_offDiag[rhead(i)+k] *= -1;
    }

    //! Check if SWAP should be performed.
    //! \pre k >= 1 && k <= n
    int should_swap(size_t k) noexcept
    {
        return m_diag[k-1] * m_diag[k+1] + sqpow(m_offDiag[rhead(k+1)+k])
            < delta * sqpow(m_diag[k]);
    }

};

/*!
 * "Swap" operation in the algorithm.
 * \pre 1 <= k < min(data->m.r, data->u.r)
 */
template<class Ops, class Derived, class...Us, class...Vs>
void swap(
    std::size_t k,
    Eigen::MatrixBase<Derived> &m,
    std::tuple<Us&...> &us, std::tuple<Vs&...> &vs,
    Lambda_t &lambda
    ) noexcept
{
    // Transform matrices
    Ops::swap(m, k, k-1);
    for_each_tuple(us, [k](auto& u){ Ops::dual_t::swap(u,k,k-1); });
    for_each_tuple(vs, [k](auto& v){ Ops::swap(v,k,k-1); });

    // Update lambda
    lambda.slice_swap(k+1, k);
    lambda.swap(k);
}

/*!
 * "Minus" operation in the algorithm.
 */
template <class Ops, class Derived, class...Us, class...Vs>
void negate(
    std::size_t k,
    Eigen::MatrixBase<Derived> &m,
    std::tuple<Us&...> &us, std::tuple<Vs&...> &vs,
    Lambda_t &lambda
    ) noexcept
{
    // Negate k-th rows
    Ops::scalar(m, k, -1);

    for_each_tuple(us, [k](auto& u){ Ops::dual_t::scalar(u, k, -1); });
    for_each_tuple(vs, [k](auto& v){ Ops::scalar(v, k, -1); });

    // Update lambda
    lambda.negate(k);
}

/*!
 * "Reduce2" operation in the algorithm.
 * \pre i < k
 * \return flag || (whether "Swap" should be performed or not).
 */
template <class Ops, class Derived, class...Us, class...Vs>
bool reduce(
    std::size_t i, std::size_t k,
    Eigen::MatrixBase<Derived> &m,
    std::tuple<Us&...> &us, std::tuple<Vs&...> &vs,
    Lambda_t &lambda,
    bool flag
    ) noexcept
{
    /* Find the first non-zero entry of the i-th and j-th vectors. */
    std::size_t l1 = Ops::find_nonzero(
        m, i,
        [&](std::size_t l,auto x) {
            if (std::signbit(x)) negate<Ops>(i, m, us, vs, lambda);
        } );
    std::size_t l2 = Ops::find_nonzero(
        m, k,
        [&](std::size_t l,auto x) {
            if (std::signbit(x)) negate<Ops>(k, m, us, vs, lambda);
        } );

    // Scalar factor of the reduction
    typename Eigen::MatrixBase<Derived>::Scalar q;

    // Compute the scalar factor
    if (l1 < Ops::size(m)) {
        q = - floor_div(m(k,l1), m(i,l1));
    }
    else if (2.0*std::abs(lambda(k,i)) > lambda(i,i))
        q = - std::round( lambda(k,i)/ lambda(i,i) );
    else
        q = 0;

    // Reduce the k-th row vector by the i-th one.
    if (q != 0) {
        Ops::axpy(m, q, i, k);
        for_each_tuple(us, [i,k,&q](auto& u){ Ops::dual_t::axpy(u,-q,k,i); });
        for_each_tuple(vs, [i,k,&q](auto& v){ Ops::axpy(v,q,i,k); });

        // Update lambda
        lambda.axpy(q, i, k);
    }

    return flag || (l1 < Ops::size(m) && l1 <= l2)
                || (l1 == Ops::size(m) && l1 == l2 && lambda.should_swap(k));
}


} // end namespace _impl_LLL

} // end namespace khover

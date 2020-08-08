/*!
 * \file hnf.hpp
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
#include "hnf_impl_lll.hpp"

//* For debug
#include <iostream>

#define DBG_MSG std::cerr << __FILE__ << ":" << __LINE__ << std::endl
// */

namespace khover {

/*!
 * Computing the Hermite normal form of a given matrix.
 * The original "pseudo-code" is found in the paper
 *   > George Havas, Bohdan S. Majewski & Keith R. Matthews (1998) Extended GCD and Hermite Normal Form Algorithms via Lattice Basis Reduction, Experimental Mathematics, 7:2, 125-136, DOI: 10.1080/10586458.1998.10504362
 * This function also applies the transformation on the target matrix and its adjoint transformation to given matrices.
 * Usage:
 * \code
 *   // Compute a row echelon form.
 *   auto u0 = u; auto m0 = m; auto v0 = v;
 *   hnf_LLL<typename khover::rowops>(m,std::tie(u),std::tie(v));
 *   assert(is_rowhnf(m));
 *   assert(u*m == u0*m0);
 *   assert(u*v == u0*v0);
 * \endcode
 *
 * \code
 *   // Compute a column echelon form.
 *   auto u0 = u; auto m0 = m; auto v0 = v;
 *   hnf_LLL<typename khover::colops>(m,std::tie(u),std::tie(v));
 *   assert(is_colhnf(m));
 *   assert(m*u == m0*u0);
 *   assert(v*u == v0*u0);
 * \endcode
 * \tparam Ops A collection of elementary operations; \see{khover::rowops}, \see{khover::colops}.
 * \param m The target matrix to be transformed into its Hermite normal form.
 * \param us A tuple of matrices subject to the adjoint transformation.
 * \param vs A tuple of matrices subject to the transformation.
 * \return Always true; (in future) it will be the flag if overflow happened.
 */
template<
    class Ops,
    class MT,int MR, int MC, int MOpt, int MRMax, int MCMax,
    class...UTs, class...VTs
    >
bool hnf_LLL(
    Eigen::Matrix<MT,MR,MC,MOpt,MRMax,MCMax> &m,
    std::tuple<UTs&...> us,
    std::tuple<VTs&...> vs
    ) noexcept
{
    static_assert(
        std::conjunction<typename khover::is_pubbase_of_template<Eigen::MatrixBase,UTs>...>::value,
        "Matrices U contain a class not derived from Eigen::MatrixBase");
    static_assert(
        std::conjunction<typename khover::is_pubbase_of_template<Eigen::MatrixBase,VTs>...>::value,
        "Matices V contain a class not derived from Eigen::MatrixBase");
    static_assert(
        std::conjunction<std::bool_constant<UTs::Flags & Eigen::LvalueBit>...>::value,
        "Matrices U contain read-only variables");
    static_assert(
        std::conjunction<std::bool_constant<VTs::Flags & Eigen::LvalueBit>...>::value,
        "Matrices V contain read-only variables");

    std::size_t nvecs = Ops::dual_t::size(m);

    // If the given matrix consists of a single row vector, then all we have to do is to ensure the first non-zero entry is positive.
    if (nvecs == 1) {
        Ops::find_nonzero(
            m, 0,
            [&](std::size_t l, auto x) {
                if (std::signbit(x)) {
                    Ops::scalar(m,0,-1);
                    for_each_tuple(us, [](auto& u){ Ops::dual_t::scalar(u,0,-1); });
                    for_each_tuple(vs, [](auto& v){ Ops::scalar(v,0,-1); });
                }
            });
        return true;
    }

    _impl_LLL::Lambda_t lambda(nvecs+1);

    // The index of the row that we currently focus on.
    // Note that the row indices may be upside-down.
    size_t cur = 1;

    // Proceed the algorithm on the first k rows.
    while (cur < nvecs) {
        if (_impl_LLL::reduce<Ops>(cur, cur-1, m, us, vs, lambda, false)) {
            _impl_LLL::swap<Ops>(cur, m, us, vs, lambda);
            if (cur > 1) --cur;
        }
        else {
            for (size_t i = 2; i <= cur; ++i)
                _impl_LLL::reduce<Ops>(cur, cur-i, m, us, vs, lambda, true);
            ++cur;
        }
    }

    return true;
}

}

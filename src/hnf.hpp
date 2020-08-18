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
#include "debug/debug.hpp"
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
 * \return If success, the rank of the given matrix over Q (the field of rationals).
 */
template<
    class Ops,
    class MT,int MR, int MC, int MOpt, int MRMax, int MCMax,
    class...UTs, class...VTs
    >
std::optional<std::size_t> hnf_LLL(
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
        std::conjunction<std::bool_constant<(UTs::Flags & Eigen::LvalueBit) != 0>...>::value,
        "Matrices U contain read-only variables");
    static_assert(
        std::conjunction<std::bool_constant<(VTs::Flags & Eigen::LvalueBit) != 0>...>::value,
        "Matrices V contain read-only variables");

    std::size_t nvecs = Ops::dual_t::size(m);

    if (!foldl_tuple(true, us, [nvecs](bool b, auto& u) { return b && Ops::size(u) >= nvecs; })) {
        ERR_MSG("Matricies U with invalid sizes.");
        return std::nullopt;
    }

    if (!foldl_tuple(true, vs, [nvecs](bool b, auto& v) { return b && Ops::dual_t::size(v) >= nvecs; })) {
        ERR_MSG("Matricies V with invalid sizes.");
        return std::nullopt;
    }

    // Nothing to do on empty matrices
    if (nvecs == 0 || Ops::size(m) == 0) {
        return std::make_optional(0);
    }

    // Ensure the pivot of the last vector to be non-negative.
    std::size_t l = Ops::find_nonzero(
        m, nvecs-1,
        [&m,&us,&vs](std::size_t l, auto x) {
            if (std::signbit(x)) {
                Ops::scalar(m,0,-1);
                for_each_tuple(us, [](auto& u){ Ops::dual_t::scalar(u,0,-1); });
                for_each_tuple(vs, [](auto& v){ Ops::scalar(v,0,-1); });
            }
        });

    // If the given matrix consists of a single vector, then all the step is finished.
    if (nvecs == 1) {
        return l < Ops::size(m) ? 1 : 0;
    }

    _impl_LLL::Lambda_t lambda(nvecs);

    // The index of the vector that we currently focus on.
    std::size_t cur = nvecs - 1;
    // The rank of the span of vectors below cursor.
    std::size_t rk = 0;
    // Flag whether the vector just below the cursor is non-zero or not.
    bool is_below_nz = false;

    // Proceed the algorithm on the first k rows.
    while (cur > 0) {
        //DBG_MSG("cur=" << cur << "\n" << "rk = " << rk << "\n" << m);

        auto howswap = _impl_LLL::reduce<Ops,false>(cur-1, cur, m, us, vs, lambda);

        if (howswap & _impl_LLL::HowSwap::ShouldSwap) {
            _impl_LLL::swap<Ops>(cur-1, m, us, vs, lambda);
            if (cur+1 < nvecs) {
                ++cur;
                if (is_below_nz) {
                    --rk;
                    is_below_nz = rk > 0;
                }
            }
        }
        else if (howswap & _impl_LLL::HowSwap::ZeroReducer) {
            --cur;
            is_below_nz = false;
        }
        else {
            for (std::size_t i = cur+1; i < nvecs; ++i)
                _impl_LLL::reduce<Ops,true>(cur-1, i, m, us, vs, lambda);
            --cur;
            ++rk;
            is_below_nz = true;
        }
    }

    if (rk > 0) {
        return rk+1;
    }
    else {
        return Ops::find_nonzero(m, 0, [](auto,auto){})
            < Ops::size(m)
              ? 1 : 0;
    }
}

}

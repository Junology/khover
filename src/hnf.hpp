/*!
 * \file hnf.hpp
 * \author Jun Yoshida
 * \copyright (c) 2020 Jun Yoshida.
 * The project is released under the 2-clause BSD License.
 * \date August, 2020: created
 */

#pragma once

#include <Eigen/Dense>

#include "utils.hpp"

namespace khover {

/*!
 * Computing the Hermite normal form of a given matrix.
 * \code
 *   auto u0 = u; auto m0 = m; auto v0 = v;
 *   hnf_row(u,m,v);
 *   assert(is_hnf(m));
 *   assert(u*m == u0*m0);
 *   assert(u*v == u0*v0);
 * \endcode
 */
template<
    class MT,int MR, int MC, int MOpt, int MRMax, int MCMax,
    class UDer, class... VTs
    >
bool hnf_row(Eigen::MatrixBase<UDer> &u,
             Eigen::Matrix<MT,MR,MC,MOpt,MRMax,MCMax> &m,
             VTs&...vs) noexcept
{
    static_assert(
        std::conjunction<typename khover::is_pubbase_of_template<Eigen::MatrixBase,VTs>::type...>::value,
        "Matices V are not derived from Eigen::MatrixBase");
    return true;
}

}

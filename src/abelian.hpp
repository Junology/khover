/*!
 * \file complex.hpp
 * \author Jun Yoshida
 * \copyright (c) 2020 Jun Yoshida.
 * The project is released under the 2-clause BSD License.
 * \date August, 2020: created
 */

#pragma once

#include <cstdint>
#include <array>
#include <queue>

#include <Eigen/Dense>

#include "hnf.hpp"

/* Debug
#include <iostream>
// */

namespace khover {

class AbelianGroup {
public:
    using integer_t = int64_t;
    using matrix_t = Eigen::Matrix<integer_t,Eigen::Dynamic,Eigen::Dynamic>;

private:
    //! The presentation matrix, which is kept to be a column echelon form.
    matrix_t m_repMat;
    //! The least rank of the free part of the abelian group.
    std::size_t m_freerk;

public:
    //! Default constructor.
    AbelianGroup() noexcept : m_repMat(0,0), m_freerk(0) {}

    //! Constructor from a presentation matrix.
    //! Represent an abelian group by a presentation matrix.
    //! \tparam IsColEchelon If IsColEchelon::value == true, then repMat must be in the column echelon form. Otherwise, its column HNF is computed.
    template <class Derived, class IsColEchelon>
    AbelianGroup(Eigen::MatrixBase<Derived> const& repMat, IsColEchelon, std::size_t freerk = 0) noexcept
        : m_repMat(repMat), m_freerk(freerk)
    {
        if constexpr(!IsColEchelon::value) {
            auto rk = hnf_LLL<khover::colops>(m_repMat, {}, {});
            if(rk) {
                m_repMat = m_repMat.leftCols(*rk).eval();
            }
            else {
                DBG_MSG("Something bad happended.");
            }
        }
    }

    //! Nothing special to do in destructor.
    ~AbelianGroup() = default;

    //! Get the number of generators.
    inline std::size_t ngens() const noexcept {
        return m_repMat.rows() + m_freerk;
    }

    //! Get the number of relations.
    inline std::size_t nrels() const noexcept { return m_repMat.cols(); }

    //! Get the presentation matrix.
    inline auto get_repmatrix() const noexcept {
        return (matrix_t(ngens(),nrels()) << m_repMat, matrix_t::Zero(m_freerk, m_repMat.cols())).finished();
    }

    //! Reduce the number of generators and relations by computing the Hermite normal form of the representation matrix.
    //! This will suffice in order to compute the rank of the free part.
    //! \param post Homomorphisms whose domains are *this*. They will be transformed so that their domains will be the reduced one.
    //! \param pre Homomorphisms whose codomains are *this*. They will be transformed so that their codomains will be the reduced one.
    //! \return If the function fails to compute HNF correctly, or if the homomorphisms are wrong, then it returns false.
    template <class...Posts, class...Pres>
    bool reduce(
        std::tuple<Posts&...> const & posts,
        std::tuple<Pres&...> const & pres
        ) noexcept
    {
        /*
         * Take columns with pivot = 1 to left
         */

        // Queue to store the indices of columns with non-unit pivots.
        std::queue<std::size_t> q_nupiv;
        // The rank of the presentation matrix over Q (the field of rationals).
        std::size_t rk = 0;
        // The number of columns with pivot = 1.
        std::size_t upivs = 0;

        // Traverse pivots
        using Ops = khover::colops;
        for(std::size_t i=0; i < Ops::size(m_repMat) && rk < Ops::dual_t::size(m_repMat); ++i) {
            // A column with pivot 1 is found.
            if (Ops::at(m_repMat, rk).coeff(i) == 1) {
                // If there is a column with non-unit pivot on left, swap with it.
                if (!q_nupiv.empty()) {
                    Ops::swap(m_repMat, rk, q_nupiv.front());
                    q_nupiv.pop();
                    q_nupiv.push(rk);
                }
                ++upivs;
                ++rk;
            }
            // A column with pivot != 1 is found.
            else if (Ops::at(m_repMat, rk).coeff(i) != 0) {
                q_nupiv.push(rk);
                ++rk;
            }
        }

        // Compute the row HNF
        if (!hnf_LLL<khover::rowops>(m_repMat,posts, pres))
            return false;

        // Reduce the presentation.
        m_freerk += static_cast<std::size_t>(m_repMat.rows() - rk);

        std::size_t rs = std::min(static_cast<std::size_t>(m_repMat.rows()),rk)-upivs;
        std::size_t cs = rk - upivs;
        m_repMat = m_repMat.block(upivs, upivs, rs, cs);

        // Reduce the homomorphisms.
        for_each_tuple(
            pres,
            [upivs](auto& u){ u = u.bottomRows(u.rows() - upivs); });
        for_each_tuple(
            posts,
            [upivs](auto& v){ v = v.rightCols(v.cols() - upivs); });

        // Keep the presentation matrix in column echelon forms.
        khover::hnf_LLL<khover::colops>(m_repMat, {}, {});

        // Finish successfully.
        return true;
    }
};

} // end namespace khover

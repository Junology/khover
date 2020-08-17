/*!
 * \file chaincomplex.hpp
 * \author Jun Yoshida
 * \copyright (c) 2020 Jun Yoshida.
 * The project is released under the 2-clause BSD License.
 * \date August, 2020: created
 */

#pragma once

#include <cstdint>
#include <deque>
#include <vector>

#include "abelian.hpp"

namespace khover{

template <class Coeff, class Homology>
struct ChainCoeffs {
    using coeff_t = Coeff;
    using homology_t = Homology;
};

/*!
 * Chain complexes over the ring of integers.
 * Homological degrees are used.
 * \warning{
 *   The class cannot represent complexes with single non-zero terms.
 *   A workaournd is formally appending the matrix of size 0 x rank.
 * }
 */
class ChainIntegral
    : public ChainCoeffs<std::int64_t,AbelianGroup>
{
public:
    using matrix_t = homology_t::matrix_t;

private:
    //! The minimal (homological) degree in which the complex is non-zero.
    int m_mindeg;

    //! The list of differentials.
    //! The compute() member function assumes this is not empty.
    std::deque<matrix_t> m_diffs;

public:
    ChainIntegral() = default;

    template<class Derived, class...Deriveds>
    ChainIntegral(
        int mindeg,
        Eigen::MatrixBase<Derived> const& diffhead,
        Eigen::MatrixBase<Deriveds> const&... diffs)
        noexcept
        : m_mindeg(mindeg), m_diffs{diffhead, diffs...}
    {}

    //! The length of the complex.
    inline std::size_t length() noexcept { return m_diffs.size() + 1; }

    inline int mindeg() noexcept { return m_mindeg; }
    inline int maxdeg() noexcept { return m_mindeg + m_diffs.size(); }

    //! Append a differential at the highest (homological) degree.
    //! Composability check is performed (if failed, the function does nothing).
    template<class Derived>
    bool prepend(Eigen::MatrixBase<Derived> const& diff) noexcept {
        if (diff.rows() == m_diffs.back().cols()) {
            m_diffs.push_back(diff);
            return true;
        }
        else {
            DBG_MSG(diff.rows() << "!=" << m_diffs.back().cols());
            return false;
        }
    }

    //! Append a differential at the highest (homological) degree.
    //! Composability check is performed (if failed, the function does nothing).
    template<class Derived>
    bool prepend(Eigen::MatrixBase<Derived> && diff) noexcept {
        if (diff.rows() == m_diffs.back().cols()) {
            m_diffs.push_back(std::move(diff));
            return true;
        }
        else {
            DBG_MSG(diff.rows() << "!=" << m_diffs.back().cols());
            return false;
        }
    }

    //! Append a differential at the lowerst (homological) degree
    //! Composability check is performed (if failed, the function does nothing).
    template<class Derived>
    bool append(Eigen::MatrixBase<Derived> const& diff) noexcept {
        if (diff.cols() == m_diffs.front().rows()) {
            m_diffs.push_front(diff);
            --m_mindeg;
            return true;
        }
        else {
            DBG_MSG(diff.rows() << "!=" << m_diffs.back().cols());
            return false;
        }
    }

    //! Append a differential at the lowerst (homological) degree
    //! Composability check is performed (if failed, the function does nothing).
    template<class Derived>
    bool append(Eigen::MatrixBase<Derived> && diff) noexcept {
        if (diff.cols() == m_diffs.front().rows()) {
            m_diffs.push_front(std::move(diff));
            --m_mindeg;
            return true;
        }
        else {
            DBG_MSG(diff.rows() << "!=" << m_diffs.back().cols());
            return false;
        }
    }

    /*!
     * Compute the homology groups.
     * \code
     *   ChainIntegral C(...);
     *   // ...
     *   auto result = C.compute;
     *   // -> H_n(C) = result[n+C.mindeg()];
     * \endcode
     */
    std::vector<homology_t> compute() {
        std::vector<matrix_t> diff_mutable(m_diffs.begin(), m_diffs.end());
        std::vector<homology_t> result{};

        for(std::size_t i = 0; i < m_diffs.size()-1; ++i) {
            auto rk = hnf_LLL<khover::colops>(
                diff_mutable[i], std::tie(diff_mutable[i+1]), {});
            if (!rk) {
                DBG_MSG("Something bad happended.");
                return {};
            }
            diff_mutable[i] = diff_mutable[i].leftCols(*rk).eval();
            diff_mutable[i+1] = diff_mutable[i+1].bottomRows(
                diff_mutable[i+1].rows() - *rk).eval();
            result.emplace_back(diff_mutable[i], std::true_type{});
        }

        //! Compute the highest homology group.
        auto rk = hnf_LLL<khover::colops>(diff_mutable.back(), {}, {});
        if (!rk) {
            DBG_MSG("Something bad happended.");
            return {};
        }

        result.emplace_back(diff_mutable.back(), std::true_type{});
        result.emplace_back(
            matrix_t(0,0), std::true_type{},
            static_cast<std::size_t>(diff_mutable.back().cols()) - *rk);

        return result;
    }
};

} // end namespace khover

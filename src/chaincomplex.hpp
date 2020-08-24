/*!
 * \file chaincomplex.hpp
 * \author Jun Yoshida
 * \copyright (c) 2020 Jun Yoshida.
 * The project is released under the 2-clause BSD License.
 * \date August, 2020: created
 */

#pragma once

#include <cstdint>
#include <limits>
#include <deque>
#include <vector>

#include "abelian.hpp"

namespace khover{

//! A set of static const variables.
template <class T>
struct constants {
static inline auto const empty_matrix
  = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>(0,0);
};

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

    //! The class for maps on graded objects.
    //! It is mainly, while not limited to, used to represent chain maps.
    class Hom {
        friend class ChainIntegral;

    private:
        //! The data for maps.
        //! - The first component carries the minimum homological degree where is non-zero map.
        //! - The second component carries matrix representations of maps.
        //! For out-of-range positions, we assume there are the zero-maps.
        std::optional<std::pair<int,std::vector<matrix_t>>> m_morphs;

        Hom() noexcept = default;

    public:
        template <class MorT>
        Hom(int mindeg, MorT const& morphs) noexcept
            : m_morphs(
                std::in_place,
                mindeg,
                std::vector<matrix_t>(std::begin(morphs), std::end(morphs)))
        {
        }

        Hom(int mindeg, std::vector<matrix_t> &&morphs) noexcept
            : m_morphs(
                std::in_place_t{},
                mindeg,
                std::move(morphs))
        {}

        Hom(Hom const&) = default;
        Hom(Hom &&) noexcept = default;

        Hom& operator=(Hom const&) = default;
        Hom& operator=(Hom &&) noexcept = default;

        inline int mindeg() const noexcept {
            return m_morphs
                ? m_morphs->first
                : std::numeric_limits<int>::max();
        }
        inline int maxdeg() const noexcept {
            return m_morphs
                ? m_morphs->first + static_cast<int>(m_morphs->second.size())-1
                : std::numeric_limits<int>::min();
        }

        inline matrix_t const& morphism(int i) const noexcept {
            return m_morphs && i-mindeg() >= 0 && i-mindeg() <= maxdeg()
                ? m_morphs->second[i-mindeg()]
                : constants<coeff_t>::empty_matrix;
        }

        inline std::pair<std::size_t,std::size_t> ranks(int i) const noexcept {
            return std::make_pair(morphism(i).cols(), morphism(i).rows());
        }

        //! Precompose another hom.
        Hom& operator*=(Hom const& other) noexcept {
            // Trivial cases
            if(!m_morphs && !other.m_morphs)
                return *this;
            else if(!m_morphs) {
                m_morphs.emplace(
                    other.m_morphs->first,
                    std::vector<matrix_t>{});
                std::transform(
                    std::begin(other.m_morphs->second),
                    std::end(other.m_morphs->second),
                    std::back_inserter(m_morphs->second),
                    [](matrix_t const& mor) {
                        return matrix_t(0, mor.cols());
                    });
                return *this;
            }
            else if(!other.m_morphs) {
                std::for_each(
                    std::begin(m_morphs->second),
                    std::end(m_morphs->second),
                    [](matrix_t& mor) {
                        mor.resize(Eigen::NoChange, 0);
                    });
                return *this;
            }

            int newmin = std::min(mindeg(), other.mindeg());
            int newmax = std::max(maxdeg(), other.maxdeg());

            // Compute the composition in each degree.
            m_morphs.emplace(newmin, std::vector<matrix_t>{});
            m_morphs->second.reserve(
                static_cast<std::size_t>(std::max(0, newmax - newmin)));
            for(int i = newmin; i <= newmax; ++i) {
                m_morphs->second.push_back(
                    morphism(i)*other.morphism(i));
            }

            return *this;
        }
    };

private:
    //! The minimal (homological) degree in which the complex is non-zero.
    int m_mindeg;

    //! The list of differentials.
    //! The compute() member function assumes this is not empty.
    std::deque<matrix_t> m_diffs;

    //! The trivial differentials in the ends of the complex.
    matrix_t m_head_zero, m_tail_zero;

public:
    ChainIntegral() = default;

    template<class Derived, class...Deriveds>
    ChainIntegral(
        int mindeg,
        Eigen::MatrixBase<Derived> const& diffhead,
        Eigen::MatrixBase<Deriveds> const&... diffs)
        noexcept
        : m_mindeg(mindeg), m_diffs{diffhead, diffs...},
          m_head_zero(0, m_diffs.front().rows()),
          m_tail_zero(m_diffs.back().cols(), 0)
    {}

    ChainIntegral(
        int mindeg,
        decltype(m_diffs) && diffs)
        noexcept
        : m_mindeg(mindeg), m_diffs(std::move(diffs)),
          m_head_zero(0, m_diffs.empty() ? 0 : m_diffs.front().rows()),
          m_tail_zero(m_diffs.empty() ? 0 : m_diffs.back().cols(), 0)
    {}

    //! The length of the complex.
    inline std::size_t length() const noexcept { return m_diffs.size() + 1; }

    inline int mindeg() const noexcept {
        return m_diffs.empty()
            ? std::numeric_limits<int>::max() : m_mindeg;
    }
    inline int maxdeg() const noexcept {
        return m_diffs.empty()
            ? std::numeric_limits<int>::min() : m_mindeg + m_diffs.size();
    }

    //! Get the rank of the complex in degree *i*.
    std::size_t rank(int i) const noexcept {
        if (!m_diffs.empty() && i==maxdeg()) {
            return m_diffs.back().cols();
        }
        else if(i >= mindeg() && i < maxdeg()) {
            return m_diffs[i-mindeg()].rows();
        }
        else {
            return 0;
        }
    }

    // Get differential as a matri.
    matrix_t const& getDiff(int i) const noexcept {
        if (i < mindeg()-1 || i > maxdeg())
            return constants<coeff_t>::empty_matrix;
        else if (i == mindeg()-1)
            return m_head_zero;
        else if (i == maxdeg())
            return m_tail_zero;
        else
            return m_diffs[i-m_mindeg];
    }

    //! Append a differential at the highest (homological) degree.
    //! Composability check is performed (if failed, the function does nothing).
    template<class Derived>
    bool prepend(Eigen::MatrixBase<Derived> const& diff) noexcept {
        if (diff.rows() == m_diffs.back().cols()) {
            m_tail_zero.resize(diff.cols(), 0);
            m_diffs.push_back(diff);
            return true;
        }
        else {
            ERR_MSG(diff.rows() << "!=" << m_diffs.back().cols());
            return false;
        }
    }

    //! Append a differential at the highest (homological) degree.
    //! Composability check is performed (if failed, the function does nothing).
    template<class Derived>
    bool prepend(Eigen::MatrixBase<Derived> && diff) noexcept {
        if (diff.rows() == m_diffs.back().cols()) {
            m_tail_zero.resize(diff.cols(), 0);
            m_diffs.push_back(std::move(diff));
            return true;
        }
        else {
            ERR_MSG(diff.rows() << "!=" << m_diffs.back().cols());
            return false;
        }
    }

    //! Append a differential at the lowerst (homological) degree
    //! Composability check is performed (if failed, the function does nothing).
    template<class Derived>
    bool append(Eigen::MatrixBase<Derived> const& diff) noexcept {
        if (diff.cols() == m_diffs.front().rows()) {
            m_head_zero.resize(0, diff.rows());
            m_diffs.push_front(diff);
            --m_mindeg;
            return true;
        }
        else {
            ERR_MSG(diff.rows() << "!=" << m_diffs.back().cols());
            return false;
        }
    }

    //! Append a differential at the lowerst (homological) degree
    //! Composability check is performed (if failed, the function does nothing).
    template<class Derived>
    bool append(Eigen::MatrixBase<Derived> && diff) noexcept {
        if (diff.cols() == m_diffs.front().rows()) {
            m_head_zero.resize(0, diff.rows());
            m_diffs.push_front(std::move(diff));
            --m_mindeg;
            return true;
        }
        else {
            ERR_MSG(diff.rows() << "!=" << m_diffs.back().cols());
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
    std::vector<homology_t> compute() const noexcept {
        std::vector<matrix_t> diff_mutable(m_diffs.begin(), m_diffs.end());
        std::vector<homology_t> result{};

        for(std::size_t i = 0; i < m_diffs.size()-1; ++i) {
            auto rk = hnf_LLL<khover::colops>(
                diff_mutable[i], std::tie(diff_mutable[i+1]), {});
            if (!rk) {
                ERR_MSG("Failed to decompose the chain complex.");
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
            ERR_MSG("Failed to decompose the chain complex.");
            return {};
        }

        result.emplace_back(diff_mutable.back(), std::true_type{});
        result.emplace_back(
            matrix_t(0,0), std::true_type{},
            static_cast<std::size_t>(diff_mutable.back().cols()) - *rk);

        return result;
    }

    /*!
     * Taking the mapping cone
     */
    static
    std::optional<ChainIntegral> cone(
        Hom const& hom,
        ChainIntegral const& dom,
        ChainIntegral const& cod
        ) noexcept
    {
        decltype(m_diffs) diffs{};
        int mindeg = std::min(dom.mindeg()+1, cod.mindeg());
        int maxdeg = std::max(dom.maxdeg()+1, cod.maxdeg());

        for(int i = mindeg; i < maxdeg; ++i) {
            // Compatibility check
            if(auto [domrk, codrk] = hom.ranks(i);
               domrk != dom.rank(i) || codrk != cod.rank(i))
            {
                ERR_MSG(
                    "Given homomorphism has incompatible (co)domains in degree "
                    << i << ":\n"
                    << domrk << "!=" << dom.rank(i) << " || "
                    << codrk << "!=" << cod.rank(i));
                return std::nullopt;
            }

            diffs.push_back(
                matrix_t::Zero(
                    cod.rank(i) + dom.rank(i-1),
                    cod.rank(i+1) + dom.rank(i)
                    ));
            diffs.back().topLeftCorner(cod.rank(i), cod.rank(i+1))
                = cod.getDiff(i);
            diffs.back().bottomRightCorner(dom.rank(i-1), dom.rank(i))
                = - dom.getDiff(i-1);
            diffs.back().topRightCorner(cod.rank(i), dom.rank(i))
                = hom.morphism(i);
        }

        return std::make_optional<ChainIntegral>(
            mindeg, std::move(diffs));
    }
};

} // end namespace khover

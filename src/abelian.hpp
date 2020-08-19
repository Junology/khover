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

//* Debug
#include "debug/debug.hpp"
//*/

namespace khover {

//! Representation of abelian groups by direct sums of cyclic groups.
struct AbGroupCyc {
    //! Rank of the free part.
    std::size_t freerank;
    //! A list of torsions.
    std::vector<int> torsions;

    //! Pretty printer
    std::string pretty() const noexcept {
        std::string str_free
            = freerank == 0 ? std::string{} : freerank == 1 ? std::string("Z") : "Z^" + std::to_string(freerank);
        std::string str_tor{};

        for(int t : torsions) {
            if (!str_tor.empty())
                str_tor += "+";
            str_tor += "Z/" + std::to_string(t);
        }

        if(str_free.empty() && str_tor.empty())
            return std::string("0");
        else if (str_free.empty())
            return str_tor;
        else if (str_tor.empty())
            return str_free;
        else
            return str_free + "+" + str_tor;
    }
};

//! Representation of abelian groups by representation matrices.
class AbelianGroup {
public:
    using integer_t = int64_t;
    using matrix_t = Eigen::Matrix<integer_t,Eigen::Dynamic,Eigen::Dynamic>;

private:
    //! The presentation matrix, which is kept to be a column HNF.
    matrix_t m_repMat;
    //! The least rank of the free part of the abelian group.
    std::size_t m_freerk;

public:
    //! Default constructor.
    AbelianGroup() noexcept : m_repMat(0,0), m_freerk(0) {}

    //! Constructor from a presentation matrix.
    //! Represent an abelian group by a presentation matrix.
    //! \tparam IsColHNF If IsColHNF::value == true, then repMat must be in the column HNF. Otherwise, its column HNF is computed.
    template <class Derived, class IsColHNF>
    AbelianGroup(Eigen::MatrixBase<Derived> const& repMat, IsColHNF, std::size_t freerk = 0) noexcept
        : m_repMat(repMat), m_freerk(freerk)
    {
        if constexpr(!IsColHNF::value) {
            auto rk = hnf_LLL<khover::colops>(m_repMat, {}, {});
            if(rk) {
                m_repMat = m_repMat.leftCols(*rk).eval();
            }
            else {
                ERR_MSG("Something bad happended.");
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
        static_assert(
            std::conjunction_v<
            std::bool_constant<Posts::ColsAtCompileTime == Eigen::Dynamic>...
            >,
            "The function may change the number of columns of matrices in the posts parameter, so make sure they have dynamic numbers of columns.");

        static_assert(
            std::conjunction_v<
            std::bool_constant<Pres::RowsAtCompileTime == Eigen::Dynamic>...
            >,
            "The function may change the number of rows of matrices in the posts parameter, so make sure they have dynamic numbers of rows.");

        /*
         * Take columns with pivot = 1 to left
         */
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
                if (rk > upivs) {
                    Ops::swap(m_repMat, rk, upivs);
                }
                // If the pivot is below diagonal, raise it.
                if (i > upivs) {
                    Ops::dual_t::swap(m_repMat, i, upivs);
                    khover::for_each_tuple(
                        pres,
                        [i,upivs](auto& u) {
                            Ops::dual_t::swap(u, i, upivs);
                        });
                    khover::for_each_tuple(
                        posts,
                        [i,upivs](auto& v) {
                            Ops::swap(v, i, upivs);
                        });
                }
                ++upivs;
                ++rk;
            }
            // A column with pivot != 1 is found.
            else if (Ops::at(m_repMat, rk).coeff(i) != 0) {
                ++rk;
            }
        }

        // Reduce homomorphisms
        matrix_t redPiv_mat
            = matrix_t::Identity(
                m_repMat.rows() + m_freerk,
                m_repMat.rows() + m_freerk
                ).bottomRows(m_repMat.rows() + m_freerk - upivs);
        redPiv_mat.block(0, 0, m_repMat.rows() - upivs, upivs).noalias()
            = - m_repMat.bottomLeftCorner(m_repMat.rows() - upivs, upivs);

        khover::for_each_tuple(
            pres,
            [&](auto& u) {
                u = redPiv_mat.bottomRows(m_repMat.rows() + m_freerk - upivs) * u;
            } );
        khover::for_each_tuple(
            posts,
            [upivs](auto& v) {
                v = v.rightCols(v.cols() - upivs); //v.block(0, upivs, v.rows(), v.cols()-upivs).eval();
            } );

        // Forget thw row/columns with unital pivots.
        m_freerk += m_repMat.rows() - rk;
        m_repMat = m_repMat.block(upivs, upivs, m_repMat.rows() - upivs, rk - upivs).eval();
        rk -= upivs;

        // Compute the row HNF
        if (!hnf_LLL<khover::rowops>(m_repMat,posts, pres))
            return false;

        // Reduce the presentation matrix again.
        m_repMat = m_repMat.topRows(rk).eval();

        // Keep the presentation matrix in column HNF.
        khover::hnf_LLL<khover::colops>(m_repMat, {}, {});

        // Finish successfully.
        return true;
    }

    //! Compute the abelian group in the form of the pair of the free rank and the list of torsions.
    //! \remark The result is not necessarily in the normal form.
    AbGroupCyc
    compute() noexcept {
        do {
            reduce({}, {});
        } while(!m_repMat.isDiagonal());

        auto diag = m_repMat.diagonal();
        std::vector<int> torsions{};
        std::size_t r = 0;

        for(int i = 0; i < diag.size(); ++i) {
            if (diag.coeff(i) == 0)
                ++r;
            else if(std::abs(diag.coeff(i)) != 1)
                torsions.push_back(std::abs(diag.coeff(i)));
        }

        return {r+m_freerk, std::move(torsions)};
    }

    //! Compute the image of a homomorphism whose codomain is this abelian group.
    //! \tparam ReturnMorphism If ReturnMorphism::value == true, then the function also returns the matrix representing the morphism from the image to this group.
    template <class Derived, class DoesReturnMorphism = std::false_type>
    auto
    image(
        Eigen::MatrixBase<Derived> const& morph,
        DoesReturnMorphism = DoesReturnMorphism{}
        ) const noexcept
        -> std::conditional_t<
            DoesReturnMorphism::value,
            std::optional<std::pair<AbelianGroup,matrix_t>>,
            std::optional<AbelianGroup>
            >
    {
        // If the number of rows doesn't agree, return immediately.
        if(morph.rows() != static_cast<int>(m_freerk) + m_repMat.rows())
            return std::nullopt;

        // Sum of the images of the homomorphism and the representation matrix.
        matrix_t sumspace(morph.rows(), morph.cols() + m_repMat.cols());
        sumspace << morph, get_repmatrix();

        // Compute the column HNF of the sumspace.
        matrix_t u = matrix_t::Identity(sumspace.cols(), sumspace.cols());
        auto rk = hnf_LLL<khover::colops>(sumspace, std::tie(u), {});

        // Error occured.
        if(!rk)
            return std::nullopt;

        // Return the resulting abelian group together with the homomorphism from it.
        if constexpr (DoesReturnMorphism::value) {
            return std::make_pair(
                AbelianGroup(
                    u.topRightCorner(*rk, m_repMat.cols()),
                    std::false_type{}),
                sumspace.leftCols(*rk)
                );
        }
        else {
            return AbelianGroup(
                u.topRightCorner(*rk, m_repMat.cols()),
                std::false_type{});
        }
    }
};

} // end namespace khover

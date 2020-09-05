/*!
 * \file states.hpp
 * \author Jun Yoshida
 * \copyright (c) 2020 Jun Yoshida.
 * The project is released under the 2-clause BSD License.
 * \date September, 2020: created
 */

#pragma once

#include "chaincomplex.hpp"
#include "cubes.hpp"

namespace khover {

//! The type carrying the data of enhancements on each state.
struct EnhancementProperty {
    //! The index of the first enhancement on a state among a certain (co)homological degree.
    ChainIntegral::matrix_t::Index headidx;

    //! The number of 'x' in the enhancements.
    //! This might be negative or exceed the number of components in case no enhancement is allowed in the q-degree.
    int xcnt = -1;
};


//! Generate the list of enhancement properties for a given state cube.
template <template <class...> class C>
auto get_enhancement_prop(
    LinkDiagram const& diagram,
    Cube<C> const& cube,
    int qdeg) noexcept
    -> std::optional<std::vector<EnhancementProperty>>
{
    std::vector<EnhancementProperty> result{};
    result.reserve(cube.size());
    for(std::size_t st = 0; st < cube.size(); ++st) {
        ChainIntegral::matrix_t::Index headidx = 0;

        if(bitsIndex(state_t{st}) > 0) {
            auto st_pred = predWithPop(state_t{st}).to_ulong();
            headidx = result[st_pred].headidx
                + binom(cube[st_pred].ncomp, result[st_pred].xcnt);
        }

        // In case where component data is not provided
        if (cube[st].arccomp.empty()) {
            result.push_back({headidx, -1});
            continue;
        }

        // Compute the degree.
        int deg = qdeg - diagram.writhe()
            - diagram.cohDegree(cube[st].state);

        // Parity check
        if ((cube[st].ncomp - deg)%2 != 0)
            return std::nullopt;

        // Append the result
        result.push_back(
            {headidx, (static_cast<int>(cube[st].ncomp) - deg)/2});
    }

    return std::make_optional(std::move(result));
}


/*!
 * \section Operations on enhancements.
 */

//! Generate matrices for multiplications
//! By the definition of Khovanov complexes, the resulting matrices consist of 0 and 1.
//! \param numx The number of 'x's in enhancements on the codomain state. It turns out that the domain has to have the same number of 'x's.
//! \param ncomps_in_cod The number of connected components in the codomain state. Since this matrix represents the multiplication, the domain has to have exactly (ncomps_in_cod + 1) components.
//! \param comp_idx1 The common index of the components in the domain and the codomain states which are involved with the multiplication.
//! \param comp_idx2 The index of the other component in the domain involved with the multiplication.
//! \remark The function assumes idx1 < idx2.
//! \return The table of positions in the matrix where coefficients are 1.
auto matrix_mult(
    int numx, int ncomps_in_cod,
    std::size_t comp_idx1, std::size_t comp_idx2)
    noexcept
    -> std::vector<
        std::pair<
            ChainIntegral::matrix_t::Index,
            ChainIntegral::matrix_t::Index
            >
        >;

//! Generate matrices for comultiplications
//! By the definition of Khovanov complexes, the resulting matrices consist of 0 and 1.
//! \param numx The number of 'x's in enhancements on the codomain state. It turns out that the domain has to have exactly (numx-1) 'x's.
//! \param ncomps_in_cod The number of connected components in the codomain state. Since this matrix represents the comultiplication, the domain has to have exactly (ncomps_in_cod - 1) components.
//! \param comp_idx1 The common index of the components in the domain and the codomain states which are involved with the comultiplication.
//! \param comp_idx2 The index of the other component in the codomain involved with the multiplication.
//! \remark The function assumes idx1 < idx2.
//! \return The table of positions in the matrix where coefficients are 1.
auto matrix_comult(
    int numx, int ncomps_in_cod,
    std::size_t comp_idx1, std::size_t comp_idx2)
    noexcept
    -> std::vector<
        std::pair<
            ChainIntegral::matrix_t::Index,
            ChainIntegral::matrix_t::Index
            >
        >;

} // end namespace khover

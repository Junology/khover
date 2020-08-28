/*!
 * \file khovanov.cpp
 * \author Jun Yoshida
 * \copyright (c) 2020 Jun Yoshida.
 * The project is released under the 2-clause BSD License.
 * \date August, 2020: created
 */

#include <vector>

#include <Eigen/Dense>

#include "khovanov.hpp"
#include "cubes.hpp"

/* Debug
#include "debug/debug.hpp"
//*/

using namespace khover;

using matrix_t = ChainIntegral::matrix_t;


/******************************
 *** Enhancement properties ***
 ******************************/
//! The type carrying the data of enhancements on each state.
struct EnhancementProperty {
    //! The index of the first enhancement on a state among a certain (co)homological degree.
    matrix_t::Index headidx;

    //! The number of 'x' in the enhancements.
    //! This might be negative or exceed the number of components in case no enhancement is allowed in the q-degree.
    int xcnt = -1;
};

template <class T>
static
std::optional<std::vector<EnhancementProperty>>
get_enhancement_prop(
    LinkDiagram const& diagram,
    //std::vector<T> const& smoothData,
    T const& smoothData,
    int qdeg)
{
    /*
    static_assert(std::is_base_of_v<LinkDiagram::smoothdata_t, T>,
                  "The template parameter T must be derived from LinkDiagram::smoothdata_t.");
    */

    std::vector<EnhancementProperty> result{};
    result.reserve(smoothData.size());
    for(std::size_t st = 0; st < smoothData.size(); ++st) {
        matrix_t::Index headidx = 0;

        if(bitsIndex(state_t{st}) > 0) {
            auto st_pred = predWithPop(state_t{st}).to_ulong();
            headidx = result[st_pred].headidx
                + binom(smoothData[st_pred].ncomp, result[st_pred].xcnt);
        }

        // In case where component data is not provided
        if (smoothData[st].arccomp.empty()) {
            result.push_back({headidx, -1});
            continue;
        }

        // Compute the degree.
        int deg = qdeg - diagram.writhe()
            - diagram.cohDegree(smoothData[st].state);

        // Parity check
        if ((smoothData[st].ncomp - deg)%2 != 0)
            return std::nullopt;

        // Append the result
        result.push_back(
            {headidx, (static_cast<int>(smoothData[st].ncomp) - deg)/2});
    }

    return std::make_optional(std::move(result));
}

/*******************************************
 *** Multiplication and comultiplication ***
 *******************************************/

//! Generate matrices for multiplications
//! By the definition of Khovanov complexes, the resulting matrices consist of 0 and 1.
//! \param numx The number of 'x's in enhancements on the codomain state. It turns out that the domain has to have the same number of 'x's.
//! \param ncomps_in_cod The number of connected components in the codomain state. Since this matrix represents the multiplication, the domain has to have exactly (ncomps_in_cod + 1) components.
//! \param comp_idx1 The common index of the components in the domain and the codomain states which are involved with the multiplication.
//! \param comp_idx2 The index of the other component in the domain involved with the multiplication.
//! \remark The function assumes idx1 < idx2.
//! \return The table of positions in the matrix where coefficients are 1.
static
std::vector<std::pair<matrix_t::Index,matrix_t::Index>>
matrix_mult(
    int numx, int ncomps_in_cod,
    std::size_t comp_idx1, std::size_t comp_idx2)
    noexcept
{
    auto enh_bound = binom(ncomps_in_cod, numx);
    std::vector<std::pair<matrix_t::Index,matrix_t::Index>> result{};

    auto dommask = low_window<max_components>(comp_idx2);

    for(std::size_t k = 0u; k < enh_bound; ++k) {
        auto enh_cod = bitsWithPop<max_components>(numx, k);
        auto enh_dom = (enh_cod & dommask) | ((enh_cod & ~dommask) << 1);

        // Codomain has 'x' at idx1.
        if(enh_cod.test(comp_idx1)) {
            result.emplace_back(k, bitsIndex(enh_dom));
            enh_dom.set(comp_idx1, false);
            enh_dom.set(comp_idx2, true);
            result.emplace_back(k, bitsIndex(enh_dom));
        }
        // Cocomain has '1' at idx1.
        else {
            result.emplace_back(k, bitsIndex(enh_dom));
        }
    }

    return result;
}

//! Generate matrices for comultiplications
//! By the definition of Khovanov complexes, the resulting matrices consist of 0 and 1.
//! \param numx The number of 'x's in enhancements on the codomain state. It turns out that the domain has to have exactly (numx-1) 'x's.
//! \param ncomps_in_cod The number of connected components in the codomain state. Since this matrix represents the comultiplication, the domain has to have exactly (ncomps_in_cod - 1) components.
//! \param comp_idx1 The common index of the components in the domain and the codomain states which are involved with the comultiplication.
//! \param comp_idx2 The index of the other component in the codomain involved with the multiplication.
//! \remark The function assumes idx1 < idx2.
//! \return The table of positions in the matrix where coefficients are 1.
static
std::vector<std::pair<int,int>>
matrix_comult(
    int numx, int ncomps_in_cod,
    std::size_t comp_idx1, std::size_t comp_idx2)
    noexcept
{
    auto enh_bound = binom(ncomps_in_cod-1, numx-1);
    std::vector<std::pair<int,int>> result{};

    auto codmask = low_window<max_components>(comp_idx2);

    for(std::size_t k = 0u; k < enh_bound; ++k) {
        auto enh_dom = bitsWithPop<max_components>(numx-1, k);
        auto enh_cod = (enh_dom & codmask) | ((enh_dom & ~codmask) << 1);

        // Domomain has 'x' at idx1.
        if(enh_dom.test(comp_idx1)) {
            enh_cod.set(comp_idx2, true);
            result.emplace_back(bitsIndex(enh_cod), k);
        }
        // Domain has '1' at idx1.
        else {
            enh_cod.set(comp_idx2, true);
            result.emplace_back(bitsIndex(enh_cod), k);
            enh_cod.set(comp_idx1, true);
            enh_cod.set(comp_idx2, false);
            result.emplace_back(bitsIndex(enh_cod), k);
        }
    }

    return result;
}

/***********************************
 *** Implementation of functions ***
 ***********************************/
// Compute Khovanov complex of a given link diagram.
std::optional<ChainIntegral>
khover::khChain(
    LinkDiagram const& diagram,
    SmoothCube const& cube,
    int qdeg)
    noexcept
{
    // Compute the enhancement data
    std::vector<EnhancementProperty> enh_prop{};
    // Check if q-degree has compatible parity
    if(auto aux = get_enhancement_prop(diagram, cube, qdeg); !aux)
        return std::nullopt;
    else
        enh_prop = std::move(*aux);

    // Compute the matrices representing matrices.
    ChainIntegral result(
        -diagram.npositive()-1,
        matrix_t(
            0, binom(cube.back().ncomp, enh_prop.back().xcnt)));
    for(int i = cube.dim(); i > 0; --i) {
        std::size_t maxst_cod
            = (low_window<max_crosses>(i)<<(diagram.ncrosses()-i)).to_ullong();
        std::size_t maxst_dom
            = (low_window<max_crosses>(i-1)<<(diagram.ncrosses()-i+1)).to_ullong();
        matrix_t diffmat = matrix_t::Zero(
            enh_prop[maxst_cod].headidx
            + binom(cube[maxst_cod].ncomp, enh_prop[maxst_cod].xcnt),
            enh_prop[maxst_dom].headidx
            + binom(cube[maxst_dom].ncomp, enh_prop[maxst_dom].xcnt));

        // Traverse all the state pairs
        for(std::size_t stidx_cod = 0;
            stidx_cod < binom(cube.dim(), i);
            ++stidx_cod)
        {
            // The state associated with the index.
            auto st_cod = bitsWithPop<max_crosses>(i, stidx_cod).to_ullong();

            for(std::size_t stidx_dom = 0;
                stidx_dom < binom(cube.dim(), i-1);
                ++stidx_dom)
            {
                // The state associated with the index.
                auto st_dom = bitsWithPop<max_crosses>(i-1, stidx_dom).to_ullong();

                // The coefficient between the domain/codomain states.
                auto coeff = diagram.stateCoeff(
                    state_t{st_dom}, state_t{st_cod});

                // Skip the case where the codomain state is not adjacent to the domain state.
                if(coeff == 0)
                    continue;

                // Find the position where saddle operation is applied.
                for(std::size_t arc = 0; arc < diagram.narcs(); ++arc) {
                    // The saddle causes multiplication.
                    if(cube[st_cod].arccomp[arc] < cube[st_dom].arccomp[arc])
                    {
                        for(auto [r,c] : matrix_mult(
                                enh_prop[st_cod].xcnt,
                                cube[st_cod].ncomp,
                                cube[st_cod].arccomp[arc],
                                cube[st_dom].arccomp[arc])) {
                            diffmat.coeffRef(
                                enh_prop[st_cod].headidx+r,
                                enh_prop[st_dom].headidx+c
                                ) += coeff;
                        }
                        break;
                    }
                    // The saddle causes comultiplication
                    else if(cube[st_cod].arccomp[arc]
                            > cube[st_dom].arccomp[arc])
                    {
                        for(auto [r,c] : matrix_comult(
                                enh_prop[st_cod].xcnt,
                                cube[st_cod].ncomp,
                                cube[st_dom].arccomp[arc],
                                cube[st_cod].arccomp[arc]))
                        {
                            diffmat.coeffRef(
                                enh_prop[st_cod].headidx+r,
                                enh_prop[st_dom].headidx+c
                                ) += coeff;
                        }
                        break;
                    }
                }
            }
        }

        // Append the matrix as a differential.
        if (!result.prepend(std::move(diffmat))) {
            std::cerr << "Wrong size matrix..." << std::endl;
            return std::nullopt;
        }
    }

    return std::make_optional(std::move(result));
}

// Compute crux complex of a given link diagram and a given crossing.
std::optional<ChainIntegral>
khover::cruxChain(
    LinkDiagram const& diagram,
    std::size_t dblpt,
    CruxCube const& cube,
    int qdeg)
    noexcept
{
    // The double point is out-of-range.
    if(dblpt >= diagram.ncrosses()) {
        return std::nullopt;
    }

    // Compute the enhancement data
    std::vector<EnhancementProperty> enh_prop{};
    // Check if q-degree has compatible parity
    if(auto aux = get_enhancement_prop(
           diagram, cube, diagram.getSign(dblpt) > 0 ? qdeg : qdeg-2); !aux)
        return std::nullopt;
    else
        enh_prop = std::move(*aux);

    // Compute the matrices representing matrices.
    ChainIntegral result(
        -diagram.npositive() - (diagram.getSign(dblpt) > 0 ? 0 : 1),
        matrix_t(0, binom(cube.back().ncomp, enh_prop.back().xcnt))
        );

    for(int i = diagram.ncrosses()-1; i > 0; --i) {
        std::size_t maxst_cod
            = (low_window<max_crosses>(i)<<(diagram.ncrosses()-i-1)).to_ullong();
        std::size_t maxst_dom
            = (low_window<max_crosses>(i-1)<<(diagram.ncrosses()-i)).to_ullong();
        matrix_t diffmat = matrix_t::Zero(
            enh_prop[maxst_cod].headidx
            + binom(cube[maxst_cod].ncomp, enh_prop[maxst_cod].xcnt),
            enh_prop[maxst_dom].headidx
            + binom(cube[maxst_dom].ncomp, enh_prop[maxst_dom].xcnt)
            );

        // Traverse all the state pairs
        for(std::size_t stidx_cod = 0;
            stidx_cod < binom(diagram.ncrosses()-1, i);
            ++stidx_cod)
        {
            // The state associated with the index.
            auto stbits_cod = bitsWithPop<max_crosses>(i, stidx_cod);
            auto st_cod = stbits_cod.to_ulong();

            // Skip states that has no enhancement in the q-degree.
            if(enh_prop[st_cod].xcnt < 0
               || (enh_prop[st_cod].xcnt
                   > static_cast<int>(cube[st_cod].ncomp)))
            {
                continue;
            }

            for(std::size_t c = 0; c < diagram.ncrosses()-1; ++c) {
                if(!stbits_cod.test(c))
                    continue;

                auto stbits_dom = stbits_cod;
                stbits_dom.set(c,false);
                auto st_dom = stbits_dom.to_ulong();

                // Skip states that has no enhancement in the q-degree.
                if(enh_prop[st_dom].xcnt < 0
                   || (enh_prop[st_dom].xcnt
                       > static_cast<int>(cube[st_cod].ncomp)))
                {
                    continue;
                }

                // The coefficient between the domain/codomain states.
                auto coeff = diagram.stateCoeff(
                    cube[st_dom].state, cube[st_cod].state);

                // Skip the case where the codomain state is not adjacent to the domain state.
                if(coeff == 0)
                    continue;

                // Find the position where saddle operation is applied.
                for(std::size_t arc = 0; arc < diagram.narcs(); ++arc) {
                    // The saddle causes multiplication.
                    if(cube[st_cod].arccomp[arc]
                       < cube[st_dom].arccomp[arc])
                    {
                        // Enabled if the operation is involved with twisted arcs.
                        std::optional<std::size_t> action_arc;
                        if ((cube[st_dom].twist
                             ^ cube[st_cod].twist).any())
                        {
                            // In that case, action_arc is a non-twisted arc that acts on the twisted arc.
                            if(cube[st_dom].twist.test(arc))
                                action_arc = cube[st_cod].arccomp[arc];
                            else
                                action_arc = cube[st_dom].arccomp[arc];
                        }

                        for(auto [r,c] : matrix_mult(
                                enh_prop[st_cod].xcnt,
                                cube[st_cod].ncomp,
                                cube[st_cod].arccomp[arc],
                                cube[st_dom].arccomp[arc]))
                        {
                            // -1 if 'x' on the act circle.
                            diffmat.coeffRef(
                                enh_prop[st_cod].headidx+r,
                                enh_prop[st_dom].headidx+c
                                ) += action_arc && bitsWithPop<max_components>(
                                    enh_prop[st_dom].xcnt, c).test(*action_arc)
                                ? -coeff : coeff;
                        }
                        break;
                    }
                    // The saddle causes comultiplication
                    else if (cube[st_cod].arccomp[arc]
                             > cube[st_dom].arccomp[arc])
                    {
                        // Enabled if the operation is involved with twisted arcs.
                        std::optional<std::size_t> coact_arc;
                        if ((cube[st_dom].twist
                             ^ cube[st_cod].twist).any())
                        {
                            // In that case, coact_arc is a non-twisted arc that coacts on the twisted arc.
                            if(cube[st_cod].twist.test(arc))
                                coact_arc = cube[st_dom].arccomp[arc];
                            else
                                coact_arc = cube[st_cod].arccomp[arc];
                        }

                        for(auto [r,c] : matrix_comult(
                                enh_prop[st_cod].xcnt,
                                cube[st_cod].ncomp,
                                cube[st_dom].arccomp[arc],
                                cube[st_cod].arccomp[arc]))
                        {
                            // -1 if '1' on the coact circle.
                            diffmat.coeffRef(
                                enh_prop[st_cod].headidx+r,
                                enh_prop[st_dom].headidx+c
                                ) += coact_arc && bitsWithPop<max_components>(
                                    enh_prop[st_cod].xcnt, r).test(*coact_arc)
                                ? coeff : -coeff;
                        }
                        break;
                    }
                }
            }
        }

        // Append the matrix as a differential.
        if (!result.prepend(std::move(diffmat))) {
            std::cerr << "Wrong size matrix..." << std::endl;
            return std::nullopt;
        }
    }

    return std::optional<ChainIntegral>(std::move(result));
}

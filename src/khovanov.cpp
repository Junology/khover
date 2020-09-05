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
#include "enhancements.hpp"

/* Debug
#include "debug/debug.hpp"
//*/

using namespace khover;

using matrix_t = ChainIntegral::matrix_t;


/***********************************
 *** Implementation of functions ***
 ***********************************/
// Compute Khovanov complex of a given link diagram.
std::optional<ChainIntegral>
khover::khChain(
    LinkDiagram const& diagram,
    SmoothCube const& cube,
    std::vector<EnhancementProperty> const& enh_prop
    //int qdeg
    ) noexcept
{
    // Compute the matrices representing matrices.
    ChainIntegral result(
        -diagram.npositive()-1,
        matrix_t(
            0, binom(cube.back().ncomp, enh_prop.back().xcnt)));
    for(int i = cube.dim(); i > 0; --i) {
        std::size_t maxst_cod = cube.maxState(i).to_ulong();
        std::size_t maxst_dom = cube.maxState(i-1).to_ulong();

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
    std::vector<EnhancementProperty> const& enh_prop
    //int qdeg
    ) noexcept
{
    // The double point is out-of-range.
    if(dblpt >= diagram.ncrosses()) {
        return std::nullopt;
    }

    // Compute the matrices representing matrices.
    ChainIntegral result(
        -diagram.npositive() - (diagram.getSign(dblpt) > 0 ? 0 : 1),
        matrix_t(0, binom(cube.back().ncomp, enh_prop.back().xcnt))
        );

    for(int i = diagram.ncrosses()-1; i > 0; --i) {
        std::size_t maxst_cod = cube.maxState(i).to_ulong();
        std::size_t maxst_dom = cube.maxState(i-1).to_ulong();
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

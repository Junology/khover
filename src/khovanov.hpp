/*!
 * \file khovanov.hpp
 * \author Jun Yoshida
 * \copyright (c) 2020 Jun Yoshida.
 * The project is released under the 2-clause BSD License.
 * \date August, 2020: created
 */

#pragma once

#include <optional>

#include "states.hpp"
#include "linkdiagram.hpp"
#include "cubes.hpp"
#include "chaincomplex.hpp"

namespace khover {

//! The type carrying the data of enhancements on each state.
struct EnhancementProperty {
    //! The index of the first enhancement on a state among a certain (co)homological degree.
    ChainIntegral::matrix_t::Index headidx;

    //! The number of 'x' in the enhancements.
    //! This might be negative or exceed the number of components in case no enhancement is allowed in the q-degree.
    int xcnt = -1;
};


//! Compute Khovanov complex of a given link diagram.
//! \param diagram The target diagram.
//! \cube This must be the value of SmoothCube::fromDiagram function for *diagram*.
//! \retval std::nullopt If the parity of qdeg is incorrect in terms of the signature of the diagram.
std::optional<ChainIntegral>
khChain(
    LinkDiagram const& diagram,
    SmoothCube const& cube,
    int qdeg) noexcept;

//! Compute crux complex of a given link diagram and a given crossing.
//! \param diagram The target diagram.
//! \param dblpt A crossing in the diagram that is regarded as a double point.
//! \param cube This must be the value of CruxCube::fromDiagram function for *diagram* and *dblpt*.
//! \retval std::nullopt If the parity of qdeg is incorrect in terms of the signature of the diagram.
std::optional<ChainIntegral>
cruxChain(
    LinkDiagram const& diagram,
    std::size_t dblpt,
    CruxCube const& cube,
    int qdeg) noexcept;

//! The morphism Xi whose mapping cone is the first Vassiliev derivative of the Khovanov homology.
std::optional<ChainIntegral::Hom>
cruxXi(LinkDiagram diagram, std::size_t dblpt, int qdeg) noexcept;

} // end namespace khover

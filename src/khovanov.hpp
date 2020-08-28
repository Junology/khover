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

} // end namespace khover

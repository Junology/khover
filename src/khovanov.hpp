/*!
 * \file khovanov.hpp
 * \author Jun Yoshida
 * \copyright (c) 2020 Jun Yoshida.
 * The project is released under the 2-clause BSD License.
 * \date August, 2020: created
 */

#pragma once

#include <optional>

#include "chaincomplex.hpp"
#include "linkdiagram.hpp"
#include "states.hpp"
#include "enhancements.hpp"
#include "cubes.hpp"

namespace khover {

//! Compute Khovanov complex of a given link diagram.
//! \param diagram The target diagram.
//! \param cube This must be the value of SmoothCube::fromDiagram function for *diagram*.
//! \param enh_prop The list of enhancement properties associated with *cube*.
//! \retval std::nullopt If the parity of qdeg is incorrect in terms of the signature of the diagram.
std::optional<ChainIntegral>
khChain(
    LinkDiagram const& diagram,
    SmoothCube const& cube,
    std::vector<EnhancementProperty> const& enh_prop
    // int qdeg
    ) noexcept;

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
    std::vector<EnhancementProperty> const& enh_prop
    //int qdeg
    ) noexcept;

//! Compute the morphism PhiHat.
std::optional<ChainIntegral::Hom>
crossPhiHat(
    LinkDiagram const& diagram,
    std::size_t crossing,
    SmoothCube const& cube_neg,
    std::vector<EnhancementProperty> const& enhprop_neg,
    SmoothCube const& cube_pos,
    std::vector<EnhancementProperty> const& enhprop_pos
    ) noexcept;

//! The morphism Xi whose mapping cone is the first Vassiliev derivative of the Khovanov homology.
std::optional<ChainIntegral::Hom>
cruxXi(
    LinkDiagram const& diagram,
    CruxCube const& cubeCrx,
    LinkDiagram const& diagV,
    SmoothCube const& cubeV,
    LinkDiagram const& diagW,
    SmoothCube const& cubeW,
    std::size_t dblpt,
    int qdeg
    ) noexcept;

} // end namespace khover

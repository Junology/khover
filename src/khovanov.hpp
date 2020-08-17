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
#include "chaincomplex.hpp"

namespace khover {

//! Compute Khovanov complex of a given link diagram.
//! \retval std::nullopt If the parity of qdeg is incorrect in terms of the signature of the diagram.
std::optional<ChainIntegral>
khChain(LinkDiagram const& diagram, int qdeg) noexcept;

} // end namespace khover

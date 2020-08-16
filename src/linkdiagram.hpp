/*!
 * \file linkdiagram.hpp
 * \author Jun Yoshida
 * \copyright (c) 2020 Jun Yoshida.
 * The project is released under the 2-clause BSD License.
 * \date August, 2020: created
 */

#pragma once

#include <utility>
#include <bitset>
#include <vector>
#include <optional>

#include "states.hpp"

//* Debug
#include <iostream>
// */

namespace khover {

//! Internal representation of link diagrams.
//! It only contains minimum information for "state-sum" methods.
class LinkDiagram {
public:
    //! The type representing each crossings.
    struct Crossing {
        //! A flag if the crossing is positive or not.
        bool is_positive;
        //! Indicating how the crossing is adjacent to arcs so that
        //!   - *adj_arc[0]* is the in-coming over-arc;
        //!   - *adj_arc[1]* is the out-going over-arc;
        //!   - *adj_arc[2]* is the in-coming under-arc;
        //!   - *adj_arc[3]* is the out-going under-arc.
        //! Here the variable assumes arcs are indexed by unsigned integers.
        unsigned int adj_arc[4];
    };

private:
    std::size_t m_numarcs;
    std::vector<Crossing> m_cross;

public:
    /*!
     * \name Constructors, Destructors, and assignment operators.
     */
    //\{
    //! This is the only chance to set the number of arcs (except copies).
    LinkDiagram(std::size_t numarcs, std::vector<Crossing> const& cross) noexcept
        : m_numarcs(numarcs), m_cross(cross)
    {}

    LinkDiagram(std::size_t numarcs, std::vector<Crossing>&& cross) noexcept
        : m_numarcs(numarcs), m_cross(std::move(cross))
    {}

    ~LinkDiagram() = default;

    //! The class is copy constructible.
    LinkDiagram(LinkDiagram const&) = default;

    //! The class is move constructible.
    LinkDiagram(LinkDiagram&&) = default;

    //! The class is copyable.
    LinkDiagram& operator=(LinkDiagram const&) = default;

    //! The class is movable.
    LinkDiagram& operator=(LinkDiagram &&) = default;
    //\}

    /*!
     * \name Getters.
     */
    //\{
    //! Get the number of arcs.
    inline
    std::size_t
    narcs() const noexcept { return m_numarcs; }

    //! Get the number of crossings.
    inline
    std::size_t
    ncrosses() const noexcept { return m_cross.size(); }

    //! Compute writhe
    int writhe() const noexcept;

    //! Get the list of crossings.
    inline
    const std::vector<Crossing>&
    crosses() const noexcept { return m_cross; }

    //\}

    /*!
     * \name State manipulation.
     */
    //\{
    //! Compute the cohomological degrees of a given state in the convention of the following article:
    //! Clark, David; Morrison, Scott; Walker, Kevin. Fixing the functoriality of Khovanov homology. Geom. Topol. 13 (2009), no.3, 1499--1582. doi:10.2140/gt.2009.13.1499.
    int cohDegree(state_t st) const noexcept;

    //! Check if a state is adjacent to the other in Khovanov's smoothing cube.
    //! \param st_before A state that should be of a lower cohomological degree.
    //! \param st_after A state that should be of a higher cohomological degree.
    //! \retval 0 Two states are not adjacent to each other.
    //! \retval 1 Two states are adjacent with positive sign.
    //! \retval -1 Two states are adjacent with negative sign.
    int stateCoeff(state_t st_before, state_t st_after) const noexcept;

    //! Compute the connected components in the smoothing of the diagram corresponding to a given state.
    //! Components are indexed by consecutive non-negative integers begining from 0 in the order so that componens with smaller indices contain arcs with smaller indices.
    //! \retval (n,c) For each 0 <= i narcs()-1, c[i] is the index of the component that contains the i-th arc, and n is the number of components.
    std::pair<std::size_t,std::vector<component_t>>
    smoothing(state_t st) const noexcept;
    //\}
};

//! Type used to represent Gauss codes.
using gcode_t = std::vector<int>;

/*!
 * Read a Gauss code.
 * One can pass a multi-component code by separating codes with '0'.
 * \param gcode A Gauss code to process.
 * \param signs A list of signs of the crossings with the minimum indices in *disjoint* components, where disjoint means they are disjoint in the 2d-plane (not in the 3d-space). Hence, only one entry suffices in usual cases.
 */
std::optional<LinkDiagram> read_gauss_code(
    gcode_t const& gcode,
    std::vector<bool> const& signs
    ) noexcept;

} // end namespace khover

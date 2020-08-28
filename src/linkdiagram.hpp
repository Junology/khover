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
#include <algorithm>

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

    using smoothdata_t = Smoothing<std::vector<component_t>>;
    using smoothtwist_t = SmoothTw<std::vector<component_t>>;

private:
    std::size_t m_numarcs;
    std::vector<Crossing> m_cross;

    // Be careful on the order of declarations.
    std::size_t m_numpositive;
    std::size_t m_numnegative;

public:
    /*!
     * \name Constructors, Destructors, and assignment operators.
     */
    //\{
    //! This is the only chance to set the number of arcs (except copies).
    LinkDiagram(std::size_t numarcs, std::vector<Crossing> const& cross) noexcept
        : m_numarcs(numarcs), m_cross(cross),
          m_numpositive(
              std::count_if(
                  std::begin(cross), std::end(cross),
                  [](auto const& crs) { return crs.is_positive; } )),
          m_numnegative(cross.size() - m_numpositive)
    {}

    LinkDiagram(std::size_t numarcs, std::vector<Crossing>&& cross) noexcept
        : m_numarcs(numarcs), m_cross(std::move(cross)),
          m_numpositive(
              std::count_if(
                  std::begin(m_cross), std::end(m_cross),
                  [](auto const& crs) { return crs.is_positive; } )),
          m_numnegative(m_cross.size() - m_numpositive)
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

    //! Get the number of negative crossings.
    inline
    std::size_t
    nnegative() const noexcept { return m_numnegative; }

    //! Get the number of positive crossings.
    inline
    std::size_t
    npositive() const noexcept { return m_numpositive; }

    //! Compute the writhe number
    inline
    int
    writhe() const noexcept { return m_numpositive - m_numnegative; }

    //! Get the list of crossings.
    inline
    const std::vector<Crossing>&
    crosses() const noexcept { return m_cross; }

    //! Get the sign of a crossing.
    //! \retval -1 The crossing is negative.
    //! \retval 1 The crossing is positive.
    //! \retval 0 Out-of-range.
    inline
    int
    getSign(std::size_t c) const noexcept {
        return c < m_cross.size()
                   ? (m_cross[c].is_positive ? 1 : -1) : 0;
    }

    //\}


    /*!
     * \name Crossing manipulation.
     */
    //\{

    //! Swap two crossings.
    inline
    void
    swapCrossings(std::size_t c0, std::size_t c1) noexcept {
        if (c0 < m_cross.size() && c1 < m_cross.size() && c0 != c1)
            std::swap(m_cross[c0], m_cross[c1]);
    }

    //! Crossing change
    inline
    void
    crossingChange(std::size_t c) noexcept {
        if (c < m_cross.size())
            m_cross[c].is_positive ^= true;
    }

    //! Make a crossing positive.
    inline
    void
    makePositive(std::size_t c) noexcept {
        if (c < m_cross.size())
            m_cross[c].is_positive = true;
    }

    //! Make a crossing negative.
    inline
    void
    makeNegative(std::size_t c) noexcept {
        if (c < m_cross.size())
            m_cross[c].is_positive = false;
    }

    //\}


    /*!
     * \name State manipulation.
     */
    //\{
    //! Compute the cohomological degrees of a given state in the convention of the following article:
    //! Clark, David; Morrison, Scott; Walker, Kevin. Fixing the functoriality of Khovanov homology. Geom. Topol. 13 (2009), no.3, 1499--1582. doi:10.2140/gt.2009.13.1499.
    inline
    int
    cohDegree(state_t st) const noexcept {
        return static_cast<int>(st.count()) - static_cast<int>(m_numnegative);
    }

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

    //! Determine which arcs are twisted in the crux complex.
    //! \param st A state on crossings except the double point. This implies that indices in the variable *st* may differ from those in the diagram in the case they have larger than that of the double point. To avoid this, one can ensure the double point has the largest index by using *swapCrossings* member function.
    //! \param dblpt The crossing which is regarded as a double point in the crux complex.
    //! \return If st is a crux state, then the returned value is a pair of
    //!   - a state on the oridinal diagram so that the double point is smoothed along the orientation; and
    //!   - a set of flags indicating whether each arcs are twisted or not.
    //! Otherwise, std::nullopt;
    std::optional<std::pair<state_t,std::bitset<max_arcs>>>
    cruxTwists(state_t st, std::size_t dblpt) const noexcept;
    //\}
};

//! Type used to represent Gauss codes.
using gcode_t = std::vector<int>;

/*!
 * Read a Gauss code.
 * One can pass a multi-component code by separating codes with '0'.
 * \param gcode A Gauss code to process.
 * \param signs A list of pairs of crossings and requested signs on them. If one contradicts to another, the latter request will be applied. If a given crossing is out-of-range, then the pair is just ignored.
 */
std::optional<LinkDiagram> read_gauss_code(
    gcode_t const& gcode,
    std::vector<std::pair<std::size_t,bool>> const& signs
    ) noexcept;

} // end namespace khover

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
    //! The type of indices of crossings;
    using crossing_t = std::size_t;

    //! The type of tables of signs.
    //! It is supposed to be convertible to the type of states.
    using signs_t = state_t;
    using arc_t = unsigned int;

    //! The type representing each crossings.
    struct Crossing {
        //! Indicating how the crossing is adjacent to arcs so that
        //!   - *adj_arc[0]* is the in-coming over-arc;
        //!   - *adj_arc[1]* is the out-going over-arc;
        //!   - *adj_arc[2]* is the in-coming under-arc;
        //!   - *adj_arc[3]* is the out-going under-arc.
        //! Here the variable assumes arcs are indexed by unsigned integers.
        arc_t adj_arc[4];
    };

private:
    //! The number of all arcs.
    std::size_t m_numarcs;

    //! The list of crossings.
    //! The length must be < khover::Limits::max_crosses.
    std::vector<Crossing> m_cross;

    //! The table of signs of crossings.
    //! 0: negative, 1: positive.
    signs_t m_signs;

    //! The list of wide edges; i.e. infinity-smoothed crossings.
    //! (see Khovanov and Rozansky, "Matrix factorizations and link homology I")
    std::vector<Crossing> m_wides{};

public:
    /*!
     * \name Constructors, Destructors, and assignment operators.
     */
    //\{
    //! This is the only chance to set the number of arcs (except copies).
    template <template <class...> class Container>
    LinkDiagram(std::size_t numarcs, Container<Crossing> const& cross, signs_t const& signs) noexcept
        : m_numarcs(numarcs), m_cross(cross), m_signs(signs)
    {}

    LinkDiagram(std::size_t numarcs, std::vector<Crossing>&& cross, signs_t const& signs) noexcept
        : m_numarcs(numarcs), m_cross(std::move(cross)), m_signs(signs)
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

    //! Get the number of positive crossings.
    inline
    std::size_t
    npositive() const noexcept { return m_signs.count(); }

    //! Get the number of negative crossings.
    inline
    std::size_t
    nnegative() const noexcept { return m_cross.size() - npositive(); }

    //! Compute the writhe number
    inline
    int
    writhe() const noexcept { return npositive() - nnegative(); }

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
    getSign(crossing_t c) const noexcept {
        return c < m_cross.size()
                   ? (m_signs[c] ? 1 : -1) : 0;
    }

    //! Get the index of the upper-right arc.
    inline
    arc_t getURArc(crossing_t c) const noexcept {
        return c < m_cross.size()
                   ? (m_signs[c]
                      ? m_cross[c].adj_arc[1]
                      : m_cross[c].adj_arc[3])
                   : std::numeric_limits<arc_t>::max();
    }

    //! Get the index of the upper-left arc.
    inline
    arc_t getULArc(crossing_t c) const noexcept {
        return c < m_cross.size()
                   ? (m_signs[c]
                      ? m_cross[c].adj_arc[3]
                      : m_cross[c].adj_arc[1])
                   : std::numeric_limits<arc_t>::max();
    }

    //! Get the index of the lower-right arc.
    inline
    arc_t getDRArc(crossing_t c) const noexcept {
        return c < m_cross.size()
                   ? (m_signs[c]
                      ? m_cross[c].adj_arc[2]
                      : m_cross[c].adj_arc[0])
                   : std::numeric_limits<arc_t>::max();
    }

    //! Get the index of the lower-left arc.
    inline
    arc_t getDLArc(crossing_t c) const noexcept {
        return c < m_cross.size()
                   ? (m_signs[c]
                      ? m_cross[c].adj_arc[0]
                      : m_cross[c].adj_arc[2])
                   : std::numeric_limits<arc_t>::max();
    }
    //\}


    /*!
     * \name Crossing manipulation.
     */
    //\{

    //! Swap two crossings.
    inline
    void
    swapCrossings(crossing_t c0, crossing_t c1) noexcept {
        if (c0 < m_cross.size() && c1 < m_cross.size() && c0 != c1)
            std::swap(m_cross[c0], m_cross[c1]);
    }

    //! Crossing change
    inline
    void
    crossingChange(crossing_t c) noexcept {
        if (c < m_cross.size())
            m_signs.flip(c);
    }

    //! Take the mirror image.
    //! i.e. Flipping the signs of all crossings.
    inline
    void
    mirroring() noexcept {
        m_signs.flip();
    }

    //! Make a crossing positive.
    inline
    void
    makePositive(crossing_t c) noexcept {
        if (c < m_cross.size())
            m_signs.set(c);
    }

    //! Make a crossing negative.
    inline
    void
    makeNegative(crossing_t c) noexcept {
        if (c < m_cross.size())
            m_signs.reset(c);
    }

    //! Vertical smoothing; i.e. the smoothing along the orientation.
    //! As some arcs are merged, the indexing on arcs, as well as those on crossings, may be changed.
    //! \remark This function may create an isolated circle.
    //! \warning This operation actually removes a crossing and hence is non-invertible.
    void makeSmoothV(crossing_t c) noexcept;

    //! Make a crossing into a "wide edge" (see [Khovanov and Rozansky]).
    //! \remark Although essentially no information is lost with this function, there is currently no way to recover the crossing.
    inline void makeWide(crossing_t c) noexcept
    {
        if (c < m_cross.size()) {
            m_wides.push_back(m_cross[c]);
            m_cross.erase(std::next(std::begin(m_cross), c));
            omit_bit(m_signs, c);
        }
    }
    //\}


    /*!
     * \name State manipulation.
     */
    //!\{
    //! Compute the cohomological degrees of a given state in the convention of the following article:
    //! Clark, David; Morrison, Scott; Walker, Kevin. Fixing the functoriality of Khovanov homology. Geom. Topol. 13 (2009), no.3, 1499--1582. doi:10.2140/gt.2009.13.1499.
    inline
    int
    cohDegree(state_t st) const noexcept {
        return static_cast<int>(st.count()) - static_cast<int>(nnegative());
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
    //!\}

protected:
    //! \name Miscellaneous functions used for implementations.
    //!\{

    //! Merge arcs
    inline bool mergeArcs(arc_t arc1, arc_t arc2) noexcept {
        if (arc1 == arc2 || arc1 >= m_numarcs || arc2 >= m_numarcs)
            return false;

        auto [arcmin,arcmax] = arc1 < arc2
            ? std::make_pair(arc1,arc2)
            : std::make_pair(arc2,arc1);
        for(auto& crs : m_cross) {
            for(std::size_t i = 0; i < 4; ++i) {
                if (crs.adj_arc[i] == arcmax)
                    crs.adj_arc[i] = arcmin;
                else if(crs.adj_arc[i] > arcmax)
                    --(crs.adj_arc[i]);
            }
        }
        --m_numarcs;
        return true;
    }
    //!\}
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

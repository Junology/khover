/*!
 * \file linkdiagram.cpp
 * \author Jun Yoshida
 * \copyright (c) 2020 Jun Yoshida.
 * The project is released under the 2-clause BSD License.
 * \date August, 2020: created
 */

#include <map>
#include <algorithm>
#include "linkdiagram.hpp"

#include "utils.hpp"

#include "debug/debug.hpp"

using namespace khover;


/*************************
 *** Utility functions ***
 *************************/

//! Merge two components containing designated two arcs.
//! \return The index of the component, if any, that was forgotten in merging.
static
std::optional<component_t>
merge_components(
    std::vector<component_t> &comptbl,
    std::size_t i1, std::size_t i2)
    noexcept
{
    // If two indices are the same, there is nothing to do.
    component_t c1 = comptbl[i1], c2 = comptbl[i2];

    if (c2 == c1)
        return std::nullopt;

    // We may assume c1 < c2
    if (c1 > c2)
        std::swap(c1,c2);

    for(auto& idx : comptbl) {
        if (idx == c2)
            idx = c1;
    }

    return c2;
}


/***************************************
 *** Implementation of class members ***
 ***************************************/

// Check if a state is adjacent to the other in Khovanov's smoothing cube.
int LinkDiagram::stateCoeff(state_t st_before, state_t st_after)
    const noexcept
{
    // Ignore higher bits.
    st_before &= low_window<max_crosses>(m_cross.size());
    st_after &= low_window<max_crosses>(m_cross.size());

    // Check if they are adjacent with respect to the cohomological degrees.
    if (st_after.count() - st_before.count() != 1)
        return 0;

    // Find different bits.
    auto st_diff = st_after ^ st_before;

    // If there are more than one differences, two states are not adjacent.
    if (st_diff.count() != 1)
        return 0;

    // The parity of the number of smoothings that "break the orientaion."
    bool sign = false;
    for(std::size_t i = 0; i < m_cross.size() && !st_diff.test(i); ++i) {
        sign ^=
            (m_cross[i].is_positive && st_before.test(i))
            || (!m_cross[i].is_positive && !st_before.test(i));
    }

    return sign ? -1 : 1;
}

// Vertical smoothing; i.e. smoothing along the orientation.
void
khover::LinkDiagram::makeSmoothV(crossing_t c) noexcept
{
    // Bound check.
    if (c >= m_cross.size())
        return;

    // Merge arcs
    mergeArcs(m_cross[c].adj_arc[0], m_cross[c].adj_arc[3]);
    mergeArcs(m_cross[c].adj_arc[1], m_cross[c].adj_arc[2]);

    // Remove the crossing
    if(m_cross[c].is_positive)
        --m_numpositive;
    else
        --m_numnegative;

    m_cross.erase(std::next(std::begin(m_cross),c));
}


// Compute the connected components in the smoothing of the diagram corresponding to a given state.
std::pair<std::size_t, std::vector<component_t>>
khover::LinkDiagram::smoothing(state_t st)
    const noexcept
{
    std::vector<component_t> result(m_numarcs);
    std::bitset<max_components> index_rmed{0u};

    // Begin with discrete spaces
    for(std::size_t i = 0; i < m_numarcs; ++i) {
        result[i] = i;
    }

    // Smoothing on wide edges.
    for(auto& w : m_wides) {
        if (auto rmc = merge_components(result, w.adj_arc[0], w.adj_arc[2]); rmc)
            index_rmed.set(*rmc);
        if (auto rmc = merge_components(result, w.adj_arc[1], w.adj_arc[3]); rmc)
            index_rmed.set(*rmc);
    }

    // Smoothing of crossings.
    for(std::size_t i = 0; i < m_cross.size(); ++i) {
        // state=1 on positive or state=0 on negative
        if(m_cross[i].is_positive == st.test(i)) {
            if (auto rmc = merge_components(
                    result, m_cross[i].adj_arc[0], m_cross[i].adj_arc[2]);
                rmc)
            {
                index_rmed.set(*rmc);
            }
            if (auto rmc = merge_components(
                    result, m_cross[i].adj_arc[1], m_cross[i].adj_arc[3]);
                rmc)
            {
                index_rmed.set(*rmc);
            }
        }
        // state=0 on positive or state=1 on negative
        else {
            if (auto rmc = merge_components(
                    result, m_cross[i].adj_arc[0], m_cross[i].adj_arc[3]);
                rmc)
            {
                index_rmed.set(*rmc);
            }
            if (auto rmc = merge_components(
                    result, m_cross[i].adj_arc[1], m_cross[i].adj_arc[2]);
                rmc)
            {
                index_rmed.set(*rmc);
            }
        }
    }

    // re-labeling the components so that indices are consecutive.
    for(auto& cind : result) {
        cind -= (index_rmed & low_window<max_components>(cind)).count();
    }

    return std::make_pair(
        m_numarcs - index_rmed.count(),
        std::move(result));
}


// Determine which arcs are twisted in the crux complex.
std::optional<std::pair<state_t,std::bitset<max_arcs>>>
LinkDiagram::cruxTwists(state_t st, std::size_t dblpt)
    const noexcept
{
    // Return immediately in case the double point is out-of-range.
    if (dblpt >= m_cross.size())
        return std::nullopt;

    // Extend the state on crossings to double point so that the double point is smoothed along the orientation.
    st = insertBit(st, dblpt, !m_cross[dblpt].is_positive);

    std::bitset<max_arcs> is_twisted{0u};
    std::size_t icur = dblpt;
    std::size_t dir = 1;

    do {
        is_twisted.set(m_cross[icur].adj_arc[dir]);

        // Find the other end of the arc.
        std::size_t i = 0, d = 0;
        for(; i < m_cross.size(); ++i) {
            for(d = 0; d < 4; ++d) {
                if (i == icur && d == dir)
                    continue;
                if(m_cross[i].adj_arc[d] == m_cross[icur].adj_arc[dir])
                    break;
            }
            if (d < 4)
                break;
        }

        // Not found => error.
        if(i >= m_cross.size()) {
            ERR_MSG("Twisted arc is not closed");
            return std::nullopt;
        }

        // Move to next.
        icur = i;
        switch(d) {
        case 0:
            dir = m_cross[i].is_positive == st.test(i) ? 2 : 3;
            break;

        case 1:
            dir = m_cross[i].is_positive == st.test(i) ? 3 : 2;
            break;

        case 2:
            dir = m_cross[i].is_positive == st.test(i) ? 0 : 1;
            break;

        case 3:
            dir = m_cross[i].is_positive == st.test(i) ? 1 : 0;
            break;

        default:
            ERR_MSG("Something bad happened.");
            return std::nullopt;
        }
    } while(icur != dblpt);

    if(dir == 1) {
        return std::nullopt;
    }
    else {
        return std::make_pair(st,is_twisted);
    }
}

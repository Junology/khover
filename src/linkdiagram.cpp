/*!
 * \file linkdiagram.cpp
 * \author Jun Yoshida
 * \copyright (c) 2020 Jun Yoshida.
 * The project is released under the 2-clause BSD License.
 * \date August, 2020: created
 */

#include <algorithm>
#include "linkdiagram.hpp"

using namespace khover;

/*************************
 *** Utility functions ***
 *************************/

//! Merge two components containing designated two arcs.
static
void merge_components(
    std::vector<LinkDiagram::component_t> &comptbl,
    std::size_t i1, std::size_t i2)
    noexcept
{
    // If two indices are the same, there is nothing to do.
    LinkDiagram::component_t c1 = comptbl[i1], c2 = comptbl[i2];

    if (c2 == c1)
        return;

    // We may assume c1 < c2
    if (c1 > c2)
        std::swap(c1,c2);

    for(auto& idx : comptbl) {
        if (idx == c2)
            idx = c1;
    }
}


/***************************************
 *** Implementation of class members ***
 ***************************************/

// Compute the writhe of the link diagram.
int khover::LinkDiagram::writhe() const noexcept
{
    int writhe = 0;
    for (auto& c : m_cross) {
        if(c.is_positive)
            ++writhe;
        else
            --writhe;
    }

    return writhe;
}

// Compute the cohomological degrees of a given state in the convention of [CMW]
int khover::LinkDiagram::cohDegree(state_t st) const noexcept
{
    int deg = 0;
    for(std::size_t i = 0; i < m_cross.size(); ++i) {
        if (m_cross[i].is_positive) {
            if (st.test(i)) ++deg;
        }
        else {
            if(!st.test(i)) --deg;
        }
    }
    return deg;
}

// Compute the connected components in the smoothing of the diagram corresponding to a given state.
std::vector<khover::LinkDiagram::component_t>
khover::LinkDiagram::smoothing(state_t st)
    const noexcept
{
    std::vector<component_t> result(m_numarcs);

    // Begin with discrete spaces
    for(std::size_t i = 0; i < m_numarcs; ++i) {
        result[i] = i;
    }

    for(std::size_t i = 0; i < m_cross.size(); ++i) {
        // state=1 on positive or state=0 on negative
        if(m_cross[i].is_positive == st.test(i)) {
            merge_components(
                result, m_cross[i].adj_arc[0], m_cross[i].adj_arc[2]);
            merge_components(
                result, m_cross[i].adj_arc[1], m_cross[i].adj_arc[3]);
        }
        // state=0 on positive or state=1 on negative
        else {
            merge_components(
                result, m_cross[i].adj_arc[0], m_cross[i].adj_arc[3]);
            merge_components(
                result, m_cross[i].adj_arc[1], m_cross[i].adj_arc[2]);
        }
    }

    return result;
}

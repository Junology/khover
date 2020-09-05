/*!
 * \file cubes.cpp
 * \author Jun Yoshida
 * \copyright (c) 2020 Jun Yoshida.
 * The project is released under the 2-clause BSD License.
 * \date August, 2020: created
 */

#include <map>

#include "cubes.hpp"

using namespace khover;

SmoothCube
SmoothCube::fromDiagram(
    LinkDiagram const& diagram,
    std::initializer_list<std::pair<std::size_t,bool>> fixed_states
    ) noexcept
{
    // Sort the fixed states; ignore all out-of-range indices.
    std::map<std::size_t,bool> fixmap{};
    std::copy_if(
        fixed_states.begin(), fixed_states.end(),
        std::inserter(fixmap, fixmap.end()),
        [ncrs = diagram.ncrosses()](std::pair<std::size_t,bool> const& p) {
            return p.first < ncrs;
        } );

    // The number of all states on crossings with no state fixed.
    SmoothCube result(
        -static_cast<int>(diagram.nnegative()),
        diagram.ncrosses() - fixmap.size());
    std::size_t nstates = cipow(2,result.dim());

    // Compute the smoothing on each state.
    for(std::size_t st = 0; st < nstates; ++st) {
        // Extend the state to all crossings including state-fixed ones.
        state_t stbit{st};
        for(auto [pos,flag] : fixmap)
            stbit = insertBit(stbit, pos, flag);

        // Get smoothing
        auto [ncomp, arccomp] = diagram.smoothing(stbit);

        // Append the data.
        result.m_smoothdata.push_back({stbit, ncomp, std::move(arccomp)});
    }

    return result;
}

CruxCube
CruxCube::fromDiagram(const LinkDiagram &diagram, std::size_t dblpt) noexcept
{
    CruxCube result = (dblpt < diagram.ncrosses() && diagram.getSign(dblpt))
        ? CruxCube(
            -static_cast<int>(diagram.nnegative())+1,
            diagram.ncrosses()-1)
        : CruxCube(
            -static_cast<int>(diagram.nnegative()),
            diagram.ncrosses());

    // The number of all states on crossings with no state fixed.
    std::size_t nstates = cipow(2,result.dim());

    // Compute the smoothing on each state.
    for(std::size_t st = 0; st < nstates; ++st) {
        if (auto iscrux = diagram.cruxTwists(st,dblpt); !iscrux) {
            result.m_smoothdata.push_back(
                {
                    {
                        insertBit(state_t{st}, dblpt, diagram.getSign(dblpt) < 0),
                        0, std::vector<component_t>{}
                    },
                    0u
                });
            continue;
        }
        // For crux states, remember the twisted arcs and compute the smoothings.
        else {
            result.m_smoothdata.push_back(
                {
                    {iscrux->first, 0, std::vector<component_t>{}},
                    iscrux->second
                });
        }

        std::tie(result.m_smoothdata.back().ncomp,
                 result.m_smoothdata.back().arccomp)
            =  diagram.smoothing(result.m_smoothdata.back().state);
    }

    return result;
}

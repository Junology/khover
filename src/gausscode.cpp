/*!
 * \file gausscode.cpp
 * \author Jun Yoshida
 * \copyright (c) 2020 Jun Yoshida.
 * The project is released under the 2-clause BSD License.
 * \date August, 2020: created
 */

#include <map>
#include <queue>
#include <stack>
#include <utility>
#include <algorithm>
#include <iterator>

#include "linkdiagram.hpp"

//* Debug
#include "debug/debug.hpp"
// */

using namespace khover;

/*************************
 *** Utility functions ***
 *************************/
static inline
std::size_t
count_crossing(gcode_t const& gcode)
    noexcept
{
    return static_cast<std::size_t>(*std::max_element(gcode.begin(), gcode.end()));
}

//! Divide Gauss code by the separator 0.
static
std::vector<gcode_t>
divide_gcode(gcode_t const& gcode)
    noexcept
{
    std::vector<gcode_t> result{};

    auto itr1 = std::begin(gcode);
    auto itr2 = itr1;

    do {
        itr2 = std::find(itr1, std::end(gcode), 0);
        result.emplace_back(itr1, itr2);

        itr1 = std::next(itr2);
        // This operation may make itr1 out-of-range if itr2 == std::end(gcode)
        // Probably this itself doesn't cause any errors/problems as long as we do not de-reference the iterator.
    } while(itr2 != std::end(gcode));

    return result;
}

/*!
 * Compute the sign enhanced version of Kauffman's dual code g* of g.
 * In the resulting code g*, the sign on a crossing means whether the out-going arc is over/under with respect to the crossing in the original code g.
 * Notice that, for every vertex in g*, the two out-going arcs are simultaneously over-arcs or under-arcs in g, so there is no ambiguity.
 * \param numcrs The number of crossings in given the Gauss code.
 * \param gcode The original Gauss code.
 * \return The dual code g* if success or the empty vector otherwise.
 */
static
std::vector<int>
compute_dual(std::size_t numcrs, gcode_t gcode)
    noexcept
{
    // Iterate over all crossings
    for(std::size_t i = 0; i < numcrs; ++i) {
        // Find the first appearance of the i-th crossing.
        auto fst = std::begin(gcode);
        auto fstbeg = fst;
        for(; fst != std::end(gcode); ++fst) {
            // Enter a new component.
            if (*fst == 0)
                fstbeg = std::next(fst);
            // Found!
            if (static_cast<std::size_t>(std::abs(*fst)) == i+1)
                break;
        }

        // If there is a leap of indices, stop the computation.
        if(fst == std::end(gcode))
            return {};

        // Find the second appearance of the i-th crossing.
        auto snd = std::next(fst);
        auto sndbeg = fstbeg;
        for(; snd != std::end(gcode); ++snd) {
            // Enter a new component.
            if (*snd == 0)
                sndbeg = std::next(snd);
            // Found!
            if (static_cast<std::size_t>(std::abs(*snd)) == i+1)
                break;
        }

        // Failed to find the second appearance of i-th crossing.
        if (snd == std::end(gcode))
            return {};

        // When the components are the same.
        if (fstbeg == sndbeg) {
            // Change the sign of the first.
            (*fst) = -(*fst);

            // Reverse the in-between sequence;
            // and change the signs of already-processed vertices in it.
            std::reverse(std::next(fst), snd);
            std::for_each(
                std::next(fst), snd,
                [n = std::abs(*fst)](auto& v) {
                    if (std::abs(v) < n) v = -v;
                });
        }
        // When the components are different.
        else {
            // Change the sign of the first.
            //(*fst) = -(*fst);

            // Find the end of the first component.
            // In this case, we can assume fstend != end.
            auto fstend = std::find(std::next(fst), std::end(gcode), 0);

            // Bring the second component to the next to the first.
            std::rotate(std::next(fstend), sndbeg, std::end(gcode));

            // Compute the second i-th crossing and its component.
            snd = std::next(fstend, std::distance(sndbeg, snd) + 1);
            sndbeg = std::next(fstend);
            auto sndend = std::find(std::next(snd), std::end(gcode), 0);

            // Make the i-th crossings on the first of the components.
            std::rotate(fstbeg, fst, fstend);
            std::rotate(sndbeg, snd, sndend);

            // Reverse the tail of the first component;
            // and change the signs of already-processed vertices in it.
            std::reverse(std::next(fstbeg), fstend);
            std::for_each(
                std::next(fstbeg), fstend,
                [n = std::abs(*fstbeg)](auto& v) {
                    if (std::abs(v) < n) v = -v;

                });

            // // Reverse the tail of the second component;
            // // and change the signs of already-processed vertices in it.
            // std::reverse(std::next(sndbeg), sndend);
            // std::for_each(
            //     std::next(sndbeg), sndend,
            //     [n = std::abs(*fst)](auto& v) {
            //         if (std::abs(v) < n) v = -v;
            //     });

            // Merge two components by removing the component separator.
            gcode.erase(fstend);
        }
    }

    return gcode;
}

//! Divide vertices (or crossings for link diagrams) accoring to g*
//! As g* yields a Jordan curve (i.e. an embedded circle in the plane) consisting of all the arcs of each disjoint component, the vertices are divided into to two sets; namely, the set of vertices inside the circle and that of vertices outside.
//! \return The list of vertices outside the circle, or std::nullopt if failed.
//! \remark The meaning of "inside"/"outside" is reversed by "Markov move I."
//! In other words, the choice is inessential in our algorithm.
static
std::optional<gcode_t>
separate_vertices(std::size_t numcrs, gcode_t const& gdual)
    noexcept
{
    // Marks put on crossings which indicate
    //  -1: the crossing is set to be "outside";
    //   1: the crossing is set to be "inside";
    //   0: undetermined.
    std::vector<int> marks(numcrs, 0);

    // Determine marks on the crossings.
    // The process repeats until all crossings are marked.
    for(auto itr = std::begin(gdual); itr != std::end(gdual);
        itr = std::find_if(
            std::begin(gdual), std::end(gdual),
            [&marks](int c) { return marks[std::abs(c)-1] == 0; }))
    {
        // Queue of setters.
        // setter_q.first : The target crossing.
        // setter_q.second: The mark to be set.
        std::queue<std::pair<int,int>> setter_q;

        // The first undetermined vertex is set to be "inside".
        setter_q.emplace(std::abs(*itr)-1,1);

        // Iterate while there is a setter in the queue.
        while(!setter_q.empty()) {
            // Alias
            auto tgt = setter_q.front().first;
            auto mrk = setter_q.front().second;

            setter_q.pop();

            // If the setter marks a crossing by the same marker, then just skip.
            if (marks[tgt] == mrk)
                continue;
            // Contradictory marking causes an error.
            else if (marks[tgt] != 0) {
                ERR_MSG("ERROR: Cannot divide vertices.");
                return std::nullopt;
            }

            // Set the mark.
            marks[tgt] = mrk;

            // Find the first appearance of the target crossing in the dual code.
            std::size_t beg = 0;
            while(beg < gdual.size()
                  && std::abs(gdual[beg]) != tgt+1)
                ++beg;

            // Count the appearance of crossings between the target crossing.
            std::vector<std::size_t> numappears(numcrs,0);
            std::size_t end = beg+1;
            while(end < gdual.size()
                  && std::abs(gdual[end]) != tgt+1)
            {
                ++numappears[std::abs(gdual[end])-1];
                ++end;
            }

            // Error if the target crossing appears less than twice.
            if(end >= gdual.size()) {
                ERR_MSG("ERROR: Every crossing must appear exactly twice.");
            }

            // Push setters that mark oddly appearing crossings
            for(std::size_t j = beg; j < end; ++j) {
                if(numappears[std::abs(gdual[j])-1]%2 != 0)
                    setter_q.emplace(std::abs(gdual[j])-1, -mrk);
            }
        }
    }

    // Pairing check on both "inside" and "outside" crossings.
    std::stack<int> pairing_out{}, pairing_inn{};
    for(auto c : gdual) {
        auto c_abs = std::abs(c);
        auto& pairing = marks[c_abs-1] == -1 ? pairing_out : pairing_inn;

        if(!pairing.empty() && pairing.top() == c_abs)
            pairing.pop();
        else
            pairing.push(c_abs);
    }
    if(!pairing_out.empty() || !pairing_inn.empty()) {
        ERR_MSG("ERROR: Pairing check failed.");
        return std::nullopt;
    }

    // Extract the "outside" crossings.
    gcode_t result{};
    for(int c = 0; static_cast<std::size_t>(c) < numcrs; ++c) {
        if (marks[c] == -1)
            result.push_back(c+1);
    }

    // Return the result.
    return result;
}


/**************************************************
 *** Implementation of read_gauss_code function ***
 **************************************************/
//! Use an enhanced version of L. Kauffman's algorithm to determine signs of crossings.
//! See:
//! Louis H. Kauffman. 1999. Virtual Knot Theory. Eur. J. Comb. 20, 7 (Oct. 1999), 663–691. DOI:https://doi.org/10.1006/eujc.1999.0314 \see{http://homepages.math.uic.edu/~kauffman/VKT.pdf}.
std::optional<LinkDiagram>
khover::read_gauss_code(
    gcode_t const& gcode,
    std::vector<std::pair<std::size_t,bool>> const& signs
    ) noexcept
{
    // Estimate the number of crossings.
    std::size_t numcrs = count_crossing(gcode);

    // Compute the dual Gauss code.
    // Note that it yields a Jordan curve (i.e. the embedding of S^1 in the plane) consisting of all the arcs of each disjoint component.
    auto gdual = compute_dual(numcrs, gcode);

    // Error occured.
    if(gdual.empty()) {
        std::cerr << "Failed to compute dual Gauss code." << std::endl;
        return std::nullopt;
    }

    auto outer = separate_vertices(numcrs, gdual);
    if(!outer) {
        std::cerr << "Failed to divide vertices." << std::endl;
        return std::nullopt;
    }

    // Determine the "counter-clockwise" ordering on edges adjacent to each vervex.
    std::vector<std::array<int,4>> eord(numcrs,{0,0,0,0});
    for(int v = 0; static_cast<std::size_t>(v) < numcrs; ++v) {
        // If v is a vertex outside the Jordan circle.
        if (std::find(std::begin(*outer), std::end(*outer), v+1)
            != std::end(*outer))
        {
            // We just believe itr != std::end(gdual).
            auto itr = std::find_if(
                std::begin(gdual), std::end(gdual),
                [v](auto n) { return std::abs(n) == v+1; });
            eord[v][3]
                = itr == std::begin(gdual)
                ? gdual.back()
                : *std::prev(itr);
            eord[v][2]
                = std::next(itr) == std::end(gdual)
                ? -(gdual.front())
                : -(*std::next(itr));

            itr = std::find_if(
                std::next(itr), std::end(gdual),
                [v](auto n) { return std::abs(n) == v+1; });
            eord[v][1]
                = itr == std::begin(gdual)
                ? gdual.back()
                : *std::prev(itr);
            eord[v][0]
                = std::next(itr) == std::end(gdual)
                ? -(gdual.front())
                : -(*std::next(itr));
        }
        // If v is a vertex inside the Jordan circle.
        else {
            // We just believe itr != std::end(gdual).
            auto itr = std::find_if(
                std::begin(gdual), std::end(gdual),
                [v](auto n) { return std::abs(n) == v+1; });
            eord[v][0]
                = itr == std::begin(gdual)
                ? gdual.back()
                : *std::prev(itr);
            eord[v][1]
                = std::next(itr) == std::end(gdual)
                ? -(gdual.front())
                : -(*std::next(itr));

            itr = std::find_if(
                std::next(itr), std::end(gdual),
                [v](auto n) { return std::abs(n) == v+1; });
            eord[v][2]
                = itr == std::begin(gdual)
                ? gdual.back()
                : *std::prev(itr);
            eord[v][3]
                = std::next(itr) == std::end(gdual)
                ? -(gdual.front())
                : -(*std::next(itr));
        }
    }

    // Create the list of crossings
    std::vector<LinkDiagram::Crossing> crosses(
        numcrs,
        LinkDiagram::Crossing{0, 0, 0, 0});

    std::size_t ixbase = 0;
    // The table remembering the source and the target of each edge.
    std::map<std::size_t,std::pair<int,int>> edge_map{};

    // Write the adjacent edges.
    // Here, edges are labeled by the indices of the source vertices in the gcode.
    for(auto comp : divide_gcode(gcode)) {
        for(std::size_t i = 0; i < comp.size(); ++i) {
            // The next vertex in the component.
            std::size_t inext = (i+1) % comp.size();
            // Set current edge as the out-going edge which is over/under according to the sign of the vertex.
            crosses[std::abs(comp[i])-1].adj_arc[comp[i] < 0 ? 3 : 1]
                = i + ixbase;
            // Set current edge as the in-coming edge which is over/under according to the sign of the next vertex.
            crosses[std::abs(comp[inext])-1].adj_arc[comp[inext] < 0 ? 2 : 0]
                = i + ixbase;
            edge_map.emplace(i+ixbase, std::make_pair(comp[i], comp[inext]));
        }
        ixbase += comp.size();
    }

    // The table of signs.
    LinkDiagram::signs_t crs_signs{};

    // Choose a sign on each crossing
    for(int v = 0; static_cast<std::size_t>(v) < numcrs; ++v) {
        auto orditr_o = std::find(
            std::begin(eord[v]),
            std::end(eord[v]),
            edge_map[crosses[v].adj_arc[1]].second);

        if (orditr_o == std::end(eord[v])) {
            ERR_MSG(edge_map[crosses[v].adj_arc[1]] << " is not an arc.");
            return std::nullopt;
        }

        auto orditr_u = std::find(
            std::begin(eord[v]),
            std::end(eord[v]),
            edge_map[crosses[v].adj_arc[3]].second);

        if (orditr_u == std::end(eord[v])) {
            ERR_MSG(edge_map[crosses[v].adj_arc[3]] << " is not an arc.");
            return std::nullopt;
        }

        crs_signs.set(
            v,
            (std::distance(orditr_o, orditr_u) + 4)%4 == 1
            );
    }

    // Take mirrors so that required signs are satisfied.
    // This step breaks gdual.
    auto disjcomps = divide_gcode(gdual);
    for(auto &sign : signs) {
        // Find the disjoint component containing the crossing in the sign request.
        auto compitr = std::find_if(
            std::begin(disjcomps), std::end(disjcomps),
            [&sign](gcode_t const& disjcomp) {
                for(auto c : disjcomp)
                    if (static_cast<std::size_t>(std::abs(c)) == sign.first)
                        return true;
                return false;
            } );

        // Component not found, or the current diagram already has the sign, skip that request.
        if (compitr == std::end(disjcomps)
            || crs_signs.test(sign.first-1) == sign.second)
            continue;

        // Remove duplicates in the disjoint component.
        std::sort(std::begin(*compitr), std::end(*compitr));
        auto lastitr = std::unique(
            std::begin(*compitr), std::end(*compitr));
        // Take the mirror image of the component.
        for(auto itr = std::begin(*compitr); itr != lastitr; ++itr) {
            crs_signs.flip(std::abs(*itr)-1);
        }
    }

    return LinkDiagram(2*numcrs, std::move(crosses), crs_signs);
}


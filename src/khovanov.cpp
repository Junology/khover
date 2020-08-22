/*!
 * \file khovanov.cpp
 * \author Jun Yoshida
 * \copyright (c) 2020 Jun Yoshida.
 * The project is released under the 2-clause BSD License.
 * \date August, 2020: created
 */

#include <vector>

#include <Eigen/Dense>

#include "khovanov.hpp"

/* Debug
#include "debug/debug.hpp"
//*/

using namespace khover;

using matrix_t = ChainIntegral::matrix_t;


/*******************************************
 *** Multiplication and comultiplication ***
 *******************************************/

//! Generate matrices for multiplications
//! By the definition of Khovanov complexes, the resulting matrices consist of 0 and 1.
//! \param numx The number of 'x's in enhancements on the codomain state. It turns out that the domain has to have the same number of 'x's.
//! \param ncomps_in_cod The number of connected components in the codomain state. Since this matrix represents the multiplication, the domain has to have exactly (ncomps_in_cod + 1) components.
//! \param comp_idx1 The common index of the components in the domain and the codomain states which are involved with the multiplication.
//! \param comp_idx2 The index of the other component in the domain involved with the multiplication.
//! \remark The function assumes idx1 < idx2.
//! \return The table of positions in the matrix where coefficients are 1.
static
std::vector<std::pair<matrix_t::Index,matrix_t::Index>>
matrix_mult(
    int numx, int ncomps_in_cod,
    std::size_t comp_idx1, std::size_t comp_idx2)
    noexcept
{
    auto enh_bound = binom(ncomps_in_cod, numx);
    std::vector<std::pair<matrix_t::Index,matrix_t::Index>> result{};

    auto dommask = low_window<max_components>(comp_idx2);

    for(std::size_t k = 0u; k < enh_bound; ++k) {
        auto enh_cod = bitsWithPop<max_components>(numx, k);
        auto enh_dom = (enh_cod & dommask) | ((enh_cod & ~dommask) << 1);

        // Codomain has 'x' at idx1.
        if(enh_cod.test(comp_idx1)) {
            result.emplace_back(k, bitsIndex(enh_dom));
            enh_dom.set(comp_idx1, false);
            enh_dom.set(comp_idx2, true);
            result.emplace_back(k, bitsIndex(enh_dom));
        }
        // Cocomain has '1' at idx1.
        else {
            result.emplace_back(k, bitsIndex(enh_dom));
        }
    }

    return result;
}

//! Generate matrices for comultiplications
//! By the definition of Khovanov complexes, the resulting matrices consist of 0 and 1.
//! \param numx The number of 'x's in enhancements on the codomain state. It turns out that the domain has to have exactly (numx-1) 'x's.
//! \param ncomps_in_cod The number of connected components in the codomain state. Since this matrix represents the comultiplication, the domain has to have exactly (ncomps_in_cod - 1) components.
//! \param comp_idx1 The common index of the components in the domain and the codomain states which are involved with the comultiplication.
//! \param comp_idx2 The index of the other component in the codomain involved with the multiplication.
//! \remark The function assumes idx1 < idx2.
//! \return The table of positions in the matrix where coefficients are 1.
static
std::vector<std::pair<int,int>>
matrix_comult(
    int numx, int ncomps_in_cod,
    std::size_t comp_idx1, std::size_t comp_idx2)
    noexcept
{
    auto enh_bound = binom(ncomps_in_cod-1, numx-1);
    std::vector<std::pair<int,int>> result{};

    auto codmask = low_window<max_components>(comp_idx2);

    for(std::size_t k = 0u; k < enh_bound; ++k) {
        auto enh_dom = bitsWithPop<max_components>(numx-1, k);
        auto enh_cod = (enh_dom & codmask) | ((enh_dom & ~codmask) << 1);

        // Domomain has 'x' at idx1.
        if(enh_dom.test(comp_idx1)) {
            enh_cod.set(comp_idx2, true);
            result.emplace_back(bitsIndex(enh_cod), k);
        }
        // Domain has '1' at idx1.
        else {
            enh_cod.set(comp_idx2, true);
            result.emplace_back(bitsIndex(enh_cod), k);
            enh_cod.set(comp_idx1, true);
            enh_cod.set(comp_idx2, false);
            result.emplace_back(bitsIndex(enh_cod), k);
        }
    }

    return result;
}

/***********************************
 *** Implementation of functions ***
 ***********************************/
// Compute Khovanov complex of a given link diagram.
std::optional<ChainIntegral>
khover::khChain(const LinkDiagram &diagram, int qdeg)
    noexcept
{
    std::size_t nstates = cipow(2,diagram.ncrosses());

    // The table of component data.
    std::vector<std::pair<std::size_t,std::vector<component_t>>> comps;

    // The number of 'x' in the enhancements.
    // An entry may be negative if there is no enhancement in the q-degree.
    std::vector<int> popcnts;

    // The table of the begining indices of states in enhanced state vectors.
    std::vector<matrix_t::Index> headidx;

    // Compute all smoothings and the number of 'x's.
    comps.reserve(nstates);
    popcnts.reserve(nstates);
    headidx.reserve(nstates);
    for(std::size_t i = 0; i < nstates; ++i) {
        // Compute smoothing.
        comps.push_back(
            diagram.smoothing(state_t{i}));

        // Compute degree.
        int deg = qdeg - diagram.writhe()
            - diagram.cohDegree(state_t{i});

        // Parity check
        if ((comps[i].first - deg)%2 != 0) {
            return std::nullopt;
        }

        // The number of 'x' in the enhancement of the degree.
        // If m is # of '1' and n is # of 'x', then
        // m + n = comps[i].first
        // m - n = deg
        popcnts.push_back((comps[i].first - deg)/2);

        // Compute the begining index
        if(auto sti = bitsIndex(state_t{i});
           sti <= 0)
        {
            headidx.push_back(0);
        }
        else {
            auto i_prev = bitsWithPop<max_components>(
                state_t{i}.count(), sti-1).to_ullong();
            headidx.push_back(
                headidx[i_prev]
                + binom(comps[i_prev].first, popcnts[i_prev]) );
        }
    }

    // for(std::size_t pcnt = 0; pcnt <= diagram.ncrosses(); ++pcnt) {
    //     DBG_MSG(
    //         "\e[32;1mCoh. degree = "
    //         << static_cast<int>(pcnt) - static_cast<int>(diagram.nnegative())
    //         << "\e[m");
    //     for(std::size_t i = 0; i < binom(diagram.ncrosses(), pcnt); ++i) {
    //         auto st = bitsWithPop<max_crosses>(pcnt, i).to_ulong();
    //         DBG_MSG(
    //             std::bitset<8>(st) << "\n"
    //             << comps[st] << "\n"
    //             << popcnts[st] << "\n"
    //             << headidx[st] );
    //     }
    // }

    // Compute the matrices representing matrices.
    ChainIntegral result(
        -diagram.npositive()-1,
        matrix_t(
            0,
            binom(comps[nstates-1].first, popcnts[nstates-1]))
        );
    for(int i = diagram.ncrosses(); i > 0; --i) {
        std::size_t maxst_cod
            = (low_window<max_crosses>(i)<<(diagram.ncrosses()-i)).to_ullong();
        std::size_t maxst_dom
            = (low_window<max_crosses>(i-1)<<(diagram.ncrosses()-i+1)).to_ullong();
        matrix_t diffmat = matrix_t::Zero(
            headidx[maxst_cod]
            + binom(comps[maxst_cod].first, popcnts[maxst_cod]),
            headidx[maxst_dom]
            + binom(comps[maxst_dom].first, popcnts[maxst_dom]));

        // Traverse all the state pairs
        for(std::size_t stidx_cod = 0;
            stidx_cod < binom(diagram.ncrosses(), i);
            ++stidx_cod)
        {
            // The state associated with the index.
            auto st_cod = bitsWithPop<max_crosses>(i, stidx_cod).to_ullong();

            for(std::size_t stidx_dom = 0;
                stidx_dom < binom(diagram.ncrosses(), i-1);
                ++stidx_dom)
            {
                // The state associated with the index.
                auto st_dom = bitsWithPop<max_crosses>(i-1, stidx_dom).to_ullong();

                // The coefficient between the domain/codomain states.
                auto coeff = diagram.stateCoeff(
                    state_t{st_dom}, state_t{st_cod});

                // Skip the case where the codomain state is not adjacent to the domain state.
                if(coeff == 0)
                    continue;

                // Find the position where saddle operation is applied.
                for(std::size_t arc = 0; arc < diagram.narcs(); ++arc) {
                    // The saddle causes multiplication.
                    if(comps[st_cod].second[arc] < comps[st_dom].second[arc]) {
                        for(auto [r,c] : matrix_mult(
                                popcnts[st_cod],
                                comps[st_cod].first,
                                comps[st_cod].second[arc],
                                comps[st_dom].second[arc])) {
                            diffmat.coeffRef(
                                headidx[st_cod]+r,
                                headidx[st_dom]+c
                                ) += coeff;
                        }
                        break;
                    }
                    // The saddle causes comultiplication
                    else if (comps[st_cod].second[arc] > comps[st_dom].second[arc]) {
                        for(auto [r,c] : matrix_comult(popcnts[st_cod], comps[st_cod].first, comps[st_dom].second[arc], comps[st_cod].second[arc])) {
                            diffmat.coeffRef(
                                headidx[st_cod]+r,
                                headidx[st_dom]+c
                                ) += coeff;
                        }
                        break;
                    }
                }
            }
        }

        // Append the matrix as a differential.
        if (!result.prepend(std::move(diffmat))) {
            std::cerr << "Wrong size matrix..." << std::endl;
            return std::nullopt;
        }
    }

    return std::optional<ChainIntegral>(std::move(result));
}

// Compute crux complex of a given link diagram and a given crossing.
std::optional<ChainIntegral>
khover::cruxChain(LinkDiagram diagram, std::size_t dblpt, int qdeg)
    noexcept
{
    // The double point is out-of-range.
    if(dblpt >= diagram.ncrosses()) {
        return std::nullopt;
    }

    // The number of states on crossings except the double point.
    std::size_t nstates = cipow(2,diagram.ncrosses()-1);

    struct CompT {
        matrix_t::Index headidx;
        state_t st_ext;
        std::size_t ncomp{0};
        int xcnt{-1};
        std::vector<component_t> arc_tbl{};
        std::bitset<max_arcs> twisted_tbl{0u};

        inline std::size_t nenhancement() const noexcept {
            return binom(ncomp,xcnt);
        }
    };

    std::vector<CompT> summand(nstates, CompT{});

    for(std::size_t st = 0; st < nstates; ++st) {
        auto stbits = state_t{st};

        // Compute the head index of the state in the complex.
        if (auto i = bitsIndex(stbits); i == 0) {
            summand[st].headidx = 0;
        }
        else {
            // Compute the previous state in the same degree.
            auto st_prev = bitsWithPop<max_crosses>(
                stbits.count(),
                bitsIndex(stbits)-1
                ).to_ulong();
            summand[st].headidx
                = summand[st_prev].headidx
                + binom(summand[st_prev].ncomp, summand[st_prev].xcnt);
        }

        // Get twisted arcs.
        // If the state is not crux, skip it.
        if (auto iscrux = diagram.cruxTwists(st,dblpt);
            !iscrux) {
            continue;
        }
        // For crux states, remember the twisted arcs and compute the smoothings.
        else {
            summand[st].st_ext = iscrux->first;
            summand[st].twisted_tbl = iscrux->second;
            std::tie(summand[st].ncomp, summand[st].arc_tbl)
                = diagram.smoothing(iscrux->first);
        }

        // Compute degree.
        // In this step, the double point is counted as a positive crossing.
        int deg
            = qdeg + 2*diagram.nnegative() - diagram.npositive()
            - stbits.count() + (diagram.crosses()[dblpt].is_positive ? 0 : -3);

        // Parity check
        if ((summand[st].ncomp - deg)%2 != 0) {
            return std::nullopt;
        }

        // The number of 'x' in the enhancement of the degree.
        // If m is # of '1' and n is # of 'x', then
        // m + n = ncomps
        // m - n = deg
        summand[st].xcnt = (summand[st].ncomp - deg)/2;
    }

    // Compute the matrices representing matrices.
    ChainIntegral result(
        -diagram.npositive() - (diagram.crosses()[dblpt].is_positive ? 0 : 1),
        matrix_t(0, summand[nstates-1].nenhancement())
        );

    for(int i = diagram.ncrosses()-1; i > 0; --i) {
        std::size_t maxst_cod
            = (low_window<max_crosses>(i)<<(diagram.ncrosses()-i-1)).to_ullong();
        std::size_t maxst_dom
            = (low_window<max_crosses>(i-1)<<(diagram.ncrosses()-i)).to_ullong();
        matrix_t diffmat = matrix_t::Zero(
            summand[maxst_cod].headidx + summand[maxst_cod].nenhancement(),
            summand[maxst_dom].headidx + summand[maxst_dom].nenhancement()
            );

        // Traverse all the state pairs
        for(std::size_t stidx_cod = 0;
            stidx_cod < binom(diagram.ncrosses()-1, i);
            ++stidx_cod)
        {
            // The state associated with the index.
            auto stbits_cod = bitsWithPop<max_crosses>(i, stidx_cod);
            auto st_cod = stbits_cod.to_ulong();

            // Skip states that has no enhancement in the q-degree.
            if(summand[st_cod].xcnt < 0
               || summand[st_cod].xcnt > static_cast<int>(summand[st_cod].ncomp))
            {
                continue;
            }

            for(std::size_t c = 0; c < diagram.ncrosses()-1; ++c) {
                if(!stbits_cod.test(c))
                    continue;

                auto stbits_dom = stbits_cod;
                stbits_dom.set(c,false);
                auto st_dom = stbits_dom.to_ulong();

                // Skip states that has no enhancement in the q-degree.
                if(summand[st_dom].xcnt < 0
                   || summand[st_dom].xcnt > static_cast<int>(summand[st_cod].ncomp))
                {
                    continue;
                }

                // The coefficient between the domain/codomain states.
                auto coeff = diagram.stateCoeff(
                    summand[st_dom].st_ext, summand[st_cod].st_ext);

                // Skip the case where the codomain state is not adjacent to the domain state.
                if(coeff == 0)
                    continue;

                // Find the position where saddle operation is applied.
                for(std::size_t arc = 0; arc < diagram.narcs(); ++arc) {
                    // The saddle causes multiplication.
                    if(summand[st_cod].arc_tbl[arc] < summand[st_dom].arc_tbl[arc])
                    {
                        // Enabled if the operation is involved with twisted arcs.
                        std::optional<std::size_t> action_arc;
                        if ((summand[st_dom].twisted_tbl
                             ^ summand[st_cod].twisted_tbl).any())
                        {
                            // In that case, action_arc is a non-twisted arc that acts on the twisted arc.
                            if(summand[st_dom].twisted_tbl.test(arc))
                                action_arc = summand[st_cod].arc_tbl[arc];
                            else
                                action_arc = summand[st_dom].arc_tbl[arc];
                        }

                        for(auto [r,c] : matrix_mult(
                                summand[st_cod].xcnt,
                                summand[st_cod].ncomp,
                                summand[st_cod].arc_tbl[arc],
                                summand[st_dom].arc_tbl[arc]))
                        {
                            // -1 if 'x' on the act circle.
                            diffmat.coeffRef(
                                summand[st_cod].headidx+r,
                                summand[st_dom].headidx+c
                                ) += action_arc && bitsWithPop<max_components>(
                                    summand[st_dom].xcnt, c).test(*action_arc)
                                ? -coeff : coeff;
                        }
                        break;
                    }
                    // The saddle causes comultiplication
                    else if (summand[st_cod].arc_tbl[arc]
                             > summand[st_dom].arc_tbl[arc])
                    {
                        // Enabled if the operation is involved with twisted arcs.
                        std::optional<std::size_t> coact_arc;
                        if ((summand[st_dom].twisted_tbl
                             ^ summand[st_cod].twisted_tbl).any())
                        {
                            // In that case, coact_arc is a non-twisted arc that coacts on the twisted arc.
                            if(summand[st_cod].twisted_tbl.test(arc))
                                coact_arc = summand[st_dom].arc_tbl[arc];
                            else
                                coact_arc = summand[st_cod].arc_tbl[arc];
                        }

                        for(auto [r,c] : matrix_comult(
                                summand[st_cod].xcnt,
                                summand[st_cod].ncomp,
                                summand[st_dom].arc_tbl[arc],
                                summand[st_cod].arc_tbl[arc]))
                        {
                            // -1 if '1' on the coact circle.
                            diffmat.coeffRef(
                                summand[st_cod].headidx+r,
                                summand[st_dom].headidx+c
                                ) += coact_arc && bitsWithPop<max_components>(
                                    summand[st_cod].xcnt, r).test(*coact_arc)
                                ? coeff : -coeff;
                        }
                        break;
                    }
                }
            }
        }

        // Append the matrix as a differential.
        if (!result.prepend(std::move(diffmat))) {
            std::cerr << "Wrong size matrix..." << std::endl;
            return std::nullopt;
        }
    }

    return std::optional<ChainIntegral>(std::move(result));
}

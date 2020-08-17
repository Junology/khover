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

using namespace khover;

using matrix_t = ChainIntegral::matrix_t;

/* Debug
template<class T>
std::ostream& operator<<(std::ostream& os, std::vector<T> const& vec) {
    if (vec.empty()) {
        os << "{}";
        return os;
    }

    os << "{" << vec.front();
    for(auto itr = std::next(std::begin(vec)); itr != std::end(vec); ++itr)
        os << ", " << *itr;
    os << "}";
    return os;
}

template <>
std::ostream& operator<<(std::ostream& os, std::vector<uint8_t> const& vec) {
    if (vec.empty()) {
        os << "{}";
        return os;
    }

    os << "{" << vec.front()+0;
    for(auto itr = std::next(std::begin(vec)); itr != std::end(vec); ++itr)
        os << ", " << (*itr)+0;
    os << "}";
    return os;
}
//*/

/*******************************************
 *** Multiplication and comultiplication ***
 *******************************************/

//! Generate matrices for multiplications
//! By the definition of Khovanov complexes, the resulting matrices consist of 0 and 1.
//! \param numx The number of 'x's in enhancements on the codomain state. It turns out that the domain has to have the same number of 'x's.
//! \param ncomps_in_cod The number of connected components in the codomain state. Since this matrix represents the multiplication, the domain has to have exactly (ncomps_in_cod + 1) components.
//! \param idx1 The common index of the components in the domain and the codomain states which are involved with the multiplication.
//! \param idx2 The index of the other component in the domain involved with the multiplication.
//! \remark The function assumes idx1 < idx2.
//! \return The table of positions in the matrix where coefficients are 1.
static
std::vector<std::pair<int,int>>
matrix_mult(
    int numx, int ncomps_in_cod,
    std::size_t idx1, std::size_t idx2)
    noexcept
{
    auto enh_bound = binom(ncomps_in_cod, numx);
    std::vector<std::pair<int,int>> result{};

    auto dommask = low_window<max_components>(idx2);

    for(std::size_t k = 0u; k < enh_bound; ++k) {
        auto enh_cod = bitsWithPop<max_components>(numx, k);
        auto enh_dom = (enh_cod & dommask) | ((enh_cod & ~dommask) << 1);

        // Codomain has 'x' at idx1.
        if(enh_cod.test(idx1)) {
            result.emplace_back(k, bitsIndex(enh_dom));
            enh_dom.set(idx1, false);
            enh_dom.set(idx2, true);
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
//! \param idx1 The common index of the components in the domain and the codomain states which are involved with the comultiplication.
//! \param idx2 The index of the other component in the codomain involved with the multiplication.
//! \remark The function assumes idx1 < idx2.
//! \return The table of positions in the matrix where coefficients are 1.
static
std::vector<std::pair<int,int>>
matrix_comult(
    int numx, int ncomps_in_cod,
    std::size_t idx1, std::size_t idx2)
    noexcept
{
    auto enh_bound = binom(ncomps_in_cod-1, numx-1);
    std::vector<std::pair<int,int>> result{};

    auto codmask = low_window<max_components>(idx2);

    for(std::size_t k = 0u; k < enh_bound; ++k) {
        auto enh_dom = bitsWithPop<max_components>(numx-1, k);
        auto enh_cod = (enh_dom & codmask) | ((enh_dom & ~codmask) << 1);

        // Domomain has 'x' at idx1.
        if(enh_dom.test(idx1)) {
            enh_cod.set(idx2, true);
            result.emplace_back(bitsIndex(enh_cod), k);
        }
        // Domain has '1' at idx1.
        else {
            enh_cod.set(idx2, true);
            result.emplace_back(bitsIndex(enh_cod), k);
            enh_cod.set(idx1, true);
            enh_cod.set(idx2, false);
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
            headidx.push_back(
                headidx[
                    bitsWithPop<max_components>(
                        state_t{i}.count(), sti-1).to_ullong()]
                + binom(comps[i].first, popcnts[i]) );
        }
    }

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
                        for(auto [r,c] : matrix_mult(popcnts[st_cod], comps[st_cod].first, comps[st_cod].second[arc], comps[st_dom].second[arc])) {
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

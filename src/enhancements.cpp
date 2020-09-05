/*!
 * \file states.hpp
 * \author Jun Yoshida
 * \copyright (c) 2020 Jun Yoshida.
 * The project is released under the 2-clause BSD License.
 * \date September, 2020: created
 */

#include "enhancements.hpp"

using namespace khover;

using matrix_t = ChainIntegral::matrix_t;

// Generate matrices for multiplications
// By the definition of Khovanov complexes, the resulting matrices consist of 0 and 1.
auto khover::matrix_mult(
    int numx, int ncomps_in_cod,
    std::size_t comp_idx1, std::size_t comp_idx2)
    noexcept
    -> std::vector<
        std::pair<matrix_t::Index,matrix_t::Index>
        >
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


// Generate matrices for comultiplications
// By the definition of Khovanov complexes, the resulting matrices consist of 0 and 1.
auto khover::matrix_comult(
    int numx, int ncomps_in_cod,
    std::size_t comp_idx1, std::size_t comp_idx2)
    noexcept
    -> std::vector<
        std::pair<matrix_t::Index,matrix_t::Index>
        >
{
    auto enh_bound = binom(ncomps_in_cod-1, numx-1);
    std::vector<std::pair<matrix_t::Index,matrix_t::Index>> result{};

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

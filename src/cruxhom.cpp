/*!
 * \file cruxhom.hpp
 * \author Jun Yoshida
 * \copyright (c) 2020 Jun Yoshida.
 * The project is released under the 2-clause BSD License.
 * \date August, 2020: created
 */

#include <thread>

#include "states.hpp"
#include "khovanov.hpp"
#include "enhancements.hpp"

#include "debug/debug.hpp"

using namespace khover;

using matrix_t = ChainIntegral::matrix_t;

//! Provide zero map for cubes
//! \param cube_dom The domain cube.
//! \param cube_cod The codomain cube.
//! \pre cube_dom.dim() == cube_cod.dim()
template <template <class...> class C1, template <class...> class C2>
static
std::vector<matrix_t>
zero_on_cubes(
    Cube<C1> const& cube_dom,
    std::vector<EnhancementProperty> enhprop_dom,
    Cube<C2> const& cube_cod,
    std::vector<EnhancementProperty> enhprop_cod
    ) noexcept
{
    auto mindeg = std::min(cube_dom.mincohdeg(), cube_cod.mincohdeg());
    auto maxdeg = std::max(cube_dom.maxcohdeg(), cube_cod.maxcohdeg());

    std::vector<matrix_t> homs{};
    for(auto i = maxdeg; i >= mindeg; --i) {
        std::size_t rk_cod, rk_dom;

        if(i >= cube_dom.mincohdeg() && i <= cube_dom.maxcohdeg()) {
            std::size_t maxst_dom = cube_dom.maxState(
                i-cube_dom.mincohdeg()).to_ulong();
            rk_dom = enhprop_dom[maxst_dom].headidx
                + binom(cube_dom[maxst_dom].ncomp, enhprop_dom[maxst_dom].xcnt);
        }
        else
            rk_dom = 0;

        if(i >= cube_cod.mincohdeg() && i <= cube_cod.maxcohdeg()) {
            std::size_t maxst_cod = cube_cod.maxState(
                i-cube_cod.mincohdeg()).to_ulong();
            rk_cod = enhprop_cod[maxst_cod].headidx
                + binom(cube_cod[maxst_cod].ncomp, enhprop_cod[maxst_cod].xcnt);
        }
        else
            rk_cod = 0;

        homs.push_back(matrix_t::Zero(rk_cod, rk_dom));
    }
    return homs;
}

/*****************
 *** Morphisms ***
 *****************/

//! The chain contratcion theta
template<template <class...> class C1, template <class...> class C2>
static
ChainIntegral::Hom
theta(
    LinkDiagram diagram,
    std::size_t dblpt,
    Cube<C1> const& cube_dom,
    std::vector<EnhancementProperty> enhprop_dom,
    Cube<C2> const& cube_cod,
    std::vector<EnhancementProperty> enhprop_cod
    ) noexcept
{
    std::vector<matrix_t> homs = zero_on_cubes(
        cube_dom, enhprop_dom, cube_cod, enhprop_cod);

    for(std::size_t st = 0; st < cube_dom.size(); ++st) {
        auto urcomp_dom = cube_dom[st].arccomp[diagram.getURArc(dblpt)];
        auto dlcomp_dom = cube_dom[st].arccomp[diagram.getDLArc(dblpt)];
        auto urcomp_cod = cube_cod[st].arccomp[diagram.getURArc(dblpt)];
        auto dlcomp_cod = cube_cod[st].arccomp[diagram.getDLArc(dblpt)];

        auto rk = binom(cube_dom[st].ncomp, enhprop_dom[st].xcnt);
        for(std::size_t e = 0; e < rk; ++e) {
            auto enh = std::make_optional<enhancement_t>(
                bitsWithPop<max_components>(
                    enhprop_dom[st].xcnt, e)
                );

            // eta \circ epsilon
            if(urcomp_dom != dlcomp_dom && urcomp_cod != dlcomp_cod) {
                if (enh->test(urcomp_dom))
                    enh->set(urcomp_dom,false);
                else
                    enh.reset();
            }
            // eta
            else if (urcomp_dom == dlcomp_dom && urcomp_cod != dlcomp_cod) {
                bool aux = enh->test(urcomp_dom);
                *enh = insertBit(
                    *enh,
                    urcomp_dom == urcomp_cod ? dlcomp_cod : urcomp_cod,
                    false);
                enh->set(dlcomp_cod, aux);
                enh->set(urcomp_cod, false);
            }
            // epsilon
            else if (urcomp_dom != dlcomp_dom && urcomp_cod == dlcomp_cod) {
                if (enh->test(urcomp_dom))
                    *enh = purgeBit(
                        *enh,
                        urcomp_dom == urcomp_cod ? dlcomp_cod : urcomp_cod
                        );
                else
                    enh.reset();
            }

            // Write
            if(enh) {
                auto e_cod = bitsIndex(*enh);
                homs[state_t{st}.count()].coeffRef(
                    enhprop_cod[st].headidx + e_cod,
                    enhprop_dom[st].headidx + e
                    ) = 1;
            }
        }
    }
    return ChainIntegral::Hom(-1,std::move(homs));
}

// Compute the morphism PhiHat.
std::optional<ChainIntegral::Hom>
khover::crossPhiHat(
    LinkDiagram const& diagram,
    std::size_t crossing,
    SmoothCube const& cube_neg,
    std::vector<EnhancementProperty> const& enhprop_neg,
    SmoothCube const& cube_pos,
    std::vector<EnhancementProperty> const& enhprop_pos
    ) noexcept
{
    std::vector<matrix_t> homs
        = zero_on_cubes(cube_neg, enhprop_neg, cube_pos, enhprop_pos);

    std::size_t nstates = 1u << (cube_neg.dim()-1);
    for(std::size_t st = 0; st < nstates; ++st) {
        auto stbit_neg = insertBit(state_t{st}, crossing, true);
        auto stbit_pos = insertBit(state_t{st}, crossing, false);
        std::size_t st_neg = stbit_neg.to_ulong();
        std::size_t st_pos = stbit_pos.to_ulong();

        auto urcomp_neg = cube_neg[st_neg].arccomp[diagram.getURArc(crossing)];
        auto dlcomp_neg = cube_neg[st_neg].arccomp[diagram.getDLArc(crossing)];

        // Compute the sign
        ChainIntegral::coeff_t sign = 1;
        for(std::size_t c = crossing+1; c < diagram.ncrosses(); ++c) {
            if((diagram.getSign(c) > 0) == stbit_neg.test(c))
                sign = -sign;
        }

        // Compute the coefficients of the matrix.
        auto rk = binom(cube_neg[st_neg].ncomp, enhprop_neg[st_neg].xcnt);
        for(std::size_t e = 0; e < rk; ++e) {
            auto enh = bitsWithPop<max_components>(enhprop_neg[st_neg].xcnt, e);

            // Single component
            if(urcomp_neg == dlcomp_neg) {
                if(!enh.test(urcomp_neg)) {
                    enh.set(urcomp_neg);
                    homs[cube_pos.dim()-stbit_pos.count()].coeffRef(
                        enhprop_pos[st_pos].headidx + bitsIndex(enh),
                        enhprop_neg[st_neg].headidx + e
                        ) = 0; // sign * 2;
                }
            }
            // Two components
            else {
                // Both '1'
                if(!enh.test(urcomp_neg) && !enh.test(dlcomp_neg)) {
                    // -1*x
                    enh.set(urcomp_neg);
                    homs[cube_pos.dim()-stbit_pos.count()].coeffRef(
                        enhprop_pos[st_pos].headidx + bitsIndex(enh),
                        enhprop_neg[st_neg].headidx + e
                        ) = -sign;
                    enh.set(urcomp_neg, false);
                    // x*1
                    enh.set(dlcomp_neg);
                    homs[cube_pos.dim()-stbit_pos.count()].coeffRef(
                        enhprop_pos[st_pos].headidx + bitsIndex(enh),
                        enhprop_neg[st_neg].headidx + e
                        ) = sign;
                }
                // '1' on left, 'x' on right
                else if(!enh.test(urcomp_neg)) {
                    enh.set(urcomp_neg);
                    enh.set(dlcomp_neg);
                    homs[cube_pos.dim()-stbit_pos.count()].coeffRef(
                        enhprop_pos[st_pos].headidx + bitsIndex(enh),
                        enhprop_neg[st_neg].headidx + e
                        ) = sign;
                }
                // 'x' on left, '1' on right
                else if(!enh.test(dlcomp_neg)) {
                    enh.set(urcomp_neg);
                    enh.set(dlcomp_neg);
                    homs[cube_pos.dim()-stbit_pos.count()].coeffRef(
                        enhprop_pos[st_pos].headidx + bitsIndex(enh),
                        enhprop_neg[st_neg].headidx + e
                        ) = -sign;
                }
            }
        }
    }

    return ChainIntegral::Hom(
        -cube_pos.maxcohdeg(),
        std::move(homs));
}

std::optional<ChainIntegral::Hom>
khover::cruxXi(
    LinkDiagram const& diagram,
    CruxCube const& cubeCrx,
    LinkDiagram const& diagV,
    SmoothCube const& cubeV,
    LinkDiagram const& diagW,
    SmoothCube const& cubeW,
    std::size_t dblpt,
    int qdeg
    ) noexcept
{
    std::optional<ChainIntegral> chV_neg, chV_pos, chW_neg, chW_pos;
    std::optional<std::vector<EnhancementProperty>> propC_dom, propC_cod;
    std::optional<std::vector<EnhancementProperty>> propV_neg, propV_pos;
    std::optional<std::vector<EnhancementProperty>> propW_neg, propW_pos;

    // Compute Khovanov homologies.
    {
        std::thread thC(
            [&diagram, &cubeCrx, &propC_dom, &propC_cod, dblpt, q=qdeg]() {
                propC_dom = get_enhancement_prop(
                    diagram, cubeCrx,
                    diagram.getSign(dblpt) > 0 ? q-2 : q-4);
                propC_cod = get_enhancement_prop(
                    diagram, cubeCrx,
                    diagram.getSign(dblpt) > 0 ? q+4 : q+2);
            });
        std::thread thV_neg(
            [&diagV, &cubeV, &chV_neg, &propV_neg, q=qdeg]() {
                propV_neg = get_enhancement_prop(diagV, cubeV, q+1);
                if(propV_neg)
                    chV_neg = khChain(diagV, cubeV, *propV_neg);
            });
        std::thread thV_pos(
            [&diagV, &cubeV, &chV_pos, &propV_pos, q=qdeg]() {
                propV_pos = get_enhancement_prop(diagV, cubeV, q-1);
                if(propV_pos)
                    chV_pos = khChain(diagV, cubeV, *propV_pos);
            });
        std::thread thW_neg(
            [&diagW, &cubeW, &chW_neg, &propW_neg, q=qdeg]() {
                propW_neg = get_enhancement_prop(diagW, cubeW, q+2);
                if(propW_neg)
                    chW_neg = khChain(diagW, cubeW, *propW_neg);
            });
        std::thread thW_pos(
            [&diagW, &cubeW, &chW_pos, &propW_pos, q=qdeg]() {
                propW_pos = get_enhancement_prop(diagW, cubeW, q-2);
                if(propW_pos)
                    chW_pos = khChain(diagW, cubeW, *propW_pos);
            });

        thC.join();
        thV_neg.join();
        thV_pos.join();
        thW_neg.join();
        thW_pos.join();
    }

    if(!propC_dom || !propC_cod
       || !propV_neg || !chV_neg || !propV_pos || !chV_pos
       || !propW_neg || !chW_neg || !propW_pos || !chW_pos
        )
    {
        return std::nullopt;
    }

    ChainIntegral::Hom result = theta(
        diagram, dblpt,
        cubeW, *propW_neg,
        cubeCrx, *propC_cod) << 3;

    result *= chW_neg->diffAsHom() << 2;
    result *= theta(diagram, dblpt, cubeV, *propV_neg, cubeW, *propW_neg) << 2;
    result *= chV_neg->diffAsHom() << 1;
    result *= theta(diagram, dblpt, cubeV, *propV_pos, cubeV, *propV_neg) << 1;
    result *= chV_pos->diffAsHom();
    result *= theta(diagram, dblpt, cubeW, *propW_pos, cubeV, *propV_pos);
    result *= chW_pos->diffAsHom() << -1;
    result *= theta(diagram, dblpt, cubeCrx, *propC_dom, cubeW, *propW_pos) << -1;

    return std::make_optional(std::move(result));
}

#include <iostream>
#include <iomanip>

#include "linkdiagram.hpp"
#include "khovanov.hpp"

#include "debug/debug.hpp"
using namespace khover;

int main (int argc, char* argv[])
{
    auto mir8_19 = read_gauss_code(
        {1,-4,5,-6,-8,7,4,-5,6,-3,2,-1,-7,8,3,-2}, {std::make_pair(1,false)});

    if(!mir8_19) {
        ERR_MSG("Failed to load the knot.");
        return EXIT_FAILURE;
    }

    mir8_19->makeNegative(6);
    auto diag6neg = *mir8_19;
    mir8_19->makePositive(6);
    auto diag6pos = *mir8_19;

    // Make cubes
    auto cube_neg = SmoothCube::fromDiagram(diag6neg);
    mir8_19->makePositive(6);
    auto cube_pos = SmoothCube::fromDiagram(diag6pos);

    // Compute Khovanov homologies with q=-13
    // Compute enhancement properties
    auto enhprop_neg = get_enhancement_prop(diag6neg, cube_neg, -13);
    auto enhprop_pos = get_enhancement_prop(diag6pos, cube_pos, -13);

    if (!enhprop_neg || !enhprop_pos) {
        ERR_MSG("Failed to compute enhancement properties.");
        return EXIT_FAILURE;
    }

    // Compute chain complexes and PhiHat
    auto ch_neg = khChain(diag6neg, cube_neg, *enhprop_neg);
    auto ch_pos = khChain(diag6pos, cube_pos, *enhprop_pos);
    auto phihat = khover::crossPhiHat(
                *mir8_19, 6,
                cube_neg, *enhprop_neg,
                cube_pos, *enhprop_pos);

    if (!ch_neg || !ch_pos) {
        ERR_MSG("Failed to construct chain complexes");
        return EXIT_FAILURE;
    }

    if (!phihat) {
        ERR_MSG("Failed to compute PhiHat.");
        return EXIT_FAILURE;
    }

    // Check if Phihat is a chain map.
    for(int i = 9; i > 0; --i) {
        auto comm = ch_pos->getDiff(i) * phihat->morphism(i+1) - phihat->morphism(i) * ch_neg->getDiff(i);
        if(!comm.isZero()) {
            ERR_MSG("Phihat does not commute with differentials at i=" << i);
            DBG_MSG("Diff of C(+):\n" << ch_pos->getDiff(i));
            DBG_MSG("Diff of C(-):\n" << ch_neg->getDiff(i));
            DBG_MSG("Phihat@i+1\n" << phihat->morphism(i+1));
            DBG_MSG("Phihat@i\n" << phihat->morphism(i));
            return EXIT_FAILURE;
        }
    }

    // Compute the derivative as the mapping cone of PhiHat.
    auto ch_der = ChainIntegral::cone(*phihat, *ch_neg, *ch_pos);
    if(!ch_der) {
        ERR_MSG("Failed to compute mapping cone at q=" << -13);
        return EXIT_FAILURE;
    }

    for(auto cplx : {*ch_der, *ch_neg, *ch_pos}) {
        for(int i = 7; i > 5; --i) {
            DBG_MSG("rank C[" << i << "] = " << cplx.rank(i));
            DBG_MSG(cplx.getDiff(i));
        }
    }

    // Check if ch_der is a chain complex;
    for(int i = 9; i > 0; --i) {
        auto dd = ch_der->getDiff(i)*ch_der->getDiff(i+1);
        if(!dd.isZero()) {
            ERR_MSG("ch_der is not a chain complex; d[i]d[i+1]!=0 at i=" << i);
            return EXIT_FAILURE;
        }
    }

    /*
    // Show the result
    for(auto cplx : {*ch_der, *ch_neg, *ch_pos}) {
        auto h = cplx.compute();

        for (auto i=cplx.mindeg(); i <= cplx.maxdeg(); ++i) {
            std::cout << std::setw(4) << std::right << (-i) << ": "
                      << std::resetiosflags(std::ios_base::adjustfield | std::ios_base::basefield)
                      << h[i-cplx.mindeg()].compute().pretty()
                      << std::endl;
        }
    }
    // */
    return EXIT_SUCCESS;
}

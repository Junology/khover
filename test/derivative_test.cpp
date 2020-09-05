#include <iostream>

#include "linkdiagram.hpp"
#include "khovanov.hpp"

using namespace khover;

int main(int, char**)
{
    // trefoil with double point.
    {
        auto trefoil = read_gauss_code(
            {1, -3, 2, -1, 3, -2}, {std::make_pair(1,false)});

        if(!trefoil) {
            ERR_MSG("Failed to load the trefoil.");
            return EXIT_FAILURE;
        }

        auto cube_neg = SmoothCube::fromDiagram(*trefoil);
        trefoil->makePositive(0);
        auto cube_pos = SmoothCube::fromDiagram(*trefoil);

        for(int q=-9; q <= 1; q+=2) {
            // Compute the Khovanov homology with the 0-th crossing negative.
            trefoil->makeNegative(0);
            auto enhprop_neg = get_enhancement_prop(*trefoil, cube_neg, q);
            auto ch_neg = enhprop_neg
                ? khChain(*trefoil, cube_neg, *enhprop_neg)
                : std::nullopt;

            if(q%2 != 0 && !ch_neg) {
                ERR_MSG("Failed to compute the Khovanov complex for q=" << q);
                return EXIT_FAILURE;
            }

            // Compute the Khovanov homology with the 0-th crossing positive.
            trefoil->makePositive(0);
            auto enhprop_pos = get_enhancement_prop(*trefoil, cube_pos, q);
            auto ch_pos = enhprop_pos
                ? khChain(*trefoil, cube_pos, *enhprop_pos)
                : std::nullopt;

            if(!ch_pos) {
                ERR_MSG("Failed to compute the Khovanov comples for q=" << q);
                return EXIT_FAILURE;
            }

            // Compute PhiHat.
            auto phihat = khover::crossPhiHat(
                *trefoil, 0, cube_neg, *enhprop_neg, cube_pos, *enhprop_pos);
            if(!phihat) {
                ERR_MSG("Failed to compute the map PhiHat for q=" << q);
                return EXIT_FAILURE;
            }

            // Compute the derivative
            auto derch = ChainIntegral::cone(*phihat, *ch_neg, *ch_pos);
            if(!derch) {
                ERR_MSG("Failed to compute the derivative for q=" << q);
                return EXIT_FAILURE;
            }

            auto h = derch->compute();
            for(auto i=-derch->maxdeg(); i >= -derch->mindeg(); ++i) {
                if(auto khder = h[-i-derch->mindeg()].compute();
                   (q==-9 && i==-4 && (khder.freerank != 1
                                      || !khder.torsions.empty()))
                   || (q==-9 && i!=-4 && (khder.freerank != 0
                                          || !khder.torsions.empty()))
                   || (q==-7 && i==-3 && (khder.freerank != 0
                                          || khder.torsions.size()!= 1
                                          || khder.torsions[0] != 2))
                   || (q==-7 && i!=-3 && (khder.freerank != 0
                                          || !khder.torsions.empty()))
                   || (q==-5 && i==-3 && (khder.freerank != 1
                                          || !khder.torsions.empty()))
                   || (q==-5 && i!=-3 && (khder.freerank != 0
                                          || !khder.torsions.empty()))
                   || (q==-3 && i==-1 && (khder.freerank != 1
                                          || !khder.torsions.empty()))
                   || (q==-3 && i!=-1 && (khder.freerank != 0
                                          || !khder.torsions.empty()))
                   || (q==-1 && i==0 && (khder.freerank != 0
                                         || khder.torsions.size() != 1
                                         || khder.torsions[0] != 2))
                   || (q==-1 && i!=0 && (khder.freerank != 0
                                         || !khder.torsions.empty()))
                   || (q==1 && i==0 && (khder.freerank != 1
                                        || !khder.torsions.empty()))
                   || (q==1 && i!=0 && (khder.freerank != 0
                                        || !khder.torsions.empty())))
                {
                    ERR_MSG("Wrong derivative at (i,q)=" << std::make_pair(i,q)
                            << ":\n"
                            << khder.pretty());
                    return EXIT_FAILURE;
                }
            }
        }
    }

    return EXIT_SUCCESS;
}

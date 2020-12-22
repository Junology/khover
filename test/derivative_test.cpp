#include <iostream>

#include "linkdiagram.hpp"
#include "khovanov.hpp"

using namespace khover;

int main(int, char**)
{
    // trefoil with double point at 1.
    {
        std::cout << "\e[34;1m---\nTrefoil with double point at 1\n---\e[m" << std::endl;
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
        std::cout << "Passed" << std::endl;
    }

    // trefoil with double point at 2.
    {
        std::cout << "\e[34;1m---\nTrefoil with double point at 2\n---\e[m" << std::endl;
        auto trefoil = read_gauss_code(
            {1, -3, 2, -1, 3, -2}, {std::make_pair(1,false)});

        if(!trefoil) {
            ERR_MSG("Failed to load the trefoil.");
            return EXIT_FAILURE;
        }

        auto cube_neg = SmoothCube::fromDiagram(*trefoil);
        trefoil->makePositive(1);
        auto cube_pos = SmoothCube::fromDiagram(*trefoil);

        for(int q=-9; q <= 1; q+=2) {
            // Compute the Khovanov homology with the 1-st crossing negative.
            trefoil->makeNegative(1);
            auto enhprop_neg = get_enhancement_prop(*trefoil, cube_neg, q);
            auto ch_neg = enhprop_neg
                ? khChain(*trefoil, cube_neg, *enhprop_neg)
                : std::nullopt;

            if(q%2 != 0 && !ch_neg) {
                ERR_MSG("Failed to compute the Khovanov complex for q=" << q);
                return EXIT_FAILURE;
            }

            // Compute the Khovanov homology with the 1-st crossing positive.
            trefoil->makePositive(1);
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
                *trefoil, 1, cube_neg, *enhprop_neg, cube_pos, *enhprop_pos);
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
        std::cout << "Passed" << std::endl;
    }

    // 4-twist knot vs 2-twist knot (aka. figure-eight)
    {
        std::cout << "\e[34;1m---\n4-twist knot with double pint at 3\n---\e[m" << std::endl;
        auto twist4 = read_gauss_code(
            {1, -3, 4, -5, 6, -1, 2, -6, 5, -4, 3, -2}, {std::make_pair(1,true)});

        if(!twist4) {
            ERR_MSG("Failed to load twist4.");
            return EXIT_FAILURE;
        }

        auto twist2 = twist4;
        twist2->crossingChange(2);

        auto cube4 = SmoothCube::fromDiagram(*twist4);
        auto cube2 = SmoothCube::fromDiagram(*twist2);

        for(int q=-9; q <= 5; q+=2) {
            auto enhprop4 = get_enhancement_prop(*twist4, cube4, q);
            auto enhprop2 = get_enhancement_prop(*twist2, cube2, q);

            if(!enhprop4 || !enhprop2) {
                ERR_MSG("Failed to compute enhancements in q=" << q);
                return EXIT_FAILURE;
            }

            // Compute Khovanov complexes.
            auto ch4 = khChain(*twist4, cube4, *enhprop4);
            auto ch2 = khChain(*twist2, cube2, *enhprop2);

            // Compute PhiHat.
            auto phihat = crossPhiHat(
                *twist2, 2, cube4, *enhprop4, cube2, *enhprop2);
            if(!phihat) {
                ERR_MSG("Failed to compute the map PhiHat for q=" << q);
                return EXIT_FAILURE;
            }

            // Compute the derivative
            auto derch = ChainIntegral::cone(*phihat, *ch4, *ch2);
            if(!derch) {
                ERR_MSG("Failed to compute the derivative for q=" << q);
                return EXIT_FAILURE;
            }

            auto h = derch->compute();
            for(auto i = -derch->maxdeg(); i >= -derch->mindeg(); ++i) {
                if(auto khder = h[-i-derch->mindeg()].compute();
                   (q==-9 && i==-5 && (khder.freerank != 1 || !khder.isTorFree()))
                   || (q==-9 && i!=-5 && !khder.isZero())
                   || (q==-7 && i==-4 && (!khder.isFinite()
                                          || khder.torsions.size() != 1
                                          || khder.torsions[0] != 2))
                   || (q==-7 && i!=-4 && !khder.isZero())
                   || (q==-5 && i==-4 && (khder.freerank != 1
                                          || !khder.isTorFree()))
                   || (q==-5 && i!=-4 && !khder.isZero())
                   || (q==-3 && i==-2 && (khder.freerank != 1
                                          || !khder.isTorFree()))
                   || (q==-3 && i!=-2 && !khder.isZero())
                   || (q==-1 && i==-1 && (!khder.isFinite()
                                         || khder.torsions.size() != 1
                                         || khder.torsions[0] != 2))
                   || (q==-1 && i!=-1 && !khder.isZero())
                   || (q==1 && i==-1 && (khder.freerank != 1
                                         || !khder.isTorFree()))
                   || (q==1 && i!=-1 && !khder.isZero())
                   || (q>1 && !khder.isZero())
                    )
                {
                    ERR_MSG("Wrong first derivatives at (i,q)==" << std::make_pair(i,q));
                    return EXIT_FAILURE;
                }
            }
        }
        std::cout << "Passed" << std::endl;
    }

    return EXIT_SUCCESS;
}

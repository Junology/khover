#include <iostream>

#include "khovanov.hpp"
#include "debug/debug.hpp"

using namespace khover;

int main(int argc, char* argv[])
{
    // Crux of figure-eight
    {
        std::cout << "\e[34;1m---\nFigure-eight knot\n---\e[m" << std::endl;

        auto figureeight = read_gauss_code(
            {1,-3,4,-1,2,-4,3,-2},
            {std::make_pair(3,false)}
            );

        if(!figureeight) {
            ERR_MSG("Failed to load the figure-eight knot.");
            return EXIT_FAILURE;
        }
        int qmin =
            - static_cast<int>(figureeight->smoothing(0u).first)
            - static_cast<int>(figureeight->nnegative())
            + figureeight->writhe();
        int qmax =
            static_cast<int>(figureeight->smoothing(~state_t(0u)).first)
            + static_cast<int>(figureeight->npositive())
            + figureeight->writhe();

        for(int q = qmin; q <= qmax; q+=2) {
            std::cout << "\e[33;1mq-degree: " << q << "\e[m" <<std::endl;
            auto cch = cruxChain(*figureeight, 2, q);

            if(!cch) {
                ERR_MSG("Error in computing the crux complex.");
                return EXIT_FAILURE;
            }

            int i = -cch->mindeg();
            for(auto h : cch->compute()) {
                auto kh = h.compute();
                if(q==-1 && i == -1) {
                    if(kh.freerank != 1 || !kh.torsions.empty()) {
                        ERR_MSG("Wrong Kh^{-1,-1}(figure-eight): " << kh.pretty());
                        return EXIT_FAILURE;
                    }
                }
                else if(q==1 && i == 0) {
                    if(kh.freerank != 0
                       || kh.torsions.size() != 1 || kh.torsions[0] != 2) {
                        ERR_MSG("Wrong Kh^{0,1}(figure-eight): " << kh.pretty());
                        return EXIT_FAILURE;
                    }
                }
                else if(q==3 && i == 0) {
                    if(kh.freerank != 1 || !kh.torsions.empty()) {
                        ERR_MSG("Wrong Kh^{0,3}(figure-eight): " << kh.pretty());
                        return EXIT_FAILURE;
                    }
                }
                else {
                    if(kh.freerank != 0 || !kh.torsions.empty()) {
                        ERR_MSG("Non-trivial Kh in wrong positions: (i,q)="
                                << std::make_pair(i,q));
                        return EXIT_FAILURE;
                    }
                }
                --i;
            }
            std::cout << "Passed" << std::endl;
        }
    }

    // Crux of mirror of figure-eight
    {
        std::cout << "\e[34;1m---\nFigure-eight knot\n---\e[m" << std::endl;

        auto mfigureeight = read_gauss_code(
            {1,-3,4,-1,2,-4,3,-2},
            {std::make_pair(3,true)}
            );

        if(!mfigureeight) {
            ERR_MSG("Failed to load the figure-eight knot.");
            return EXIT_FAILURE;
        }
        int qmin =
            - static_cast<int>(mfigureeight->smoothing(0u).first)
            - static_cast<int>(mfigureeight->nnegative())
            + mfigureeight->writhe();
        int qmax =
            static_cast<int>(mfigureeight->smoothing(~state_t(0u)).first)
            + static_cast<int>(mfigureeight->npositive())
            + mfigureeight->writhe();

        for(int q = qmin; q <= qmax; q+=2) {
            std::cout << "\e[33;1mq-degree: " << q << "\e[m" <<std::endl;
            auto cch = cruxChain(*mfigureeight, 2, q);

            if(!cch) {
                ERR_MSG("Error in computing the crux complex.");
                return EXIT_FAILURE;
            }

            int i = -cch->mindeg();
            for(auto h : cch->compute()) {
                auto kh = h.compute();
                if(q==-1 && i == 0) {
                    if(kh.freerank != 1 || !kh.torsions.empty()) {
                        ERR_MSG("Wrong Kh^{0,-1}(mfigure-eight): " << kh.pretty());
                        return EXIT_FAILURE;
                    }
                }
                else if(q==1 && i == 1) {
                    if(kh.freerank != 0
                       || kh.torsions.size() != 1 || kh.torsions[0] != 2) {
                        ERR_MSG("Wrong Kh^{1,1}(mfigure-eight): " << kh.pretty());
                        return EXIT_FAILURE;
                    }
                }
                else if(q==3 && i == 1) {
                    if(kh.freerank != 1 || !kh.torsions.empty()) {
                        ERR_MSG("Wrong Kh^{0,3}(mfigure-eight): " << kh.pretty());
                        return EXIT_FAILURE;
                    }
                }
                else {
                    if(kh.freerank != 0 || !kh.torsions.empty()) {
                        ERR_MSG("Non-trivial Kh in wrong positions: (i,q)="
                                << std::make_pair(i,q));
                        return EXIT_FAILURE;
                    }
                }
                --i;
            }
            std::cout << "Passed" << std::endl;
        }
    }

    // Crux of the mirror of true lover's knot
    {
        std::cout << "\e[34;1m---\nTrue-lover's knot\n---\e[m" << std::endl;
        auto truelover = read_gauss_code(
            {1,-2,-3,4,-5,6,-7,8,2,-1,-6,7,-8,3,-4,5},
            {std::make_pair(2,false)});

        if(!truelover) {
            ERR_MSG("Failed to load the true-lover's knot.");
            return EXIT_FAILURE;
        }

        int qmin =
            - static_cast<int>(truelover->smoothing(0u).first)
            - static_cast<int>(truelover->nnegative())
            + truelover->writhe();
        int qmax =
            static_cast<int>(truelover->smoothing(~state_t(0u)).first)
            + static_cast<int>(truelover->npositive())
            + truelover->writhe();

        for(int q = qmin; q <= qmax; q+=2) {
            std::cout << "\e[33;1mq-degree: " << q << "\e[m" <<std::endl;
            auto cch = cruxChain(*truelover, 1, q);

            if(!cch) {
                ERR_MSG("Error in computing the crux complex.");
                return EXIT_FAILURE;
            }

            int i = -cch->mindeg();
            for(auto h : cch->compute()) {
                auto kh = h.compute();
                if(q==-13 && i==-4) {
                    if(kh.freerank != 2 || !kh.torsions.empty()) {
                        ERR_MSG("Wrong Kh^{-4,-13}(mir. of true-lovers): "
                                << kh.pretty());
                        return EXIT_FAILURE;
                    }
                }
                else if(q==-11 && i==-3) {
                    if(kh.freerank != 0
                       || kh.torsions.size() != 2
                       || kh.torsions[0] != 2
                       || kh.torsions[1] != 2) {
                        ERR_MSG("Wrong Kh^{-3,-11}(mir. of true-lovers): "
                                << kh.pretty());
                        return EXIT_FAILURE;
                    }
                }
                else if(q==-9 && i==-3) {
                    if(kh.freerank != 2 || !kh.torsions.empty()) {
                        ERR_MSG("Wrong Kh^{-3,-9}(mir. of true-lovers): "
                                << kh.pretty());
                        return EXIT_FAILURE;
                    }
                }
                else if(q==-9 && i==-2) {
                    if(kh.freerank != 1 || !kh.torsions.empty()) {
                        ERR_MSG("Wrong Kh^{-2,-9}(mir. of true-lovers): "
                                << kh.pretty());
                        return EXIT_FAILURE;
                    }
                }
                else if(q==-7 && i==-1) {
                    if(kh.freerank != 0
                       || kh.torsions.size() != 1 || kh.torsions[0] != 2) {
                        ERR_MSG("Wrong Kh^{-1,-7}(mir. of true-lovers): "
                                << kh.pretty());
                        return EXIT_FAILURE;
                    }
                }
                else if(q==-5 && i==-1) {
                    if(kh.freerank != 1 || !kh.torsions.empty()) {
                        ERR_MSG("Wrong Kh^{-1,-5}(mir. of true-lovers): "
                                << kh.pretty());
                        return EXIT_FAILURE;
                    }
                }
                else {
                    if(kh.freerank != 0 || !kh.torsions.empty()) {
                        ERR_MSG("Non-trivial Kh in wrong positions: (i,q)="
                                << std::make_pair(i,q)
                                << kh.pretty());
                        return EXIT_FAILURE;
                    }
                }
                --i;
            }
            std::cout << "Passed" << std::endl;
        }
    }

    return EXIT_SUCCESS;
}

#include <iostream>

#include "khovanov.hpp"

#include "debug/debug.hpp"

using namespace khover;

int main (int argc, char* argv[])
{
    // Trefoil knot
    {
        auto trefoil = read_gauss_code(
            {1, -3, 2, -1, 3, -2}, {std::make_pair(1,false)});
        auto trefoilcube = SmoothCube::fromDiagram(*trefoil);

        for(int q = -9; q < 0 ; ++q) {
            auto enh_prop = get_enhancement_prop(*trefoil, trefoilcube, q);
            auto ch = enh_prop
                ? khover::khChain(*trefoil, trefoilcube, *enh_prop)
                : std::nullopt;

            if(q%2 == 0) {
                if (ch) {
                    ERR_MSG("Non-trivial complex in illegal quantum degrees.");
                    return EXIT_FAILURE;
                }
                continue;
            }

            if(!ch) {
                ERR_MSG("Failed to compute Khovanov homology of the trefiol");
                return EXIT_FAILURE;
            }

            auto h = ch->compute();

            if(h.size() != 5) {
                ERR_MSG("Wrong number of cohomology groups.");
                return EXIT_FAILURE;
            }

            if(auto h0 = h[0].compute();
               h0.freerank != 0 || !h0.torsions.empty())
            {
                ERR_MSG(
                    "Incorrect 1st cohomology groups: (deg,qdeg)=("
                    << -ch->mindeg() << "," << q
                    << ")");
                return EXIT_FAILURE;
            }
            if(auto h1 = h[1].compute();
               ((q==-1 || q==-3) && (h1.freerank != 1 || !h1.torsions.empty()))
               || (q!=-1 && q!=-3 && (h1.freerank != 0 || !h1.torsions.empty())))
            {
                ERR_MSG(
                    "Incorrect 0th cohomology groups: (deg,qdeg)=("
                    << -ch->mindeg()-1 << "," << q
                    << ")");
                return EXIT_FAILURE;
            }
            if(auto h2 = h[2].compute();
               h2.freerank != 0 || !h2.torsions.empty())
            {
                ERR_MSG(
                    "Incorrect -1st cohomology groups: (deg,qdeg)=("
                    << -ch->mindeg()-2 << "," << q
                    << ")");
                return EXIT_FAILURE;
            }
            if(auto h3 = h[3].compute();
               (q==-5 && (h3.freerank != 1 || !h3.torsions.empty()))
               || (q==-7
                   && (h3.freerank != 0
                       || h3.torsions.size()!=1
                       || h3.torsions[0]!=2))
               || (q!=-5 && q!=-7 && (h3.freerank != 0 || !h3.torsions.empty())))
            {
                ERR_MSG(
                    "Incorrect -2nd cohomology groups: (deg,qdeg)=("
                    << -ch->mindeg()-3 << "," << q
                    << ")\n"
                    << h3.pretty()
                    );
                return EXIT_FAILURE;
            }
            if(auto h4 = h[4].compute();
               (q==-9 && (h4.freerank != 1 || !h4.torsions.empty()))
               || (q!=-9 && (h4.freerank != 0 || !h4.torsions.empty())))
            {
                ERR_MSG(
                    "Incorrect -3rd cohomology groups: (deg,qdeg)=("
                    << -ch->mindeg()-4 << "," << q
                    << ")");
                return EXIT_FAILURE;
            }
        }
    }

    // Unknot obtained by change the 0-th crossing in the trefoil.
    {
        auto trefoil = read_gauss_code(
            {1, -3, 2, -1, 3, -2}, {std::make_pair(1,false)});
        trefoil->makePositive(0);
        auto trefoilcube = SmoothCube::fromDiagram(*trefoil);

        for(int q = -9; q <= 1 ; ++q) {
            auto enh_prop = get_enhancement_prop(*trefoil, trefoilcube, q);
            auto ch = enh_prop
                ? khover::khChain(*trefoil, trefoilcube, *enh_prop)
                : std::nullopt;

            if(q%2 == 0) {
                if (ch) {
                    ERR_MSG("Non-trivial complex in q=" << q);
                    return EXIT_FAILURE;
                }
                continue;
            }

            if(!ch) {
                ERR_MSG("Failed to compute Khovanov homology of the unknot at q=" << q);
                return EXIT_FAILURE;
            }

            auto h = ch->compute();
            for(int i = -ch->maxdeg(); i <= -ch->mindeg(); ++i) {
                if(auto kh = h[-i-ch->mindeg()].compute();
                   (std::abs(q)==1 && i == 0
                    && (kh.freerank != 1 || !kh.torsions.empty()))
                   || ((std::abs(q)!=1 || i!=0)
                       && (kh.freerank != 0 || !kh.torsions.empty())))
                {
                    ERR_MSG("Wrong homology at (q,i)="
                            << std::make_pair(q,i) << ":\n"
                            << kh.pretty());
                    return EXIT_FAILURE;
                }
            }
        }
    }

    // Figure-eight knot
    {
        auto fig8 = read_gauss_code(
            {1,-3,4,-1,2,-4,3,-2},
            {std::make_pair(3,false)}
            );

        if(!fig8) {
            ERR_MSG("Failed to load the figure-eight knot.");
            return EXIT_FAILURE;
        }

        auto fig8cube = SmoothCube::fromDiagram(*fig8);

        for(int q = -5; q <= 5; ++q) {
            auto enh_prop = get_enhancement_prop(*fig8, fig8cube, q);
            auto ch = enh_prop
                ? khChain(*fig8, fig8cube, *enh_prop)
                : std::nullopt;
            if((q%2==0 && ch) || (q%2!=0 && !ch)) {
                ERR_MSG("Illegal cohomology: q=" << q);
                return EXIT_FAILURE;
            }
            else if (q%2==0)
                continue;

            auto h = ch->compute();

            for(int i = -ch->maxdeg(); i <= -ch->mindeg(); ++i) {
                if(auto kh = h[-i-ch->mindeg()].compute();
                   (i==-2 && q==-5 && (kh.freerank != 1 ||
                                       !kh.torsions.empty()))
                   || (i==-2 && q!=-5
                       && (kh.freerank != 0
                                          || !kh.torsions.empty()))
                   || (i==-1 && q==-3 && (kh.freerank != 0
                                          || kh.torsions.size() != 1
                                          || kh.torsions[0] != 2))
                   || (i==-1 && q==-1 && (kh.freerank != 1
                                          || !kh.torsions.empty()))
                   || (i==-1 && q!=-3 && q!=-1 && (kh.freerank != 0
                                                   || !kh.torsions.empty()))
                   || (i==0 && std::abs(q)==1 && (kh.freerank != 1
                                                  || !kh.torsions.empty()))
                   || (i==0 && std::abs(q)!=1 && (kh.freerank != 0
                                                  || !kh.torsions.empty()))
                   || (i==1 && q==1 && (kh.freerank != 1 || !kh.torsions.empty()))
                   || (i==1 && q!=1 && (kh.freerank != 0 || !kh.torsions.empty()))
                   || (i==2 && q==3 && (kh.freerank != 0
                                        || kh.torsions.size() != 1
                                        || kh.torsions[0] != 2))
                   || (i==2 && q==5 && (kh.freerank != 1 || !kh.torsions.empty()))
                   || (i==2 && q!=3 && q!= 5 && (kh.freerank != 0
                                                 || !kh.torsions.empty())))
                {
                    ERR_MSG("Wrong cohomology group (i,q)="
                            << std::make_pair(i,q) << ":\n"
                            << kh.pretty());
                    return EXIT_FAILURE;
                }
            }
        }
    }

    // Figure-eight knot with V-smoothings.
    {
        auto fig8V = read_gauss_code(
            {1,-3,4,-1,2,-4,3,-2},
            {std::make_pair(3,false)}
            );

        if(!fig8V) {
            ERR_MSG("Failed to load the figure-eight knot.");
            return EXIT_FAILURE;
        }

        fig8V->makeSmoothV(3);
        fig8V->makeSmoothV(2);
        // Now, fig8V = (positive Hopf link) + (unknot)

        auto fig8Vcube = SmoothCube::fromDiagram(*fig8V);

        for(int q = -1; q <= 7; ++q) {
            auto enh_prop = get_enhancement_prop(*fig8V, fig8Vcube, q);
            auto ch = enh_prop
                ? khChain(*fig8V, fig8Vcube, *enh_prop)
                : std::nullopt;
            if((q%2==0 && ch) || (q%2!=0 && !ch)) {
                ERR_MSG("Illegal homology: q=" << q);
                return EXIT_FAILURE;
            }
            else if (q%2==0)
                continue;

            auto h = ch->compute();
            for(int i=-ch->maxdeg(); i <= -ch->mindeg(); ++i) {
                if(auto kh = h[-i-ch->mindeg()].compute();
                   !kh.torsions.empty()
                   || (i==0 && (q==-1 || q==3) && kh.freerank != 1)
                   || (i==0 && q==1 && kh.freerank != 2)
                   || (i==0 && q!=-1 && q!=1 && q!=3 && kh.freerank != 0)
                   || (i==1 && kh.freerank != 0)
                   || (i==2 && (q==3 || q==7) && kh.freerank != 1)
                   || (i==2 && q==5 && kh.freerank != 2)
                   || (i==2 && q!=3 && q!= 5 && q!= 7 && kh.freerank != 0))
                {
                    ERR_MSG("Wrong cohomology group (i,q)="
                            << std::make_pair(i,q) << ":\n"
                            << kh.pretty());
                    return EXIT_FAILURE;
                }
            }
        }
    }

    // Figure-eight knot with some crossings replaced with wide edges.
    {
        auto fig8W = read_gauss_code(
            {1,-3,4,-1,2,-4,3,-2},
            {std::make_pair(3,false)}
            );

        if(!fig8W) {
            ERR_MSG("Failed to load the figure-eight knot.");
            return EXIT_FAILURE;
        }

        fig8W->makeWide(1);
        fig8W->makeWide(0);
        // Now, fig8W = (unknot)

        auto fig8Wcube = SmoothCube::fromDiagram(*fig8W);

        for(int q = -1; q <= 1; ++q) {
            auto enh_prop = get_enhancement_prop(*fig8W, fig8Wcube, q);
            auto ch = enh_prop
                ? khChain(*fig8W, fig8Wcube, *enh_prop)
                : std::nullopt;
            if((q==0 && ch) || (q!=0 && !ch)) {
                ERR_MSG("Illegal homology: q=" << q);
                return EXIT_FAILURE;
            }
            else if (q==0)
                continue;

            auto h = ch->compute();
            auto kh = h[0-ch->mindeg()].compute();
            if(kh.freerank != 1 || !kh.torsions.empty())
            {
                ERR_MSG("Wrong cohomology group (i,q)="
                        << std::make_pair(0,q) << ":\n"
                        << kh.pretty());
                return EXIT_FAILURE;
            }
        }
    }

    // 6_2 knot
    {
        auto six_two = read_gauss_code(
            {-1,4,-3,1,-2,6,-5,3,-4,2,-6,5},
            {std::make_pair(5,false)});
        auto sixtwocube = SmoothCube::fromDiagram(*six_two);

        int qmin =
            - static_cast<int>(six_two->smoothing(0u).first)
            - static_cast<int>(six_two->nnegative())
            + six_two->writhe();
        int qmax =
            static_cast<int>(six_two->smoothing(~state_t(0u)).first)
            + static_cast<int>(six_two->npositive())
            + six_two->writhe();

        for(int q = qmin; q <= qmax; q+=2) {
            std::cout << "q-degree: " << q << std::endl;
            auto enh_prop = get_enhancement_prop(*six_two, sixtwocube, q);
            auto ch = enh_prop
                ? khChain(*six_two, sixtwocube, *enh_prop)
                : std::nullopt;
            if(!ch) {
                ERR_MSG("Cannot compute homology: q=" << q);
                return EXIT_FAILURE;
            }

            auto h = ch->compute();

            if(auto h0 = h[0].compute();
               h0.freerank != 0 || !h0.torsions.empty())
            {
                ERR_MSG(
                    "Incorrect 3rd cohomology groups: (deg,qdeg)=("
                    << -ch->mindeg() << "," << q
                    << ")");
                return EXIT_FAILURE;
            }
            if(auto h1 = h[1].compute();
               (q==3 && (h1.freerank != 1 || !h1.torsions.empty()))
               || (q==1
                   && (h1.freerank != 0
                       || h1.torsions.size() != 1 || h1.torsions[0] != 2))
               || (q!=3 && q != 1 && (h1.freerank != 0 || !h1.torsions.empty())))
            {
                ERR_MSG(
                    "Incorrect 2nd cohomology groups: (deg,qdeg)=("
                    << -ch->mindeg()-1 << "," << q
                    << ")");
                return EXIT_FAILURE;
            }
            if(auto h2 = h[2].compute();
               (q==-1 && (h2.freerank != 1 || !h2.torsions.empty()))
               || (q!=-1 && (h2.freerank != 0 || !h2.torsions.empty())))
            {
                ERR_MSG(
                    "Incorrect 1st cohomology groups: (deg,qdeg)=("
                    << -ch->mindeg()-2 << "," << q
                    << ")");
                return EXIT_FAILURE;
            }
            if(auto h3 = h[3].compute();
               (q==-1 && (h3.freerank != 2 || !h3.torsions.empty()))
               || (q==-3
                   && (h3.freerank != 1
                       || h3.torsions.size() != 1 || h3.torsions[0] != 2))
               || (q!=-1 && q!=-3 && (h3.freerank != 0 || !h3.torsions.empty())))
            {
                ERR_MSG(
                    "Incorrect 0th cohomology groups: (deg,qdeg)=("
                    << -ch->mindeg()-3 << "," << q
                    << ")");
                return EXIT_FAILURE;
            }
            if(auto h4 = h[4].compute();
               (q==-3 && (h4.freerank != 1 || !h4.torsions.empty()))
               || (q==-5
                   && (h4.freerank != 1
                       || h4.torsions.size() != 1 || h4.torsions[0] != 2))
               || (q!=-3 && q!=-5 && (h4.freerank != 0 || !h4.torsions.empty())))
            {
                ERR_MSG(
                    "Incorrect -1st cohomology groups: (deg,qdeg)=("
                    << -ch->mindeg()-4 << "," << q
                    << ")");
                return EXIT_FAILURE;
            }
            if(auto h5 = h[5].compute();
               (q==-5 && (h5.freerank != 1 || !h5.torsions.empty()))
               || (q==-7
                   && (h5.freerank != 1
                       || h5.torsions.size() != 1 || h5.torsions[0] != 2))
               || (q!=-5 && q!=-7 && (h5.freerank != 0 || !h5.torsions.empty())))
            {
                ERR_MSG(
                    "Incorrect -2nd cohomology groups: (deg,qdeg)=("
                    << -ch->mindeg()-5 << "," << q
                    << ")");
                return EXIT_FAILURE;
            }
            if(auto h6 = h[6].compute();
               (q==-7 && (h6.freerank != 1 || !h6.torsions.empty()))
               || (q==-9
                   && (h6.freerank != 1
                       || h6.torsions.size() != 1 || h6.torsions[0] != 2))
               || (q!=-7 && q!=-9 && (h6.freerank != 0 || !h6.torsions.empty())))
            {
                ERR_MSG(
                    "Incorrect -3rd cohomology groups: (deg,qdeg)=("
                    << -ch->mindeg()-6 << "," << q
                    << ")");
                return EXIT_FAILURE;
            }
            if(auto h7 = h[7].compute();
               (q==-11 && (h7.freerank != 1 || !h7.torsions.empty()))
               || (q!=-11 && (h7.freerank != 0 || !h7.torsions.empty())))
            {
                ERR_MSG(
                    "Incorrect -4th cohomology groups: (deg,qdeg)=("
                    << -ch->mindeg()-7 << "," << q
                    << ")");
                return EXIT_FAILURE;
            }
        }
    }

    return EXIT_SUCCESS;
}

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
            auto ch = khover::khChain(*trefoil, trefoilcube, q);

            if(q%2 == 0) {
                if (ch) {
                    ERR_MSG("Non-trivial complex in illegal quantum degrees.");
                    return -1;
                }
                continue;
            }

            if(!ch) {
                ERR_MSG("Failed to compute Khovanov homology of the trefiol");
                return -1;
            }

            auto h = ch->compute();

            if(h.size() != 5) {
                ERR_MSG("Wrong number of cohomology groups.");
                return -1;
            }

            if(auto h0 = h[0].compute();
               h0.freerank != 0 || !h0.torsions.empty())
            {
                ERR_MSG(
                    "Incorrect 1st cohomology groups: (deg,qdeg)=("
                    << -ch->mindeg() << "," << q
                    << ")");
                return -1;
            }
            if(auto h1 = h[1].compute();
               ((q==-1 || q==-3) && (h1.freerank != 1 || !h1.torsions.empty()))
               || (q!=-1 && q!=-3 && (h1.freerank != 0 || !h1.torsions.empty())))
            {
                ERR_MSG(
                    "Incorrect 0th cohomology groups: (deg,qdeg)=("
                    << -ch->mindeg()-1 << "," << q
                    << ")");
                return -1;
            }
            if(auto h2 = h[2].compute();
               h2.freerank != 0 || !h2.torsions.empty())
            {
                ERR_MSG(
                    "Incorrect -1st cohomology groups: (deg,qdeg)=("
                    << -ch->mindeg()-2 << "," << q
                    << ")");
                return -1;
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
                return -1;
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
            auto ch = khChain(*six_two, sixtwocube, q);
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

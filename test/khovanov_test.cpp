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

        for(int q = -9; q < 0 ; ++q) {
            auto ch = khover::khChain(*trefoil, q);

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
               h0.first != 0 || !h0.second.empty())
            {
                ERR_MSG(
                    "Incorrect 1st cohomology groups: (deg,qdeg)=("
                    << -ch->mindeg() << "," << q
                    << ")");
                return -1;
            }
            if(auto h1 = h[1].compute();
               ((q==-1 || q==-3) && (h1.first != 1 || !h1.second.empty()))
               || (q!=-1 && q!=-3 && (h1.first != 0 || !h1.second.empty())))
            {
                ERR_MSG(
                    "Incorrect 0th cohomology groups: (deg,qdeg)=("
                    << -ch->mindeg()-1 << "," << q
                    << ")");
                return -1;
            }
            if(auto h2 = h[2].compute();
               h2.first != 0 || !h2.second.empty())
            {
                ERR_MSG(
                    "Incorrect -1st cohomology groups: (deg,qdeg)=("
                    << -ch->mindeg()-2 << "," << q
                    << ")");
                return -1;
            }
            if(auto h3 = h[3].compute();
               (q==-5 && (h3.first != 1 || !h3.second.empty()))
               || (q==-7
                   && (h3.first != 0
                       || h3.second.size()!=1
                       || h3.second[0]!=2))
               || (q!=-5 && q!=-7 && (h3.first != 0 || !h3.second.empty())))
            {
                ERR_MSG(
                    "Incorrect -2nd cohomology groups: (deg,qdeg)=("
                    << -ch->mindeg()-3 << "," << q
                    << ")\n"
                    << h3
                    );
                return -1;
            }
            if(auto h4 = h[4].compute();
               (q==-9 && (h4.first != 1 || !h4.second.empty()))
               || (q!=-9 && (h4.first != 0 || !h4.second.empty())))
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
            auto ch = khChain(*six_two, q);
            if(!ch) {
                ERR_MSG("Cannot compute homology: q=" << q);
                return EXIT_FAILURE;
            }

            auto h = ch->compute();

            if(auto h0 = h[0].compute();
               h0.first != 0 || !h0.second.empty())
            {
                ERR_MSG(
                    "Incorrect 3rd cohomology groups: (deg,qdeg)=("
                    << -ch->mindeg() << "," << q
                    << ")");
                return EXIT_FAILURE;
            }
            if(auto h1 = h[1].compute();
               (q==3 && (h1.first != 1 || !h1.second.empty()))
               || (q==1
                   && (h1.first != 0
                       || h1.second.size() != 1 || h1.second[0] != 2))
               || (q!=3 && q != 1 && (h1.first != 0 || !h1.second.empty())))
            {
                ERR_MSG(
                    "Incorrect 2nd cohomology groups: (deg,qdeg)=("
                    << -ch->mindeg()-1 << "," << q
                    << ")");
                return EXIT_FAILURE;
            }
            if(auto h2 = h[2].compute();
               (q==-1 && (h2.first != 1 || !h2.second.empty()))
               || (q!=-1 && (h2.first != 0 || !h2.second.empty())))
            {
                ERR_MSG(
                    "Incorrect 1st cohomology groups: (deg,qdeg)=("
                    << -ch->mindeg()-2 << "," << q
                    << ")");
                return EXIT_FAILURE;
            }
            if(auto h3 = h[3].compute();
               (q==-1 && (h3.first != 2 || !h3.second.empty()))
               || (q==-3
                   && (h3.first != 1
                       || h3.second.size() != 1 || h3.second[0] != 2))
               || (q!=-1 && q!=-3 && (h3.first != 0 || !h3.second.empty())))
            {
                ERR_MSG(
                    "Incorrect 0th cohomology groups: (deg,qdeg)=("
                    << -ch->mindeg()-3 << "," << q
                    << ")");
                return EXIT_FAILURE;
            }
            if(auto h4 = h[4].compute();
               (q==-3 && (h4.first != 1 || !h4.second.empty()))
               || (q==-5
                   && (h4.first != 1
                       || h4.second.size() != 1 || h4.second[0] != 2))
               || (q!=-3 && q!=-5 && (h4.first != 0 || !h4.second.empty())))
            {
                ERR_MSG(
                    "Incorrect -1st cohomology groups: (deg,qdeg)=("
                    << -ch->mindeg()-4 << "," << q
                    << ")");
                return EXIT_FAILURE;
            }
            if(auto h5 = h[5].compute();
               (q==-5 && (h5.first != 1 || !h5.second.empty()))
               || (q==-7
                   && (h5.first != 1
                       || h5.second.size() != 1 || h5.second[0] != 2))
               || (q!=-5 && q!=-7 && (h5.first != 0 || !h5.second.empty())))
            {
                ERR_MSG(
                    "Incorrect -2nd cohomology groups: (deg,qdeg)=("
                    << -ch->mindeg()-5 << "," << q
                    << ")");
                return EXIT_FAILURE;
            }
            if(auto h6 = h[6].compute();
               (q==-7 && (h6.first != 1 || !h6.second.empty()))
               || (q==-9
                   && (h6.first != 1
                       || h6.second.size() != 1 || h6.second[0] != 2))
               || (q!=-7 && q!=-9 && (h6.first != 0 || !h6.second.empty())))
            {
                ERR_MSG(
                    "Incorrect -3rd cohomology groups: (deg,qdeg)=("
                    << -ch->mindeg()-6 << "," << q
                    << ")");
                return EXIT_FAILURE;
            }
            if(auto h7 = h[7].compute();
               (q==-11 && (h7.first != 1 || !h7.second.empty()))
               || (q!=-11 && (h7.first != 0 || !h7.second.empty())))
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

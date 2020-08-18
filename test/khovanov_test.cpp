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
                    << ")");
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
                return -1;
            }
        }
    }

    return 0;
}

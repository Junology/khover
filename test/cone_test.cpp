#include <iostream>
#include "chaincomplex.hpp"
#include "debug/debug.hpp"

using namespace khover;

int main(int argc, char* argv[])
{
    ChainIntegral D2(
        -1,
        ChainIntegral::matrix_t::Zero(2,2), // -1
        ChainIntegral::matrix_t::Zero(2,2), // 0
        (ChainIntegral::matrix_t(2,2) <<    // 1
         -2, -2,
         2, 2).finished(),
        ChainIntegral::matrix_t::Zero(2,2), // 2
        (ChainIntegral::matrix_t(2,2) <<    // 3
         -2, -2,
         2, 2).finished()
        );
    ChainIntegral D6(
        1,
        (ChainIntegral::matrix_t(6,6) <<    // 1
         -6, -6, -6, -6, -6, -6,
         30, 30, 30, 30, 30, 30,
         -60,-60,-60,-60,-60,-60,
         60, 60, 60, 60, 60, 60,
         -30,-30,-30,-30,-30,-30,
         6, 6, 6, 6, 6, 6).finished(),
        ChainIntegral::matrix_t::Zero(6,6), // 2
        (ChainIntegral::matrix_t(6,6) <<    // 3
         -6, -6, -6, -6, -6, -6,
         30, 30, 30, 30, 30, 30,
         -60,-60,-60,-60,-60,-60,
         60, 60, 60, 60, 60, 60,
         -30,-30,-30,-30,-30,-30,
         6, 6, 6, 6, 6, 6).finished(),
        ChainIntegral::matrix_t::Zero(6,6)  // 4
        );

    // Mapping cone of the zero morphism
    {
        ChainIntegral::Hom hom0(-1, {
                ChainIntegral::matrix_t::Zero(0,2),  // -1
                ChainIntegral::matrix_t::Zero(0,2),  // 0
                ChainIntegral::matrix_t::Zero(6,2),  // 1
                ChainIntegral::matrix_t::Zero(6,2),  // 2
                ChainIntegral::matrix_t::Zero(6,2),  // 3
                ChainIntegral::matrix_t::Zero(6,2),  // 4
                ChainIntegral::matrix_t::Zero(6,0)   // 5
                }
            );

        auto dsum = ChainIntegral::cone(hom0, D2, D6);

        if(!dsum) {
            ERR_MSG("Failed to construct the mapping cone of the zero map.");
            return EXIT_FAILURE;
        }

        int i = dsum->mindeg();
        for(auto h : dsum->compute()) {
            auto hi = h.compute();
            if (i == 0 &&
                (hi.freerank != 2 || !hi.torsions.empty())
                )
            {
                ERR_MSG("Wrong 0th homology: " << hi.pretty());
                return EXIT_FAILURE;
            }
            if (i == 1 &&
                (hi.freerank != 7
                 || hi.torsions.size() != 1 || hi.torsions[0] != 6)
                )
            {
                ERR_MSG("Wrong 1st homology: " << hi.pretty());
                return EXIT_FAILURE;
            }
            if (i == 2 &&
                (hi.freerank != 6
                 || hi.torsions.size() != 1 || hi.torsions[0] != 2)
                )
            {
                ERR_MSG("Wrong 2nd homology: " << hi.pretty());
                return EXIT_FAILURE;
            }
            if (i == 3 &&
                (hi.freerank != 6
                 || hi.torsions.size() != 1 || hi.torsions[0] != 6)
                )
            {
                ERR_MSG("Wrong 3rd homology: " << hi.pretty());
                return EXIT_FAILURE;
            }
            if (i == 4 &&
                (hi.freerank != 6
                 || hi.torsions.size() != 1 || hi.torsions[0] != 2)
                )
            {
                ERR_MSG("Wrong 4th homology: " << hi.pretty());
                return EXIT_FAILURE;
            }
            if (i == 5 &&
                (hi.freerank != 7 || !hi.torsions.empty())
                )
            {
                ERR_MSG("Wrong 5th homology: " << hi.pretty());
                return EXIT_FAILURE;
            }
            ++i;
        }
    }

    // Mapping cone of a non-zero map
    {
        ChainIntegral::Hom hom0(-1, {
                ChainIntegral::matrix_t::Zero(0,2),  // -1
                ChainIntegral::matrix_t::Zero(0,2),  // 0
                (ChainIntegral::matrix_t(6,2) <<
                 1, -2,
                 0, 15,
                 0,-30,
                 0, 30,
                 0,-15,
                 0,  3
                    ).finished(),  // 1
                (ChainIntegral::matrix_t(6,2) <<
                 1, -2,
                 0, 15,
                 0,-30,
                 0, 30,
                 0,-15,
                 0,  3
                    ).finished(),  // 2
                (ChainIntegral::matrix_t(6,2) <<
                 1, -2,
                 0, 15,
                 0,-30,
                 0, 30,
                 0,-15,
                 0,  3
                    ).finished(),  // 3
                (ChainIntegral::matrix_t(6,2) <<
                 1, -2,
                 0, 15,
                 0,-30,
                 0, 30,
                 0,-15,
                 0,  3
                    ).finished(),  // 4
                ChainIntegral::matrix_t::Zero(6,0)   // 5
                }
            );

        auto cn = ChainIntegral::cone(hom0, D2, D6);

        if(!cn) {
            ERR_MSG("Failed to construct the mapping cone of the zero map.");
            return EXIT_FAILURE;
        }

        int i = cn->mindeg();
        for(auto h : cn->compute()) {
            auto hi = h.compute();
            if (i == 0 &&
                (hi.freerank != 2 || !hi.torsions.empty())
                )
            {
                ERR_MSG("Wrong 0th homology: " << hi.pretty());
                return EXIT_FAILURE;
            }
            if (i == 1 &&
                (hi.freerank != 6
                 || hi.torsions.size() != 1 || hi.torsions[0] != 3)
                )
            {
                ERR_MSG("Wrong 1st homology: " << hi.pretty());
                return EXIT_FAILURE;
            }
            if (i == 2 &&
                (hi.freerank != 4
                 || hi.torsions.size() != 1 || hi.torsions[0] != 3)
                )
            {
                ERR_MSG("Wrong 2nd homology: " << hi.pretty());
                return EXIT_FAILURE;
            }
            if (i == 3 &&
                (hi.freerank != 4
                 || hi.torsions.size() != 1 || hi.torsions[0] != 3)
                )
            {
                ERR_MSG("Wrong 3rd homology: " << hi.pretty());
                return EXIT_FAILURE;
            }
            if (i == 4 &&
                (hi.freerank != 4
                 || hi.torsions.size() != 1 || hi.torsions[0] != 3)
                )
            {
                ERR_MSG("Wrong 4th homology: " << hi.pretty());
                return EXIT_FAILURE;
            }
            if (i == 5 &&
                (hi.freerank != 6 || !hi.torsions.empty())
                )
            {
                ERR_MSG("Wrong 5th homology: " << hi.pretty());
                return EXIT_FAILURE;
            }
            ++i;
        }
    }
    return EXIT_SUCCESS;
}

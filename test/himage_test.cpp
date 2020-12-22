#include <iostream>
#include "chaincomplex.hpp"
#include "debug/debug.hpp"

using namespace khover;

int main(int argc, char* argv[])
{
    /*!
      Consider the inclusion
        i:\partial\Delta^2 \hookrightarrow \Delta^2.
      In this test, we check if ChainIntegral::compute works correctly for the maps \Delta^2\to Cone(i) and Cone(i)\to \partial\Delta^2.
     */
    ChainIntegral BD2(
        0,
        (ChainIntegral::matrix_t(3,3) <<    // 0
         0, -1, -1,
         -1, 0, 1,
         1, 1, 0).finished()
        );
    ChainIntegral D2(
        0,
        (ChainIntegral::matrix_t(3,3) <<    // 0
         0, -1, -1,
         -1, 0,  1,
         1,  1,  0).finished(),
        (ChainIntegral::matrix_t(3,1) <<    // 1
         1, -1, 1).finished()
        );
    ChainIntegral C(
        0,
        (ChainIntegral::matrix_t(3,6) <<    // 0
         0, -1, -1, 1, 0, 0,
         -1, 0,  1, 0, 1, 0,
         1,  1,  0, 0, 0, 1).finished(),
        (ChainIntegral::matrix_t(6,4) <<    // 1
         1,  1, 0, 0,
         -1, 0, 1, 0,
         1,  0, 0, 1,
         0,  0, 1, 1,
         0,  1, 0, -1,
         0, -1, -1, 0).finished()
        );

    // Check first that the complexes give correct homology groups.
    {
        std::cout << "\e[34;1m---\nHomology groups of \\partial\\Delta^2\n---\e[m" << std::endl;

        auto hs = BD2.compute();

        if (hs.size() != 2) {
            ERR_MSG("Wrong number of homology groups.");
            return EXIT_FAILURE;
        }

        if (auto h0 = hs[0].compute();
            h0.freerank != 1 || !h0.torsions.empty()) {
            ERR_MSG("Wrong 0th homology.");
            return EXIT_FAILURE;
        }

        if (auto h1 = hs[1].compute();
            h1.freerank != 1 || !h1.torsions.empty()) {
            ERR_MSG("Wrong 1st homology.");
            return EXIT_FAILURE;
        }

        std::cout << "Passed" << std::endl;
    }

    {
        std::cout << "\e[34;1m---\nHomology groups of \\Delta^2\n---\e[m" << std::endl;

        auto hs = D2.compute();

        if (hs.size() != 3) {
            ERR_MSG("Wrong number of homology groups.");
            return EXIT_FAILURE;
        }

        if (auto h0 = hs[0].compute();
            h0.freerank != 1 || !h0.torsions.empty()) {
            ERR_MSG("Wrong 0th homology.");
            return EXIT_FAILURE;
        }

        if (auto h1 = hs[1].compute();
            h1.freerank != 0 || !h1.torsions.empty()) {
            ERR_MSG("Wrong 1st homology.");
            return EXIT_FAILURE;
        }

        if (auto h2 = hs[2].compute();
            h2.freerank != 0 || !h2.torsions.empty()) {
            ERR_MSG("Wrong 2nd homology.");
            return EXIT_FAILURE;
        }

        std::cout << "Passed" << std::endl;
    }

    {
        std::cout << "\e[34;1m---\nHomology groups of \\Delta^2//\\partial\\Delta^2\n---\e[m" << std::endl;

        auto hs = C.compute();

        if (hs.size() != 3) {
            ERR_MSG("Wrong number of homology groups.");
            return EXIT_FAILURE;
        }

        if (auto h0 = hs[0].compute();
            h0.freerank != 0 || !h0.torsions.empty()) {
            ERR_MSG("Wrong 0th homology.");
            return EXIT_FAILURE;
        }

        if (auto h1 = hs[1].compute();
            h1.freerank != 0 || !h1.torsions.empty()) {
            ERR_MSG("Wrong 1st homology.");
            return EXIT_FAILURE;
        }

        if (auto h2 = hs[2].compute();
            h2.freerank != 1 || !h2.torsions.empty()) {
            ERR_MSG("Wrong 2nd homology.");
            return EXIT_FAILURE;
        }

        std::cout << "Passed" << std::endl;
    }

    // hom induced by \\Delta^2 \to Cone,
    // which has the trivial image on the homology
    {
        std::cout << "\e[34;1m---\n\\Delta^2\\to Cone\n---\e[m" << std::endl;
        ChainIntegral::Hom inc(0, {
                ChainIntegral::matrix_t::Identity(3,3),  // 0
                ChainIntegral::matrix_t::Identity(6,3),  // 1
                ChainIntegral::matrix_t::Identity(4,1),  // 2
                }
            );

        auto res_dom = D2.compute({},inc);
        auto res_cod = C.compute(inc,{});

        for(int i = 0; i < 3; ++i) {
            auto im = res_cod[i].image(inc.morphism(i));
            if (!im) {
                ERR_MSG("Failed to compute the image: " << i);
                return EXIT_FAILURE;
            }
            auto grp = im->compute();
            if (grp.freerank != 0 || !grp.torsions.empty()) {
                ERR_MSG("Wrong image at " << i);
                return EXIT_FAILURE;
            }
        }
        std::cout << "Passed" << std::endl;
    }

    // hom induced by Cone \to \\partial\\Delta^2[-1]
    {
        std::cout << "\e[34;1m---\nCone\\to\\partial\\Delta^2[-1]\n---\e[m" << std::endl;
        ChainIntegral::Hom proj(0, {
                ChainIntegral::matrix_t::Zero(0,3),  // 0
                (ChainIntegral::matrix_t(3,6) <<     // 1
                 ChainIntegral::matrix_t::Zero(3,3),
                 ChainIntegral::matrix_t::Identity(3,3)).finished(),
                (ChainIntegral::matrix_t(3,4) <<     // 2
                 ChainIntegral::matrix_t::Zero(3,1),
                 ChainIntegral::matrix_t::Identity(3,3)).finished()
            });
        auto SBD2 = BD2 << 1;

        auto res_dom = C.compute({},proj);
        auto res_cod = SBD2.compute(proj,{});

        auto im1 = res_cod[0].image(proj.morphism(1));
        if (!im1) {
            ERR_MSG("Failed to compute the image: " << 1);
            return EXIT_FAILURE;
        }
        auto grp1 = im1->compute();
        if (grp1.freerank != 0 || !grp1.torsions.empty()) {
            ERR_MSG("Wrong image at " << 1);
            return EXIT_FAILURE;
        }

        auto res = res_cod[1].image(proj.morphism(2), std::true_type{});
        if (!res) {
            ERR_MSG("Failed to compute the image: " << 2);
            return EXIT_FAILURE;
        }
        auto [im2,mat] = *res;

        auto grp2 = im2.compute();
        if (grp2.freerank != 1 || !grp2.torsions.empty()) {
            ERR_MSG("Wrong image at " << 2);
            return EXIT_FAILURE;
        }

        if(std::abs(mat.cast<double>().determinant()) != 1) {
            ERR_MSG("Matrix is not unimodular");
            return EXIT_FAILURE;
        }

        std::cout << "Passed" << std::endl;
    }
    return EXIT_SUCCESS;
}

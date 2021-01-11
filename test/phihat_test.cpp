#include <iostream>
#include <iomanip>

#include "linkdiagram.hpp"
#include "khovanov.hpp"

#include "debug/debug.hpp"
using namespace khover;

template <class T, std::size_t n>
constexpr std::size_t lengthof(T const(&)[n]) {
    return n;
}

int test_phihat(LinkDiagram const& diagram, LinkDiagram::crossing_t c)
{
    auto diag_neg = diagram, diag_pos = diagram;
    diag_neg.makeNegative(c);
    diag_pos.makePositive(c);

    // Make cubes
    auto cube_neg = SmoothCube::fromDiagram(diag_neg);
    auto cube_pos = SmoothCube::fromDiagram(diag_pos);

    // Compute the bounds of q-gradings
    int qmin =
        std::min(
            - static_cast<int>(diag_pos.smoothing(0u).first)
            - static_cast<int>(diag_pos.nnegative())
            + diag_pos.writhe(),
            - static_cast<int>(diag_neg.smoothing(0u).first)
            - static_cast<int>(diag_neg.nnegative())
            + diag_neg.writhe()
            );
    int qmax =
        std::max(
            static_cast<int>(diag_pos.smoothing(~state_t(0u)).first)
            + static_cast<int>(diag_pos.npositive())
            + diag_pos.writhe(),
            static_cast<int>(diag_neg.smoothing(~state_t(0u)).first)
            + static_cast<int>(diag_neg.npositive())
            + diag_neg.writhe()
            );
    for(int q = qmin; q <= qmax; q+=2) {
        // Compute enhancement properties
        auto enhprop_neg = get_enhancement_prop(diag_neg, cube_neg, q);
        auto enhprop_pos = get_enhancement_prop(diag_pos, cube_pos, q);

        if (!enhprop_neg || !enhprop_pos) {
            ERR_MSG("Failed to compute enhancement properties.");
            return EXIT_FAILURE;
        }

        // Compute Khovanov complexes
        auto ch_neg = khChain(diag_neg, cube_neg, *enhprop_neg);
        auto ch_pos = khChain(diag_pos, cube_pos, *enhprop_pos);

        if (!ch_neg || !ch_pos) {
            ERR_MSG("Failed to construct chain complexes");
            return EXIT_FAILURE;
        }

        // Compute PhiHat
        auto phihat = khover::crossPhiHat(
                diag_pos, c,
                cube_neg, *enhprop_neg,
                cube_pos, *enhprop_pos);

        if (!phihat) {
            ERR_MSG("Failed to compute PhiHat.");
            return EXIT_FAILURE;
        }

        // Check if Phihat is a chain map.
        for(int i = std::max(ch_neg->maxdeg(), ch_pos->maxdeg()); i > std::min(ch_neg->mindeg(), ch_pos->mindeg()); --i) {
            auto comm = ch_pos->getDiff(i) * phihat->morphism(i+1) - phihat->morphism(i) * ch_neg->getDiff(i);
            if(!comm.isZero()) {
                ERR_MSG("Phihat does not commute with differentials at (q,i)=(" << q << "," << i << ")");
                DBG_MSG("Diff of C(+):\n" << ch_pos->getDiff(i));
                DBG_MSG("Diff of C(-):\n" << ch_neg->getDiff(i));
                DBG_MSG("Phihat@i+1\n" << phihat->morphism(i+1));
                DBG_MSG("Phihat@i\n" << phihat->morphism(i));
                return EXIT_FAILURE;
            }
        }
    }

    return EXIT_SUCCESS;
}

int main (int argc, char* argv[])
{
    std::vector<std::pair<std::size_t,bool>> sgn = {std::make_pair(0,false)};
    // Diagrams of prime knots with <= 7 crossings
    // (from https://knotinfo.math.indiana.edu/)
    std::optional<LinkDiagram> diagrams[] = {
        read_gauss_code({1, -2, 3, -1, 2, -3}, sgn),
        read_gauss_code({-1, 2, -3, 1, -4, 3, -2, 4}, sgn),
        read_gauss_code({-1, 2, -3, 4, -5, 1, -2, 3, -4, 5}, sgn),
        read_gauss_code({1, -2, 3, -1, 4, -5, 2, -3, 5, -4}, sgn),
        read_gauss_code({1, -2, 3, -4, 2, -1, 5, -6, 4, -3, 6, -5}, sgn),
        read_gauss_code({1, -2, 3, -4, 5, -6, 2, -1, 6, -3, 4, -5}, sgn),
        read_gauss_code({-1, 2, -3, 1, -4, 5, -2, 3, -6, 4, -5, 6}, sgn),
        read_gauss_code({1, -2, 3, -4, 5, -6, 7, -1, 2, -3, 4, -5, 6, -7}, sgn),
        read_gauss_code({-1, 2, -3, 4, -5, 6, -7, 1, -2, 7, -6, 5, -4, 3}, sgn),
        read_gauss_code({1, -2, 3, -4, 5, -6, 7, -1, 2, -3, 4, -7, 6, -5}, sgn),
        read_gauss_code({-1, 2, -3, 4, -5, 6, -7, 3, -2, 1, -4, 7, -6, 5}, sgn),
        read_gauss_code({-1, 2, -3, 1, -4, 5, -6, 7, -2, 3, -7, 4, -5, 6}, sgn),
        read_gauss_code({1, -2, 3, -4, 5, -6, 7, -3, 2, -7, 6, -1, 4, -5}, sgn),
        read_gauss_code({1, -2, 3, -4, 5, -6, 4, -7, 2, -1, 7, -3, 6, -5}, sgn)
    };


    // Check if all the diagrams are loaded correctly.
    for(std::size_t i = 0; i < lengthof(diagrams); ++i) {
        if (!diagrams[i]) {
            ERR_MSG("Failed to load the diagram@i=" << i);
            return EXIT_FAILURE;
        }
    }

    for(std::size_t i = 0; i < lengthof(diagrams); ++i) {
        for(std::size_t c = 0; c < diagrams[i]->ncrosses(); ++c) {
            std::cout << "\e[34;1m---\nCheck: (" << i << "," << c << ")\n---\e[m" << std::endl;

            std::cout << "direct" << std::endl;
            if (test_phihat(*(diagrams[i]), c) == EXIT_FAILURE)
                return EXIT_FAILURE;

            std::cout << "mirror" << std::endl;
            diagrams[i]->mirroring();
            if (test_phihat(*(diagrams[i]), c) == EXIT_FAILURE)
                return EXIT_FAILURE;

            DBG_MSG("Passed");
        }
    }

    return EXIT_SUCCESS;
}

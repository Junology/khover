#include <iostream>
#include <algorithm>

#include "linkdiagram.hpp"

#include "debug/debug.hpp"

using namespace khover;

int main(int argc, char* argv[])
{
    auto six_two = read_gauss_code(
        {1,-2, 5,-6, 3,-1, 2,-4, 6,-5, 4,-3},
        {std::make_pair(1,false)});

    if(!six_two) {
        ERR_MSG("Failed to load a knot.");
        return EXIT_FAILURE;
    }

    if (six_two->nnegative() != 4) {
        ERR_MSG("Wrong number of negative crossings: "
                << six_two->nnegative() << " should be " << 4);
        return EXIT_FAILURE;
    }

    if (six_two->npositive() != 2) {
        ERR_MSG("Wrong number of positive crossings: "
                << six_two->npositive() << " should be " << 2);
        return EXIT_FAILURE;
    }

    for(unsigned long s = 0; s < 0b1000000 ; ++s) {
        if(six_two->cohDegree(s) != static_cast<int>(std::bitset<6>(s).count()) - 4)
        {
            ERR_MSG("Wrong cohomological degree for " << std::bitset<6>(s));
            return EXIT_FAILURE;
        }
    }

    // Smoothing test: all 0
    {
        auto [ncomp, comptbl] = six_two->smoothing(0b000000u);
        std::vector<unsigned char> shouldbe = {0, 1, 2, 3, 4, 0, 1, 3, 2, 1, 3, 4};

        if (ncomp != 5
            || !std::equal(comptbl.begin(), comptbl.end(), shouldbe.begin()))
        {
            ERR_MSG("Wrong components:" << comptbl);
            return EXIT_FAILURE;
        }
    }

    // Smoothing test: 0b000001
    {
        auto [ncomp, comptbl] = six_two->smoothing(0b000001u);
        std::vector<unsigned char> shouldbe = {0, 1, 2, 3, 0, 0, 1, 3, 2, 1, 3, 0};

        if (ncomp != 4
            || !std::equal(comptbl.begin(), comptbl.end(), shouldbe.begin()))
        {
            ERR_MSG("Wrong components:" << comptbl);
            return EXIT_FAILURE;
        }
    }

    // Smoothing test: 0b000010
    {
        auto [ncomp, comptbl] = six_two->smoothing(0b000010u);
        std::vector<unsigned char> shouldbe = {0, 0, 1, 2, 3, 0, 0, 2, 1, 0, 2, 3};

        if (ncomp != 4
            || !std::equal(comptbl.begin(), comptbl.end(), shouldbe.begin()))
        {
            ERR_MSG(
                "Wrong components:" << comptbl << "\n"
                << "Shouldbe:" << shouldbe);
            return EXIT_FAILURE;
        }
    }

    // Smoothing test: 0b000011
    {
        auto [ncomp,comptbl] = six_two->smoothing(0b000011u);
        std::vector<unsigned char> shouldbe = {0, 0, 1, 2, 0, 0, 0, 2, 1, 0, 2, 0};

        if (ncomp != 3
            || !std::equal(comptbl.begin(), comptbl.end(), shouldbe.begin()))
        {
            ERR_MSG("Wrong components:" << comptbl);
            return EXIT_FAILURE;
        }
    }

    // Smoothing after V-smoothing
    {
        LinkDiagram aux = *six_two;
        aux.makeSmoothV(5);
        aux.makeSmoothV(4);
        if(aux.ncrosses() != 4 || aux.narcs() != 9) {
            ERR_MSG("Wrong number of crossings/components:\n"
                    << "(ncrx,narc)="
                    << std::make_pair(aux.ncrosses(), aux.narcs()));
            return EXIT_FAILURE;
        }
        auto [ncomp, comptbl] = aux.smoothing(0b000000u);
        std::vector<unsigned char> shouldbe = {0, 1, 2, 3, 4, 0, 1, 3, 4};

        if (ncomp != 5
            || !std::equal(comptbl.begin(), comptbl.end(), shouldbe.begin())) {
            ERR_MSG("Wrong components:" << comptbl);
            return EXIT_FAILURE;
        }
    }

    // Smoothing with wide edges
    {
        LinkDiagram aux = *six_two;

        aux.makeWide(2);
        aux.makeWide(1);
        auto [ncomp, comptbl] = aux.smoothing(0b0000u);
        std::vector<unsigned char> shouldbe = {0, 1, 2, 3, 4, 0, 1, 3, 2, 1, 3, 4};

        if (ncomp != 5
            || !std::equal(comptbl.begin(), comptbl.end(), shouldbe.begin())) {
            ERR_MSG("Wrong components:" << comptbl);
            return EXIT_FAILURE;
        }
    }

    // Adjacent state test:
    {
        state_t st0 = 0b010100;
        state_t st1 = 0b010101;
        state_t st2 = 0b011101;
        state_t st3 = 0b111101;

        if(six_two->stateCoeff(st0, st1) != 1
           || six_two->stateCoeff(st1, st0) != 0)
        {
            std::cerr << "Wrong state adjacency: " << std::endl
                      << st0 << std::endl << st1 << std::endl
                      << "Coeffs = ("
                      << six_two->stateCoeff(st0, st1) << ","
                      << six_two->stateCoeff(st1, st0) << ")"
                      << std::endl;
            return -1;
        }

        if(six_two->stateCoeff(st1, st2) != -1
           || six_two->stateCoeff(st2, st1) != 0)
        {
            std::cerr << "Wrong state adjacency: " << std::endl
                      << st1 << std::endl << st2 << std::endl
                      << "Coeffs = ("
                      << six_two->stateCoeff(st1, st2) << ","
                      << six_two->stateCoeff(st2, st1) << ")"
                      << std::endl;
            return -1;
        }

        if(six_two->stateCoeff(st2, st3) != 1
           || six_two->stateCoeff(st3, st2) != 0)
        {
            std::cerr << "Wrong state adjacency: " << std::endl
                      << st2 << std::endl << st3 << std::endl
                      << "Coeffs = ("
                      << six_two->stateCoeff(st2, st3) << ","
                      << six_two->stateCoeff(st3, st2) << ")"
                      << std::endl;
            return -1;
        }
    }

    // Skew-commutativity of cube.
    for(std::size_t k = 0; k < 32; ++k) {
        state_t st00 = k;

        for(std::size_t i = 0; i < 5; ++i) {
            if(st00.test(i))
                continue;

            for(std::size_t j = i+1; j < 6; ++j) {
                if(st00.test(j))
                    continue;

                state_t st01 = st00, st10 = st00, st11 = st00;
                st01.set(i);
                st10.set(j);
                st11.set(i);
                st11.set(j);

                if(six_two->stateCoeff(st00, st01)*six_two->stateCoeff(st01, st11)
                   != -six_two->stateCoeff(st00, st10)*six_two->stateCoeff(st10, st11))
                {
                    std::cerr << "Cube is not skew-commutative:" << std::endl
                              << st00 << std::endl
                              << "(i,j)=(" << i << "," << j << ")" << std::endl;
                    return -1;
                }
            }
        }
    }
    {
        auto trefoil = read_gauss_code(
            {1, -3, 2, -1, 3, -2}, {std::make_pair(1,false)});
        if(auto comps = trefoil->smoothing(state_t{0b011});
           comps.first != 1
           || comps.second[0] != 0
           || comps.second[1] != 0
           || comps.second[2] != 0
           || comps.second[3] != 0
           || comps.second[4] != 0
           || comps.second[5] != 0
            )
        {
            std::cerr << "Wrong components:" << std::endl;
            for(auto c : comps.second) {
                std::cout << static_cast<int>(c) << "; ";
            }
            std::cout << std::endl;
            return -1;
        }
    }

    // Crux smoothing test
    {
        for(std::size_t st = 0; st < cipow(2,5); ++st) {
            auto stbit = std::bitset<max_crosses>(st);
            auto crux_result = six_two->cruxTwists(st, 4);

            if (stbit.test(4) && stbit.test(3))
            {
                if (!crux_result) {
                    ERR_MSG("Wrong detection of non-crux state:\n" << stbit);
                    return EXIT_FAILURE;
                }

                if(auto [stbit,is_twisted] = *crux_result;
                   is_twisted.to_ulong() != 0b001010000100)
                {
                    ERR_MSG("Wrong twisted arcs:\n" << stbit << "\n" << is_twisted);
                    return EXIT_FAILURE;
                }
            }
            else if (stbit.test(4) && stbit.test(0) && stbit.test(1) && stbit.test(2))
            {
                if (!crux_result) {
                    ERR_MSG("Wrong detection of non-crux state:\n" << stbit);
                    return EXIT_FAILURE;
                }

                if(auto [stbit,is_twisted] = *crux_result;
                   st == 0b100111 && is_twisted.to_ulong() != 0b011011010101ul)
                {
                    ERR_MSG("Wrong twisted arcs:\n" << stbit << "\n" << is_twisted);
                    return EXIT_FAILURE;
                }
            }
            else {
                if (crux_result) {
                    ERR_MSG("Failed to detect non-crux states:" << stbit);
                    return EXIT_FAILURE;
                }
            }
        }
    }

    auto six_two2 = read_gauss_code(
        {-1,4,-3,1,-2,6,-5,3,-4,2,-6,5},
        {std::make_pair(5,false)});

    if(!six_two2) {
        ERR_MSG("Failed to load the knot 6_2.");
        return -1;
    }

    if (six_two2->getSign(0) >= 0
        || six_two2->getSign(1) >= 0
        || six_two2->getSign(2) <= 0
        || six_two2->getSign(3) <= 0
        || six_two2->getSign(4) >= 0
        || six_two2->getSign(5) >= 0
        )
    {
        ERR_MSG("Wrong signs on crossings.");
        return -1;
    }

    for(unsigned long k = 0; k < cipow(2,6); ++k) {
        if (auto comps = six_two2->smoothing(state_t{k});
            (k == 0b000000ul && comps.first != 5)
            || (k == 0b000001ul && comps.first != 4)
            || (k == 0b000010ul && comps.first != 4)
            || (k == 0b000011ul && comps.first != 3)
            || (k == 0b000100ul && comps.first != 4)
            || (k == 0b000101ul && comps.first != 3)
            || (k == 0b000110ul && comps.first != 3)
            || (k == 0b000111ul && comps.first != 2)
            || (k == 0b001000ul && comps.first != 4)
            || (k == 0b001001ul && comps.first != 3)
            || (k == 0b001010ul && comps.first != 3)
            || (k == 0b001011ul && comps.first != 2)
            || (k == 0b001100ul && comps.first != 3)
            || (k == 0b001101ul && comps.first != 4)
            || (k == 0b001110ul && comps.first != 2)
            || (k == 0b001111ul && comps.first != 3)
            )
        {
            ERR_MSG("Wrong smoothings");
            return -1;
        }
    }

    return 0;
}

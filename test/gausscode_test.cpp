#include <iostream>

#include "linkdiagram.hpp"
#include "debug/debug.hpp"

using namespace khover;

int main(int argc, char* argv[])
{
    // Trefoil.
    {
        std::cout << "trefoil:" << std::endl;
        auto trefoil = read_gauss_code(
            {1, -3, 2, -1, 3, -2}, {std::make_pair(1,false)});

        if (!trefoil
            || trefoil->ncrosses() != 3
            || trefoil->narcs() != 6
            || trefoil->getSign(0) >= 0
            || trefoil->crosses()[0].adj_arc[0] != 5
            || trefoil->crosses()[0].adj_arc[1] != 0
            || trefoil->crosses()[0].adj_arc[2] != 2
            || trefoil->crosses()[0].adj_arc[3] != 3
            || trefoil->getSign(1) >= 0
            || trefoil->crosses()[1].adj_arc[0] != 1
            || trefoil->crosses()[1].adj_arc[1] != 2
            || trefoil->crosses()[1].adj_arc[2] != 4
            || trefoil->crosses()[1].adj_arc[3] != 5
            || trefoil->getSign(2) >= 0
            || trefoil->crosses()[2].adj_arc[0] != 3
            || trefoil->crosses()[2].adj_arc[1] != 4
            || trefoil->crosses()[2].adj_arc[2] != 0
            || trefoil->crosses()[2].adj_arc[3] != 1
            )
        {
            std::cerr << "Failed to load trefoil." << std::endl;
            return -1;
        }

        std::cout << "Pass" << std::endl;
    }

    // 6_2
    {
        std::cout << "6_2:" << std::endl;
        auto six_two = read_gauss_code(
            {1,-2, 5,-6, 3,-1, 2,-4, 6,-5, 4,-3},
            {std::make_pair(1,false)});

        if(!six_two
           || six_two->narcs() != 12
           || six_two->ncrosses() != 6
           || six_two->getSign(0) >= 0
           || six_two->crosses()[0].adj_arc[0] != 11
           || six_two->crosses()[0].adj_arc[1] != 0
           || six_two->crosses()[0].adj_arc[2] != 4
           || six_two->crosses()[0].adj_arc[3] != 5
           || six_two->getSign(1) >= 0
           || six_two->crosses()[1].adj_arc[0] != 5
           || six_two->crosses()[1].adj_arc[1] != 6
           || six_two->crosses()[1].adj_arc[2] != 0
           || six_two->crosses()[1].adj_arc[3] != 1
           || six_two->getSign(2) >= 0
           || six_two->crosses()[2].adj_arc[0] != 3
           || six_two->crosses()[2].adj_arc[1] != 4
           || six_two->crosses()[2].adj_arc[2] != 10
           || six_two->crosses()[2].adj_arc[3] != 11
           || six_two->getSign(3) >= 0
           || six_two->crosses()[3].adj_arc[0] != 9
           || six_two->crosses()[3].adj_arc[1] != 10
           || six_two->crosses()[3].adj_arc[2] != 6
           || six_two->crosses()[3].adj_arc[3] != 7
           || six_two->getSign(4) <= 0
           || six_two->crosses()[4].adj_arc[0] != 1
           || six_two->crosses()[4].adj_arc[1] != 2
           || six_two->crosses()[4].adj_arc[2] != 8
           || six_two->crosses()[4].adj_arc[3] != 9
           || six_two->getSign(5) <= 0
           || six_two->crosses()[5].adj_arc[0] != 7
           || six_two->crosses()[5].adj_arc[1] != 8
           || six_two->crosses()[5].adj_arc[2] != 2
           || six_two->crosses()[5].adj_arc[3] != 3
            ) {
            std::cerr << "Failed to load 6_2." << std::endl;
            return -1;
        }

        std::cout << "Pass" << std::endl;
    }

    // 6_2 (mirror)
    {
        std::cout << "6_2 (mirror):" << std::endl;
        auto six_two = read_gauss_code(
            {1,-2, 5,-6, 3,-1, 2,-4, 6,-5, 4,-3},
            {std::make_pair(1,true)});

        if(!six_two
           || six_two->narcs() != 12
           || six_two->ncrosses() != 6
           || six_two->getSign(0) <= 0
           || six_two->crosses()[0].adj_arc[0] != 11
           || six_two->crosses()[0].adj_arc[1] != 0
           || six_two->crosses()[0].adj_arc[2] != 4
           || six_two->crosses()[0].adj_arc[3] != 5
           || six_two->getSign(1) <= 0
           || six_two->crosses()[1].adj_arc[0] != 5
           || six_two->crosses()[1].adj_arc[1] != 6
           || six_two->crosses()[1].adj_arc[2] != 0
           || six_two->crosses()[1].adj_arc[3] != 1
           || six_two->getSign(2) <= 0
           || six_two->crosses()[2].adj_arc[0] != 3
           || six_two->crosses()[2].adj_arc[1] != 4
           || six_two->crosses()[2].adj_arc[2] != 10
           || six_two->crosses()[2].adj_arc[3] != 11
           || six_two->getSign(3) <= 0
           || six_two->crosses()[3].adj_arc[0] != 9
           || six_two->crosses()[3].adj_arc[1] != 10
           || six_two->crosses()[3].adj_arc[2] != 6
           || six_two->crosses()[3].adj_arc[3] != 7
           || six_two->getSign(4) >= 0
           || six_two->crosses()[4].adj_arc[0] != 1
           || six_two->crosses()[4].adj_arc[1] != 2
           || six_two->crosses()[4].adj_arc[2] != 8
           || six_two->crosses()[4].adj_arc[3] != 9
           || six_two->getSign(5) >= 0
           || six_two->crosses()[5].adj_arc[0] != 7
           || six_two->crosses()[5].adj_arc[1] != 8
           || six_two->crosses()[5].adj_arc[2] != 2
           || six_two->crosses()[5].adj_arc[3] != 3
            ) {
            std::cerr << "Failed to load 6_2 (mirror)." << std::endl;
            return -1;
        }

        std::cout << "Pass" << std::endl;
    }

    // 6_2 (mirror)
    {
        std::cout << "6_2 (mirror):" << std::endl;
        auto six_two = read_gauss_code(
            {1,-2, 5,-6, 3,-1, 2,-4, 6,-5, 4,-3},
            {std::make_pair(5,false)});

        if(!six_two
           || six_two->narcs() != 12
           || six_two->ncrosses() != 6
           || six_two->getSign(0) <= 0
           || six_two->crosses()[0].adj_arc[0] != 11
           || six_two->crosses()[0].adj_arc[1] != 0
           || six_two->crosses()[0].adj_arc[2] != 4
           || six_two->crosses()[0].adj_arc[3] != 5
           || six_two->getSign(1) <= 0
           || six_two->crosses()[1].adj_arc[0] != 5
           || six_two->crosses()[1].adj_arc[1] != 6
           || six_two->crosses()[1].adj_arc[2] != 0
           || six_two->crosses()[1].adj_arc[3] != 1
           || six_two->getSign(2) <= 0
           || six_two->crosses()[2].adj_arc[0] != 3
           || six_two->crosses()[2].adj_arc[1] != 4
           || six_two->crosses()[2].adj_arc[2] != 10
           || six_two->crosses()[2].adj_arc[3] != 11
           || six_two->getSign(3) <= 0
           || six_two->crosses()[3].adj_arc[0] != 9
           || six_two->crosses()[3].adj_arc[1] != 10
           || six_two->crosses()[3].adj_arc[2] != 6
           || six_two->crosses()[3].adj_arc[3] != 7
           || six_two->getSign(4) >= 0
           || six_two->crosses()[4].adj_arc[0] != 1
           || six_two->crosses()[4].adj_arc[1] != 2
           || six_two->crosses()[4].adj_arc[2] != 8
           || six_two->crosses()[4].adj_arc[3] != 9
           || six_two->getSign(5) >= 0
           || six_two->crosses()[5].adj_arc[0] != 7
           || six_two->crosses()[5].adj_arc[1] != 8
           || six_two->crosses()[5].adj_arc[2] != 2
           || six_two->crosses()[5].adj_arc[3] != 3
            ) {
            std::cerr << "Failed to load 6_2 (mirror)." << std::endl;
            return -1;
        }

        std::cout << "Pass" << std::endl;
    }

    // Borromean ring
    {
        std::cout << "Borromean ring (L6a4):" << std::endl;
        auto borromean = read_gauss_code(
            {1,-2, 3,-4, 0,-1, 5,-3, 6, 0, 2,-5, 4,-6},
            {std::make_pair(1,false)});

        if(!borromean
           || borromean->narcs() != 12
           || borromean->ncrosses() != 6
           || borromean->getSign(0) >= 0
           || borromean->crosses()[0].adj_arc[0] != 3
           || borromean->crosses()[0].adj_arc[1] != 0
           || borromean->crosses()[0].adj_arc[2] != 7
           || borromean->crosses()[0].adj_arc[3] != 4
           || borromean->getSign(1) <= 0
           || borromean->crosses()[1].adj_arc[0] != 11
           || borromean->crosses()[1].adj_arc[1] != 8
           || borromean->crosses()[1].adj_arc[2] != 0
           || borromean->crosses()[1].adj_arc[3] != 1
           || borromean->getSign(2) <= 0
           || borromean->crosses()[2].adj_arc[0] != 1
           || borromean->crosses()[2].adj_arc[1] != 2
           || borromean->crosses()[2].adj_arc[2] != 5
           || borromean->crosses()[2].adj_arc[3] != 6
           || borromean->getSign(3) >= 0
           || borromean->crosses()[3].adj_arc[0] != 9
           || borromean->crosses()[3].adj_arc[1] != 10
           || borromean->crosses()[3].adj_arc[2] != 2
           || borromean->crosses()[3].adj_arc[3] != 3
           || borromean->getSign(4) >= 0
           || borromean->crosses()[4].adj_arc[0] != 4
           || borromean->crosses()[4].adj_arc[1] != 5
           || borromean->crosses()[4].adj_arc[2] != 8
           || borromean->crosses()[4].adj_arc[3] != 9
           || borromean->getSign(5) <= 0
           || borromean->crosses()[5].adj_arc[0] != 6
           || borromean->crosses()[5].adj_arc[1] != 7
           || borromean->crosses()[5].adj_arc[2] != 10
           || borromean->crosses()[5].adj_arc[3] != 11
            ) {
            ERR_MSG("Failed to load Borromean ring.");
            return EXIT_FAILURE;
        }

        std::cout << "Pass" << std::endl;
    }

    // True-lover's knot (8_19) as a non-alternating knot.
    {
        auto truelover = read_gauss_code(
            {-1, 2, -3, -4, 5, -6, 7, -8, 4, -5, 6, 1, -2, -7, 8, 3},
            {std::make_pair(1,false)}
            );

        if(!truelover) {
            ERR_MSG("Failed to load true lover's knot.");
            return EXIT_FAILURE;
        }

        // The diagram has only negative crossings.
        for(std::size_t c = 0; c < truelover->ncrosses(); ++c) {
            if (truelover->getSign(c) >= 0) {
                ERR_MSG("Wrong sign of crossings in truelover.");
                return EXIT_FAILURE;
            }
        }
        std::cout << std::endl;
    }

    // Another code for true-lover's knot.
    {
        auto truelover2 = read_gauss_code(
            {1, -8, 2, -1, -4, 5, 8, -2, -3, 7, -6, 4, -5, 3, -7, 6},
            {std::make_pair(1,false)}
            );

        if(!truelover2) {
            ERR_MSG("Failed to load true lover's knot.");
            return EXIT_FAILURE;
        }

        // The diagram has only negative crossings.
        for(std::size_t c = 0; c < truelover2->ncrosses(); ++c) {
            if (truelover2->getSign(c) >= 0) {
                ERR_MSG("Wrong sign of crossings in truelover2.");
                return EXIT_FAILURE;
            }
        }
        std::cout << std::endl;
    }

    return EXIT_SUCCESS;
}

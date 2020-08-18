#include <iostream>

#include "linkdiagram.hpp"

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
            || trefoil->crosses()[0].is_positive != false
            || trefoil->crosses()[0].adj_arc[0] != 5
            || trefoil->crosses()[0].adj_arc[1] != 0
            || trefoil->crosses()[0].adj_arc[2] != 2
            || trefoil->crosses()[0].adj_arc[3] != 3
            || trefoil->crosses()[1].is_positive != false
            || trefoil->crosses()[1].adj_arc[0] != 1
            || trefoil->crosses()[1].adj_arc[1] != 2
            || trefoil->crosses()[1].adj_arc[2] != 4
            || trefoil->crosses()[1].adj_arc[3] != 5
            || trefoil->crosses()[2].is_positive != false
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
           || six_two->crosses()[0].is_positive != false
           || six_two->crosses()[0].adj_arc[0] != 11
           || six_two->crosses()[0].adj_arc[1] != 0
           || six_two->crosses()[0].adj_arc[2] != 4
           || six_two->crosses()[0].adj_arc[3] != 5
           || six_two->crosses()[1].is_positive != false
           || six_two->crosses()[1].adj_arc[0] != 5
           || six_two->crosses()[1].adj_arc[1] != 6
           || six_two->crosses()[1].adj_arc[2] != 0
           || six_two->crosses()[1].adj_arc[3] != 1
           || six_two->crosses()[2].is_positive != false
           || six_two->crosses()[2].adj_arc[0] != 3
           || six_two->crosses()[2].adj_arc[1] != 4
           || six_two->crosses()[2].adj_arc[2] != 10
           || six_two->crosses()[2].adj_arc[3] != 11
           || six_two->crosses()[3].is_positive != false
           || six_two->crosses()[3].adj_arc[0] != 9
           || six_two->crosses()[3].adj_arc[1] != 10
           || six_two->crosses()[3].adj_arc[2] != 6
           || six_two->crosses()[3].adj_arc[3] != 7
           || six_two->crosses()[4].is_positive != true
           || six_two->crosses()[4].adj_arc[0] != 1
           || six_two->crosses()[4].adj_arc[1] != 2
           || six_two->crosses()[4].adj_arc[2] != 8
           || six_two->crosses()[4].adj_arc[3] != 9
           || six_two->crosses()[5].is_positive != true
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
           || six_two->crosses()[0].is_positive != true
           || six_two->crosses()[0].adj_arc[0] != 11
           || six_two->crosses()[0].adj_arc[1] != 0
           || six_two->crosses()[0].adj_arc[2] != 4
           || six_two->crosses()[0].adj_arc[3] != 5
           || six_two->crosses()[1].is_positive != true
           || six_two->crosses()[1].adj_arc[0] != 5
           || six_two->crosses()[1].adj_arc[1] != 6
           || six_two->crosses()[1].adj_arc[2] != 0
           || six_two->crosses()[1].adj_arc[3] != 1
           || six_two->crosses()[2].is_positive != true
           || six_two->crosses()[2].adj_arc[0] != 3
           || six_two->crosses()[2].adj_arc[1] != 4
           || six_two->crosses()[2].adj_arc[2] != 10
           || six_two->crosses()[2].adj_arc[3] != 11
           || six_two->crosses()[3].is_positive != true
           || six_two->crosses()[3].adj_arc[0] != 9
           || six_two->crosses()[3].adj_arc[1] != 10
           || six_two->crosses()[3].adj_arc[2] != 6
           || six_two->crosses()[3].adj_arc[3] != 7
           || six_two->crosses()[4].is_positive != false
           || six_two->crosses()[4].adj_arc[0] != 1
           || six_two->crosses()[4].adj_arc[1] != 2
           || six_two->crosses()[4].adj_arc[2] != 8
           || six_two->crosses()[4].adj_arc[3] != 9
           || six_two->crosses()[5].is_positive != false
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
           || six_two->crosses()[0].is_positive != true
           || six_two->crosses()[0].adj_arc[0] != 11
           || six_two->crosses()[0].adj_arc[1] != 0
           || six_two->crosses()[0].adj_arc[2] != 4
           || six_two->crosses()[0].adj_arc[3] != 5
           || six_two->crosses()[1].is_positive != true
           || six_two->crosses()[1].adj_arc[0] != 5
           || six_two->crosses()[1].adj_arc[1] != 6
           || six_two->crosses()[1].adj_arc[2] != 0
           || six_two->crosses()[1].adj_arc[3] != 1
           || six_two->crosses()[2].is_positive != true
           || six_two->crosses()[2].adj_arc[0] != 3
           || six_two->crosses()[2].adj_arc[1] != 4
           || six_two->crosses()[2].adj_arc[2] != 10
           || six_two->crosses()[2].adj_arc[3] != 11
           || six_two->crosses()[3].is_positive != true
           || six_two->crosses()[3].adj_arc[0] != 9
           || six_two->crosses()[3].adj_arc[1] != 10
           || six_two->crosses()[3].adj_arc[2] != 6
           || six_two->crosses()[3].adj_arc[3] != 7
           || six_two->crosses()[4].is_positive != false
           || six_two->crosses()[4].adj_arc[0] != 1
           || six_two->crosses()[4].adj_arc[1] != 2
           || six_two->crosses()[4].adj_arc[2] != 8
           || six_two->crosses()[4].adj_arc[3] != 9
           || six_two->crosses()[5].is_positive != false
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
           || borromean->crosses()[0].is_positive != false
           || borromean->crosses()[0].adj_arc[0] != 3
           || borromean->crosses()[0].adj_arc[1] != 0
           || borromean->crosses()[0].adj_arc[2] != 7
           || borromean->crosses()[0].adj_arc[3] != 4
           || borromean->crosses()[1].is_positive != true
           || borromean->crosses()[1].adj_arc[0] != 11
           || borromean->crosses()[1].adj_arc[1] != 8
           || borromean->crosses()[1].adj_arc[2] != 0
           || borromean->crosses()[1].adj_arc[3] != 1
           || borromean->crosses()[2].is_positive != true
           || borromean->crosses()[2].adj_arc[0] != 1
           || borromean->crosses()[2].adj_arc[1] != 2
           || borromean->crosses()[2].adj_arc[2] != 5
           || borromean->crosses()[2].adj_arc[3] != 6
           || borromean->crosses()[3].is_positive != false
           || borromean->crosses()[3].adj_arc[0] != 9
           || borromean->crosses()[3].adj_arc[1] != 10
           || borromean->crosses()[3].adj_arc[2] != 2
           || borromean->crosses()[3].adj_arc[3] != 3
           || borromean->crosses()[4].is_positive != false
           || borromean->crosses()[4].adj_arc[0] != 4
           || borromean->crosses()[4].adj_arc[1] != 5
           || borromean->crosses()[4].adj_arc[2] != 8
           || borromean->crosses()[4].adj_arc[3] != 9
           || borromean->crosses()[5].is_positive != true
           || borromean->crosses()[5].adj_arc[0] != 6
           || borromean->crosses()[5].adj_arc[1] != 7
           || borromean->crosses()[5].adj_arc[2] != 10
           || borromean->crosses()[5].adj_arc[3] != 11
            ) {
            std::cerr << "Failed to load Borromean ring." << std::endl;
            return -1;
        }

        std::cout << "Pass" << std::endl;
    }

    return 0;
}

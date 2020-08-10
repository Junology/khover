#include <iostream>

#include "abelian.hpp"

using namespace khover;

int main(int argc, char* argv[])
{
    {
        AbelianGroup a(
            (AbelianGroup::matrix_t(3,4) << 0,1,2,3,4,5,6,7,8,9,10,11).finished(),
            std::false_type{},
            3);
        auto shouldbe = (Eigen::Matrix<int64_t,5,1>() << 4,0,0,0,0).finished();

        a.reduce(std::tuple<>{}, std::tuple<>{});
        if(a.ngens() != 5
           || a.nrels() != 1
           || a.get_repmatrix() != shouldbe)
        {
            std::cout << a.get_repmatrix() << std::endl;
            return -1;
        }
    }
    {
        auto mat = (Eigen::Matrix<int64_t,4,5>() <<
                    6, 0,-4,-4, 2,
                    0, 1, 1, 2, 1,
                    2, 1, 5, 6, 7,
                    8, 2, 2, 4,10 ).finished();
        auto shouldbe = (Eigen::Matrix<int64_t,5,2>() <<
                         2, 0, 0, 16, 0, 0, 0, 0, 0, 0).finished();

        AbelianGroup a(mat, std::false_type{}, 2);

        a.reduce(std::tuple<>{}, std::tuple<>{});

        if(a.ngens() != 5
           || a.nrels() != 2
           || a.get_repmatrix() != shouldbe)
        {
            std::cout << a.get_repmatrix() << std::endl;
            return -1;
        }
    }
    return 0;
}

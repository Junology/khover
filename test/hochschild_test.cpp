#include <iostream>
#include <map>

#include "chaincomplex.hpp"
#include "debug/debug.hpp"

/*! \file hochschild_test.cpp
 * In this test, the program computes the Hochshild homology of the ring Z[x]/(x-1)^n.
 * As it is isomorphic to Z[x]/(x)^n, the result can be found in pp.304 Exercise 9.1.4 in [Weibel1994].
 * Namely, HH_p(Z[x]/(x)^n) is
 *   - Z[x]/(x)^n if p = 0;
 *   - Z[x]/(x^n,nx^{n-1}) if p:odd;
 *   - <x,x^2,...,x^{n-1}> if p:even > 0.
 */

using namespace khover;

using matrix_t = Eigen::Matrix<int64_t,Eigen::Dynamic,Eigen::Dynamic>;
using coeff_t = int64_t;

//! Miscellaneous functions
template <class T>
constexpr T binom(T n, T k)
{
    if (k < 0 || k > n)
        return 0;

    if (k > n-k)
        k = n-k;

    if (k == 0)
        return 1;

    T result = 1;

    for (T i = 0; i < k; ++i)
    {
        result *= n-i;
        result /= (i+1);
    }

    return result;
}

template <std::size_t n>
std::uint64_t adic_digit(std::size_t k, std::uint64_t num)
{
    return (num / cipow(n,k)) % n;
}

//! The representation matrix of the multiplication by x.
template <std::size_t n>
Eigen::Matrix<coeff_t,n,n> mulX() noexcept
{
    Eigen::Matrix<coeff_t,n,n> result = decltype(result)::Zero();

    result.block(1,0,n-1,n-1) = Eigen::Matrix<coeff_t,n-1,n-1>::Identity();
    for(size_t i = 0; i < n; ++i) {
        result.coeffRef(i, n-1) = ((n-i) & 0x1 ? 1 : -1)*binom(n,i);
    }

    return result;
}

template <std::size_t n>
std::map<uint64_t,coeff_t> face(std::size_t deg, std::size_t kth, uint64_t num)
{
    const static Eigen::Matrix<coeff_t,n,n> matX = mulX<n>();

    std::map<uint64_t,coeff_t> result;
    std::size_t knext, num_base;
    if (kth >= deg) {
        knext = 0;
        num_base = (num % cipow(n,deg) / n) * n;
    }
    else {
        knext = kth+1;
        num_base = (num / cipow(n,kth+2)) * cipow(n,kth+1) + num % cipow(n,kth);
    }

    uint64_t r = adic_digit<n>(kth,num) + adic_digit<n>(knext,num);

    if (r < n) {
        result.emplace(num_base + r*cipow(n,kth % deg), 1);
    }
    else {
        auto m = cipow(matX,r-n+1,decltype(matX)::Identity().eval());
        for(std::size_t i = 0; i < n; ++i) {
            result.emplace(num_base + i*cipow(n,kth % deg), m.coeff(i,n-1));
        }
    }

    return result;
}

//! Compute the i-th differential of the Hochschild chain.
template <std::size_t n>
matrix_t get_diff(std::size_t deg)
{
    matrix_t result = matrix_t::Zero(cipow(n,deg+1),cipow(n,deg+2));

    // Traverse (deg+1)-chains
    for(std::size_t s = 0; s < cipow(n,deg+2); ++s) {
        for(std::size_t k = 0; k <= deg+1; ++k) {
            coeff_t sign = (k & 0x1) ? -1 : 1;
            for(auto coeff : face<n>(deg+1, k, s)) {
                result(coeff.first, s) += sign * coeff.second;
            }
        }
    }

    return result;
}

int main(int argc, char* argv[])
{
    constexpr const std::size_t n = 3;

    ChainIntegral ch(
        0,
        get_diff<n>(0),
        get_diff<n>(1),
        get_diff<n>(2),
        get_diff<n>(3)
        );

    auto h = ch.compute();

    std::cout << "---" << std::endl;
    std::cout << "Computing H_0..." << std::endl;
    auto h0 = h[0].compute();
    std::cout << "Free rank: " << h0.first << std::endl;
    std::cout << "Torsions: " << h0.second << std::endl;
    if (h0.first != 3 || !h0.second.empty() != 0) {
        ERR_MSG("Wrong 0th homology");
        return -1;
    }

    std::cout << "---" << std::endl;
    std::cout << "Computing H_1..." << std::endl;
    auto h1 = h[1].compute();
    std::cout << "Free rank: " << h1.first << std::endl;
    std::cout << "Torsions: " << h1.second << std::endl;
    if (h1.first != 2 || h1.second.size() != 1 || h1.second[0] != n) {
        ERR_MSG("Wrong 1st homology");
        return -1;
    }

    std::cout << "---" << std::endl;
    std::cout << "Computing H_2..." << std::endl;
    auto h2 = h[2].compute();
    std::cout << "Free rank: " << h2.first << std::endl;
    std::cout << "Torsions: " << h2.second << std::endl;
    if (h2.first != 2 || !h2.second.empty()) {
        ERR_MSG("Wrong 2nd homology");
        return -1;
    }

    std::cout << "---" << std::endl;
    std::cout << "Computing H_3..." << std::endl;
    auto h3 = h[3].compute();
    std::cout << "Free rank: " << h3.first << std::endl;
    std::cout << "Torsions: " << h3.second << std::endl;
    if (h3.first != 2 || h3.second.size() != 1 || h3.second[0] != n) {
        ERR_MSG("Wrong 3rd homology");
        return -1;
    }

    return 0;
}

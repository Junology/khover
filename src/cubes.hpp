/*!
 * \file cubes.hpp
 * \author Jun Yoshida
 * \copyright (c) 2020 Jun Yoshida.
 * The project is released under the 2-clause BSD License.
 * \date August, 2020: created
 */

#pragma once

#include <vector>

#include "states.hpp"
#include "linkdiagram.hpp"

namespace khover {

//! Base class for cubes of smoothings.
template <template <class...> class C>
class Cube {
public:
    using smoothdata_t = C<std::vector<component_t>>;
    using degree_t = std::int64_t;

protected:
    degree_t m_mincohdeg;
    std::size_t m_dim;
    std::vector<smoothdata_t> m_smoothdata{};

    Cube(degree_t mincohdeg, std::size_t dim)
        : m_mincohdeg(mincohdeg), m_dim(dim)
    {
        m_smoothdata.reserve(cipow(2,dim));
    }

public:
    Cube(Cube const&) = default;
    Cube(Cube&&) noexcept = default;
    ~Cube() noexcept = default;

    Cube& operator=(Cube const&) = default;
    Cube& operator=(Cube&&) noexcept = default;

    //! Array index access.
    inline decltype(auto) operator[](std::size_t i) { return m_smoothdata[i]; }

    //! Array index access.
    inline decltype(auto) operator[](std::size_t i) const { return m_smoothdata[i]; }
    inline decltype(auto) front() { return m_smoothdata.front(); }
    inline decltype(auto) front() const { return m_smoothdata.front(); }
    inline decltype(auto) back() { return m_smoothdata.back(); }
    inline decltype(auto) back() const { return m_smoothdata.back(); }

    inline auto size() const noexcept { return m_smoothdata.size(); }
    inline auto dim() const noexcept { return m_dim; }

    inline state_t minState(std::size_t popcnt) const noexcept {
        popcnt = std::min(popcnt, m_dim);
        return low_window<max_crosses>(popcnt > m_dim ? m_dim : popcnt);
    }

    inline state_t maxState(std::size_t popcnt) const noexcept {
        popcnt = std::min(popcnt, m_dim);
        return low_window<max_crosses>(popcnt)<<(m_dim-popcnt);
    }

    inline degree_t mincohdeg() const noexcept { return m_mincohdeg; }
    inline degree_t maxcohdeg() const noexcept {
        return m_mincohdeg + static_cast<degree_t>(m_dim);
    }

    inline degree_t cohdegree(state_t st) const noexcept {
        return m_mincohdeg + static_cast<degree_t>(st.count());
    }
};

//! The class for cubes of ordinary smoothings.
class SmoothCube : public Cube<Smoothing>
{
protected:
    using Cube<Smoothing>::Cube;

public:
    //! Generate the smoothing cube from a diagram.
    //! \param fixed_states A list of fixed states on crossings.
    //! \retval smt For a *state* st, smt[st] carries the data of components in the smoothing according to st together with fixed states; here, by *state* we mean states on all crossings but ones with fixed states.
    static
    SmoothCube
    fromDiagram(
        LinkDiagram const& diagram,
        std::initializer_list<std::pair<std::size_t,bool>> fixed_states= {}
        ) noexcept;
};

//! The class for crux cubes.
class CruxCube : public Cube<SmoothTw>
{
protected:
    using Cube<SmoothTw>::Cube;

public:
    //! Generate the crux cube from a diagram.
    //! In contrast to SmoothCube::fromDiagram function, this function does not compute components for non-crux states and keep it empty.
    //! \param dblpt The index of a double point.
    //! \retval smt For a *state* st, smt[st] carries the data of components in the smoothing according to st together with that of twisted arcs.
    static
    CruxCube
    fromDiagram(LinkDiagram const& diagram, std::size_t dblpt)
        noexcept;
};

} // end namespace khover


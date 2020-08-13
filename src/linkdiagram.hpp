/*!
 * \file linkdiagram.hpp
 * \author Jun Yoshida
 * \copyright (c) 2020 Jun Yoshida.
 * The project is released under the 2-clause BSD License.
 * \date August, 2020: created
 */

#include <vector>
#include <optional>

//* Debug
#include <iostream>
// */

namespace khover {

//! Internal representation of link diagrams.
//! It only contains minimum information for "state-sum" methods.
class LinkDiagram {
public:
    struct Crossing {
        //! A flag if the crossing is positive or not.
        bool is_positive;
        //! Indicating how the crossing is adjacent to arcs so that
        //!   - *adj_arc[0]* is the in-coming over-arc;
        //!   - *adj_arc[1]* is the out-going over-arc;
        //!   - *adj_arc[2]* is the in-coming under-arc;
        //!   - *adj_arc[3]* is the out-going under-arc.
        //! Here the variable assumes arcs are indexed by unsigned integers.
        unsigned int adj_arc[4];
    };

private:
    std::size_t m_numarcs;
    std::vector<Crossing> m_cross;

public:
    //! Constructor.
    //! This is the only chance to set the number of arcs (except copies).
    LinkDiagram(std::size_t numarcs, std::vector<Crossing> const& cross) noexcept
        : m_numarcs(numarcs), m_cross(cross)
    {}

    LinkDiagram(std::size_t numarcs, std::vector<Crossing>&& cross) noexcept
        : m_numarcs(numarcs), m_cross(std::move(cross))
    {}

    ~LinkDiagram() = default;

    //! The class is copy constructible.
    LinkDiagram(LinkDiagram const&) = default;
    //! The class is move constructible.
    LinkDiagram(LinkDiagram&&) = default;
    //! The class is copyable.
    LinkDiagram& operator=(LinkDiagram const&) = default;
    //! The class is movable.
    LinkDiagram& operator=(LinkDiagram &&) = default;

    std::size_t get_numarcs() const noexcept { return m_numarcs; }
    std::size_t get_numcrosses() const noexcept { return m_cross.size(); }

    const std::vector<Crossing>& crosses() const noexcept { return m_cross; }
};

//! Type used to represent Gauss codes.
using gcode_t = std::vector<int>;

/*!
 * Read a Gauss code.
 * One can pass a multi-component code by separating codes with '0'.
 * \param gcode A Gauss code to process.
 * \param signs A list of signs of the crossings with the minimum indices in *disjoint* components, where disjoint means they are disjoint in the 2d-plane (not in the 3d-space). Hence, only one entry suffices in usual cases.
 */
std::optional<LinkDiagram> read_gauss_code(
    gcode_t const& gcode,
    std::vector<bool> const& signs
    ) noexcept;

} // end namespace khover

#ifndef JGAP_DESCRIPTOR_HPP
#define JGAP_DESCRIPTOR_HPP

#include <array>


namespace jgap {

    /// @brief A wrapper around a fixed-sized array, aimed at emphasizing that the array contains descriptor info.
    template<size_t Dim>
    requires (Dim > 0)
    using Descriptor = std::array<double, Dim>;

}

#endif

#ifndef JGAP_POTENTIALLOADER_HPP
#define JGAP_POTENTIALLOADER_HPP

#include <string>
#include <vector>
#include "jgap/core/ValuePtr.hpp"
#include "jgap/core/potentials/Potential.hpp"

namespace jgap {
    /// Loads a potential from a file path or base name (handling .jgap.h5, .h5, .tabgap.h5, .eam.fs).
    ValuePtr<Potential> loadPotential(const std::string& pot_path);

    /// Loads a potential from a list of file paths (e.g. tabgap.h5 + eam.fs).
    ValuePtr<Potential> loadPotential(const std::vector<std::string>& pot_paths);
}

#endif

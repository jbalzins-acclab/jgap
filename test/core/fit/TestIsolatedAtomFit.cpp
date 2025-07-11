#include <gtest/gtest.h>

#include <utility>

#include "core/fit/IsolatedAtomFit.hpp"
#include "core/potentials/IsolatedAtomPotential.hpp"
#include "utils/Utils.hpp"

jgap::AtomicStructure makeIsolated(jgap::Species species) {
    return jgap::AtomicStructure{
        .atoms = {
            jgap::AtomData{
                .species = std::move(species),
                .position = {0,0,0}
            }
        },
        .lattice = {
            jgap::Vector3{30.0, 0.0, 0.0},
            jgap::Vector3{0.0, 30.0, 0.0},
            jgap::Vector3{0.0, 0.0, 30.0},
        }
    };
}

TEST(TestIsolatedAtom, testCrMnFeNi650) {

    string xyzDataFn = "test/resources/xyz-samples/iter-3-3-train.xyz";
    auto structs = jgap::readXyz(xyzDataFn);

    auto fit = jgap::IsolatedAtomFit(nlohmann::json::parse("{}"));
    auto pot = fit.fit(structs);

    ASSERT_NEAR(pot->predict(makeIsolated("Cr")).energy.value(), -5.4323394, 1e-9);
    ASSERT_NEAR(pot->predict(makeIsolated("Mn")).energy.value(), -5.14356606, 1e-9);
    ASSERT_NEAR(pot->predict(makeIsolated("Fe")).energy.value(), -3.16884796, 1e-9);
    ASSERT_NEAR(pot->predict(makeIsolated("Ni")).energy.value(), -0.70201695, 1e-9);
}

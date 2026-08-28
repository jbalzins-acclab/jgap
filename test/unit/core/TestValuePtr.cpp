#include <gtest/gtest.h>
#include "jgap/core/ValuePtr.hpp"
#include "jgap/core/cutoff/CosCutoff.hpp"
#include "jgap/core/cutoff/CutoffFunction.hpp"
#include "jgap/core/cutoff/PerriotPolynomialCutoff.hpp"

using namespace jgap;

TEST(TestValuePtr, DirectRValueAssignment) {
    ValuePtr<CutoffFunction> c;
    EXPECT_FALSE(bool(c));

    // Direct r-value assignment
    c = CosCutoff(4.5, 1.0);
    EXPECT_TRUE(bool(c));
    EXPECT_DOUBLE_EQ(c->getCutoff(), 4.5);

    auto* cos_ptr = c.as<CosCutoff>();
    ASSERT_NE(cos_ptr, nullptr);
    EXPECT_DOUBLE_EQ(cos_ptr->getCutoff(), 4.5);

    // Direct assignment with another derived type
    c = PerriotPolynomialCutoff(6.0, 0.8);
    EXPECT_DOUBLE_EQ(c->getCutoff(), 6.0);
    EXPECT_NE(c.as<PerriotPolynomialCutoff>(), nullptr);
    EXPECT_EQ(c.as<CosCutoff>(), nullptr);
}

TEST(TestValuePtr, CrossTypeMoveAndCopyAssignment) {
    ValuePtr<CosCutoff> cos_ptr = CosCutoff(5.0, 2.0);
    ValuePtr<CutoffFunction> base_ptr;

    // Cross-type copy assignment
    base_ptr = cos_ptr;
    EXPECT_TRUE(bool(base_ptr));
    EXPECT_TRUE(bool(cos_ptr));
    EXPECT_DOUBLE_EQ(base_ptr->getCutoff(), 5.0);
    EXPECT_NE(base_ptr.get(), cos_ptr.get()); // deep clone

    // Cross-type move assignment
    ValuePtr<CutoffFunction> moved_base_ptr;
    CosCutoff* original_raw = cos_ptr.get();
    moved_base_ptr = std::move(cos_ptr);

    EXPECT_FALSE(bool(cos_ptr));
    EXPECT_TRUE(bool(moved_base_ptr));
    EXPECT_EQ(moved_base_ptr.get(), original_raw); // stole pointer without reallocating
}

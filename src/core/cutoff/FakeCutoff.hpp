#ifndef JGAP_FAKECUTOFF_HPP
#define JGAP_FAKECUTOFF_HPP

#include "CutoffFunction.hpp"
#include "utils/Utils.hpp"

namespace jgap {

    class FakeCutoff : public CutoffFunction {
    public:
        FakeCutoff(const double r) : r_(r) {}

        double evaluate(double r) override { return 1.0; }
        double differentiate(double r) override { return 0.0; }

        std::string getType() override { IO_NOT_INTENDED(FakeCutoff.getType); }
        DataNode serialize() override { IO_NOT_INTENDED(FakeCutoff.serialize); }
        double getCutoff() override { return r_; }

    private:
        double r_;
    };
}

#endif
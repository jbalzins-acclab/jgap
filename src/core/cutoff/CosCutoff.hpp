#ifndef JGAP_DEFAULTCUTOFFFUNCTION_HPP
#define JGAP_DEFAULTCUTOFFFUNCTION_HPP

#include "CutoffFunction.hpp"
#include "io/parse/ParserRegistry.hpp"

namespace jgap {
    class CosCutoff : public CutoffFunction {
    public:
        CosCutoff(double cutoff, double cutoff_transition_width);
        explicit CosCutoff(const DataNode& params);

        static constexpr const char* TYPE = "coscutoff";

        std::string getType() override {return TYPE;};
        DataNode serialize() override;
        double getCutoff() override { return cutoff_; }

        ~CosCutoff() override = default;

        double evaluate(double r) override;
        double differentiate(double r) override;

    private:
        double cutoff_;
        double cutoff_transition_width_;

        double cutoff_transition_width_inverse_;
    };

    SETUP_PARSER(CutoffFunction, CosCutoff)
}
#endif

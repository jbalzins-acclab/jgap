#ifndef JGAP_CUTOFFFUNCTION_HPP
#define JGAP_CUTOFFFUNCTION_HPP

#include <string>
#include "../DataNode.hpp"

namespace jgap {

    class CutoffFunction {
    public:
        virtual ~CutoffFunction() = default;

        virtual double evaluate(double r) = 0;
        virtual double differentiate(double r) = 0;

        virtual double getCutoff() = 0;
    };
}

#endif

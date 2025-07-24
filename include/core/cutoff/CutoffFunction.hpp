#ifndef CUTOFFFUNCTION_HPP
#define CUTOFFFUNCTION_HPP

#include <string>
#include <nlohmann/json.hpp>

using namespace std;

namespace jgap {

    class CutoffFunction {
    public:
        virtual ~CutoffFunction() = default;

        virtual double evaluate(double r) = 0;
        virtual double differentiate(double r) = 0;

        virtual string getType() = 0;
        virtual nlohmann::json serialize() = 0;
        virtual double getCutoff() = 0;
    };
}

#endif //CUTOFFFUNCTION_HPP

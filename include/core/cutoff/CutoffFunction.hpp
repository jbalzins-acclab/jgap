#ifndef CUTOFFFUNCTION_HPP
#define CUTOFFFUNCTION_HPP

namespace jgap {

    class CutoffFunction {
    public:
        virtual ~CutoffFunction() = default;
        virtual double evaluate(double r) = 0;
        virtual double differentiate(double r) = 0;
    };

    class DefaultCutoffFunction : public CutoffFunction {
    public:
        DefaultCutoffFunction(const double cutoff, const double cutoffTransitionWidth)
            : _cutoff(cutoff),
              _cutoffTransitionWidth(cutoffTransitionWidth),
              _cutoffTransitionWidthInverse(1.0 / cutoffTransitionWidth) {
        }
        ~DefaultCutoffFunction() override = default;

        double evaluate(double r) override;
        double differentiate(double r) override;

    private:
        double _cutoff;
        double _cutoffTransitionWidth;
        double _cutoffTransitionWidthInverse;
    };
}

#endif //CUTOFFFUNCTION_HPP

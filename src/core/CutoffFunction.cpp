#include "../../include/core/cutoff/CutoffFunction.hpp"

#include <cmath>

namespace jgap {
    double DefaultCutoffFunction::differentiate(const double r) {
        if (r <= _cutoff - _cutoffTransitionWidth || r >= _cutoff) {
            return 0;
        }
        return -0.5 * M_PI * _cutoffTransitionWidthInverse
                * sin(M_PI*(r - _cutoff + _cutoffTransitionWidth) * _cutoffTransitionWidthInverse) ;
    }

    double DefaultCutoffFunction::evaluate(const double r) {
        /*
          elemental function coordination_function_upper(r,cutoff_in,transition_width)

              real(dp)             :: coordination_function_upper
              real(dp), intent(in) :: r, cutoff_in, transition_width

              if( r > cutoff_in ) then
                  coordination_function_upper = 0.0_dp
              elseif( r > (cutoff_in-transition_width) ) then
                  coordination_function_upper = 0.5_dp * ( cos(PI*(r-cutoff_in+transition_width)/transition_width) + 1.0_dp )
              else
                  coordination_function_upper = 1.0_dp
              endif

           endfunction coordination_function_upper
        */
        if (r <= _cutoff - _cutoffTransitionWidth) {
            return 1;
        }
        if (r >= _cutoff) {
            return 0;
        }
        return 0.5 * (cos(M_PI * (r - _cutoff + _cutoffTransitionWidth) * _cutoffTransitionWidthInverse) + 1);
    }
}

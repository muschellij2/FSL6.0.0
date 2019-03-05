/*
 * Data structure returned by the Base::LogDeterminant method.
 * Copied more or less verbatim from newmat8.cpp.
 */

#ifndef __LOGANDSIGN_HPP__
#define __LOGANDSIGN_HPP__

#include <math.h>

namespace armawrap {

  template<typename elem_type>
  class AWLogAndSign {

    elem_type log_value;
    int       sign;

  public:

    inline AWLogAndSign() { 
      log_value = 0; 
      sign      = 1; 
    }

    inline AWLogAndSign(elem_type f) {

      if (f == 0.0) { 
        log_value = 0;
        sign      = 0; 
        return; 
      }
      else if (f < 0) { 
        sign = -1; 
        f    = -f; 
      }
      else 
        sign = 1;
   
      log_value = log(f);
    }

    inline AWLogAndSign(elem_type f, int s) {
      log_value = f;
      sign      = s;
    }

    inline void operator*=(elem_type x) {
      if (x > 0.0) { 
        log_value += log(x); 
      }
      else if (x < 0.0) { 
        log_value += log(-x); sign = -sign; 
      }
      else 
        sign = 0;
    }

    inline void PowEq(int k) {
      if (sign) {
        log_value *= k;
        if ( (k & 1) == 0 ) sign = 1;
      }
    }

    inline void ChangeSign() { 
      sign = -sign; 
    }

    inline elem_type LogValue() const { 
      return log_value; 
    }

    inline int Sign() const { 
      return sign; 
    }

    inline elem_type Value() const {
      //Tracer et("LogAndSign::Value");
      //if (log_value >= FloatingPointPrecision::LnMaximum())
      //  Throw(OverflowException("Overflow in exponential"));
      return sign * exp(log_value);
    }
  };
}

#endif /* __LOGANDSIGN_HPP__ */

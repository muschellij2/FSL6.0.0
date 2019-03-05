#ifndef __FUNCTION_QR_HPP__
#define __FUNCTION_QR_HPP__

namespace armawrap {
  template<typename T1, typename T2>
  inline void QRZ(T1 &X, T2 &U) {

    arma::Mat<typename T1::elem_type> O;

    arma::qr_econ(O, U, X);

    X = O;
  }
}

#endif /* __FUNCTION_QR_HPP__ */

#ifndef __FUNCTION_FFT_HPP__
#define __FUNCTION_FFT_HPP__

namespace armawrap {
  template<typename xT, typename yT, typename fT, typename gT>
  inline void FFT(const xT &X, 
                  const yT &Y,
                  fT       &F, 
                  gT       &G) {

    arma::Mat<arma::cx_double> XY(X, Y);
    arma::Mat<arma::cx_double> FG = arma::fft(XY);
    F = arma::real(FG);
    G = arma::imag(FG);
  }

  template<typename xT, typename yT, typename fT, typename gT>
  inline void FFTI(const fT &F, 
                   const gT &G,
                   xT       &X, 
                   yT       &Y) {

    arma::Mat<arma::cx_double> FG(F, G);
    arma::Mat<arma::cx_double> XY = arma::ifft(FG);
    X = arma::real(XY);
    Y = arma::imag(XY);
  }

  template<typename xT, typename fT, typename gT>
  inline void RealFFT(const xT &X, 
                      fT       &F, 
                      gT       &G) {

    arma::Col<typename xT::elem_type> Xc = X;
    arma::Col<arma::cx_double> FG = arma::fft(Xc);

    // Newmat only returns the first (X / 2) + 1 elements
    F = arma::real(FG).eval().submat(0, 0, Xc.n_elem / 2, 0);
    G = arma::imag(FG).eval().submat(0, 0, Xc.n_elem / 2, 0);
  }

  template<typename xT, typename fT, typename gT>
  inline void RealFFTI(const fT &F, 
                       const gT &G, 
                       xT       &X) {

    arma::uword nvals = (F.n_elem - 1) * 2;

    arma::Col<typename fT::elem_type> Fc = F;
    arma::Col<typename gT::elem_type> Gc = G;

    arma::Col<typename fT::elem_type> Ffull(nvals);
    arma::Col<typename gT::elem_type> Gfull(nvals);

    // Because the newmat RealFFT function only returns the 
    // first (X / 2) + 1 elements, we need to reconstruct 
    // the full complex vectors
    Ffull.subvec(0,         Fc.n_elem - 1) = Fc;
    Gfull.subvec(0,         Gc.n_elem - 1) = Gc;
    Ffull.subvec(Fc.n_elem, nvals     - 1) =  arma::flipud(Fc.subvec(1, Fc.n_elem - 2));
    Gfull.subvec(Gc.n_elem, nvals     - 1) = -arma::flipud(Gc.subvec(1, Gc.n_elem - 2));

    arma::Mat<arma::cx_double> FG(Ffull, Gfull);
    arma::Mat<arma::cx_double> XY = arma::ifft(FG);

    X = arma::real(XY);
  }
}


#endif /* __FUNCTION_FFT_HPP__ */

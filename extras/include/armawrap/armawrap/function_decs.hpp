#ifndef __FUNCTION_DECS_HPP__
#define __FUNCTION_DECS_HPP__

/*
 * Declarations of standalone armawrap functions.
 */

/* Functions which accept armawrap:: types */

namespace armawrap{ 
  template<typename wT> inline typename enable_if<is_armawrap_type<wT>::value,              typename wT::elem_type>  ::type Minimum(             wT &a);
  template<typename wT> inline typename enable_if<is_armawrap_type<wT>::value,              typename wT::elem_type>  ::type Maximum(             wT &a);
  template<typename wT> inline typename enable_if<is_armawrap_type<wT>::value,              typename wT::elem_type>  ::type MinimumAbsoluteValue(wT &a);
  template<typename wT> inline typename enable_if<is_armawrap_type<wT>::value,              typename wT::elem_type>  ::type MaximumAbsoluteValue(wT &a);
  template<typename wT> inline typename enable_if<is_armawrap_type<wT>::value,              typename wT::elem_type>  ::type SumSquare(           wT &a);
  template<typename wT> inline typename enable_if<is_armawrap_type<wT>::value,              typename wT::elem_type>  ::type SumAbsoluteValue(    wT &a);
  template<typename wT> inline typename enable_if<is_armawrap_type<wT>::value,              typename wT::elem_type>  ::type Sum(                 wT &a);
  template<typename wT> inline typename enable_if<is_armawrap_type<wT>::value,              typename wT::elem_type>  ::type Norm1(               wT &a);
  template<typename wT> inline typename enable_if<is_armawrap_type<wT>::value,              typename wT::elem_type>  ::type NormInfinity(        wT &a);
  template<typename wT> inline typename enable_if<is_armawrap_type<wT>::value,              typename wT::elem_type>  ::type NormFrobenius(       wT &a);
  template<typename wT> inline typename enable_if<is_armawrap_type<wT>::value,              typename wT::elem_type>  ::type Trace(               wT &a);
  template<typename wT> inline typename enable_if<is_armawrap_type<wT>::value,              typename wT::elem_type>  ::type Determinant(         wT &a);
  template<typename wT> inline typename enable_if<is_armawrap_type<wT>::value, AWLogAndSign<typename wT::elem_type> >::type LogDeterminant(      wT &a);
  template<typename wT> inline typename enable_if<is_armawrap_type<wT>::value,              typename wT::elem_type>  ::type IsZero(              wT &a);

  /* Those same functions, but now accepting arma:: types */
  template<typename elem_type, typename T> inline elem_type               Minimum(             arma::Base<elem_type, T> const &a);
  template<typename elem_type, typename T> inline elem_type               Maximum(             arma::Base<elem_type, T> const &a);
  template<typename elem_type, typename T> inline elem_type               MinimumAbsoluteValue(arma::Base<elem_type, T> const &a);
  template<typename elem_type, typename T> inline elem_type               MaximumAbsoluteValue(arma::Base<elem_type, T> const &a);
  template<typename elem_type, typename T> inline elem_type               SumSquare(           arma::Base<elem_type, T> const &a);
  template<typename elem_type, typename T> inline elem_type               SumAbsoluteValue(    arma::Base<elem_type, T> const &a);
  template<typename elem_type, typename T> inline elem_type               Sum(                 arma::Base<elem_type, T> const &a);
  template<typename elem_type, typename T> inline elem_type               Trace(               arma::Base<elem_type, T> const &a);
  template<typename elem_type, typename T> inline elem_type               Norm1(               arma::Base<elem_type, T> const &a);
  template<typename elem_type, typename T> inline elem_type               NormInfinity(        arma::Base<elem_type, T> const &a);
  template<typename elem_type, typename T> inline elem_type               NormFrobenius(       arma::Base<elem_type, T> const &a);
  template<typename elem_type, typename T> inline elem_type               Determinant(         arma::Base<elem_type, T> const &a);
  template<typename elem_type, typename T> inline AWLogAndSign<elem_type> LogDeterminant(      arma::Base<elem_type, T> const &a);
  template<typename elem_type, typename T> inline bool                    IsZero(              arma::Base<elem_type, T> const &a);

  /*
   * As* functions, used by the AWBase class methods of the
   * same name in base.hpp, and defined in functions_as.hpp.
   */
  template<typename wT> inline typename armawrap_asrow_type_map< wT>::type AsRow(     const wT &a);
  template<typename wT> inline typename armawrap_ascol_type_map< wT>::type AsColumn(  const wT &a);
  template<typename wT> inline typename armawrap_asdiag_type_map<wT>::type AsDiagonal(const wT &a);
  template<typename wT> inline typename armawrap_asmat_type_map< wT>::type AsMatrix(  const wT &a, const int nrows, const int ncols);

  /*
   * Miscellaneous internally-used functions.
   */
  inline void FixSubviewIndices(int &first_row, int &last_row, int &first_col, int &last_col);
}

#endif /* __FUNCTION_DECS_HPP__ */

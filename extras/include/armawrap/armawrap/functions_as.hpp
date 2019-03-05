#ifndef __FUNCTION_AS_HPP__
#define __FUNCTION_AS_HPP__


namespace armawrap {


  template<typename wT>
  struct AsRowS {
    static typename armawrap_asrow_type_map<wT>::type AsRow(const wT &a) {
      return arma::vectorise<
        typename armawrap_type_map<wT>::type>(
          a.get_at_ref(), 1);
    } 
  };
 
  template<typename wT>
  struct AsColS {
    static typename armawrap_ascol_type_map<wT>::type AsColumn(const wT &a) {
      return arma::vectorise<
        typename armawrap_type_map<wT>::type>(
          a.get_at_ref(), 0);
    } 
  }; 


  template<typename wT>
  struct AsDiagS {
  
    static typename armawrap_asdiag_type_map<wT>::type
    AsDiagonal(const wT &a) {

      // The inner expression is created on the
      // heap, and deleted by the outer AWOp
      // object - see the AWOp class in exprs.hpp
      arma::Op<typename armawrap_type_map<wT>::type,
               arma::op_vectorise_all> *inner =
        new arma::Op<typename armawrap_type_map<wT>::type,
                     arma::op_vectorise_all>(
                       a.get_at_ref(), 1);

      return AWOp<
        arma::Op<typename armawrap_type_map<wT>::type, arma::op_vectorise_all>,
        arma::op_diagmat>(true, *inner);
    }
  };


  
  template<typename wT>
  struct AsMatS {
    static typename armawrap_asmat_type_map<wT>::type
    AsMatrix(const wT &a, const int nrows, const int ncols) {

      arma::Op<typename armawrap_type_map<wT>::type,
               arma::op_vectorise_all> *inner =
        new arma::Op<typename armawrap_type_map<wT>::type,
                     arma::op_vectorise_all>(
                       a.get_at_ref(), 1);

      return AWOp<
        arma::Op<typename armawrap_type_map<wT>::type, arma::op_vectorise_all>,
        arma::op_reshape>(true, *inner, nrows, ncols); 
    }
  };


  template<typename elem_type>
  struct AsRowS<AWDiagonalMatrix<elem_type> > {
    
    static inline
    typename armawrap_asrow_type_map<AWDiagonalMatrix<elem_type> >::type
    AsRow(const AWDiagonalMatrix<elem_type> &a) {

      arma::Col<elem_type> *inner = new arma::Col<elem_type>(a.get_at_ref().diag());
      return AWOp<arma::Col<elem_type>, arma::op_htrans>(true, *inner);
    }
  };

  
  template<typename elem_type>
  struct AsColS<AWDiagonalMatrix<elem_type> > { 
    static inline typename armawrap_ascol_type_map<AWDiagonalMatrix<elem_type> >::type
    AsColumn(const AWDiagonalMatrix<elem_type> &a) {
      return a.get_at_ref().diag();
    }
  };

  template<typename elem_type>
  struct AsMatS<AWDiagonalMatrix<elem_type> > {

    static inline typename armawrap_asmat_type_map<AWDiagonalMatrix<elem_type> >::type
    AsMatrix(const AWDiagonalMatrix<elem_type> &a, const int nrows, const int ncols) {
      arma::Col<elem_type> *inner = new arma::Col<elem_type>(a.get_at_ref().diag());
      return AWOp<arma::Col<elem_type>, arma::op_reshape>(true, *inner, nrows, ncols);
    }
  };
  

  template<typename wT> inline
  typename armawrap_asrow_type_map< wT>::type AsRow(const wT &a) {
    return AsRowS<wT>::AsRow(a);
  }
  
  template<typename wT> inline
  typename armawrap_ascol_type_map< wT>::type AsColumn(const wT &a) {
    return AsColS<wT>::AsColumn(a);
  }
  
  template<typename wT> inline
  typename armawrap_asdiag_type_map<wT>::type AsDiagonal(const wT &a) {
    return AsDiagS<wT>::AsDiagonal(a);
  }
  
  template<typename wT> inline
  typename armawrap_asmat_type_map< wT>::type AsMatrix(
    const wT &a, const int nrows, const int ncols) {
    return AsMatS<wT>::AsMatrix(a, nrows, ncols);
  } 
}

#endif /* __FUNCTION_AS_HPP__ */

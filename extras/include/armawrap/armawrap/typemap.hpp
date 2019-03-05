#ifndef __TYPEMAP_HPP__
#define __TYPEMAP_HPP__

/*
 * Type mappings from arma::types to corresponding armawrap::types.
 */


namespace armawrap {

  /*
   * Equivalent of std::enable_if (only available in C++11) for 
   * template function matching.
   */
  template<bool, typename result_type> struct enable_if                    {                           };
  template<      typename result_type> struct enable_if<true, result_type> { typedef result_type type; };

  /*
   * Template macros for determining whether one/two 
   * specified types are armawrap::types.
   */
  template<typename T> 
  struct is_armawrap_type { 
    static const bool value = false;
  };

  template<typename elem_type> 
  struct is_armawrap_type<AWMatrix<elem_type> > {
    static const bool value = true;
  };

  template<typename elem_type> 
  struct is_armawrap_type<AWColVector<elem_type> > {
    static const bool value = true;
  };

  template<typename elem_type> 
  struct is_armawrap_type<AWRowVector<elem_type> > {
    static const bool value = true;
  };  

  template<typename elem_type> 
  struct is_armawrap_type<AWDiagonalMatrix<elem_type> > {
    static const bool value = true;
  };

  template<typename elem_type> 
  struct is_armawrap_type<AWIdentityMatrix<elem_type> > {
    static const bool value = true;
  };

  template<typename elem_type> 
  struct is_armawrap_type<AWSymmetricMatrix<elem_type> > {
    static const bool value = true;
  };

  template<typename elem_type> 
  struct is_armawrap_type<AWUpperTriangularMatrix<elem_type> > {
    static const bool value = true;
  };

  template<typename elem_type> 
  struct is_armawrap_type<AWLowerTriangularMatrix<elem_type> > {
    static const bool value = true;
  };

  template<typename elem_type> 
  struct is_armawrap_type<AWBandMatrix<elem_type> > {
    static const bool value = true;
  };

  template<typename elem_type> 
  struct is_armawrap_type<AWLowerBandMatrix<elem_type> > {
    static const bool value = true;
  };
  
  template<typename elem_type> 
  struct is_armawrap_type<AWUpperBandMatrix<elem_type> > {
    static const bool value = true;
  }; 

  template<typename elem_type> 
  struct is_armawrap_type<AWSymmetricBandMatrix<elem_type> > {
    static const bool value = true;
  };

  template<typename elem_type> 
  struct is_armawrap_type<AWInvertedCroutMatrix<elem_type> > {
    static const bool value = true;
  }; 

  template<typename elem_type, typename aw_type> 
  struct is_armawrap_type<AWSubView<elem_type, aw_type> > {
    static const bool value = true;
  };

  template<typename T1, typename T2, typename glue_type> 
  struct is_armawrap_type<AWGlue<T1, T2, glue_type> > {
    static const bool value = true;
  };

  template<typename T1, typename T2, typename eglue_type> 
  struct is_armawrap_type<AWEGlue<T1, T2, eglue_type> > {
    static const bool value = true;
  };

  template<typename T, typename eop_type> 
  struct is_armawrap_type<AWEOp<T, eop_type> > {
    static const bool value = true;
  };

  template<typename T, typename op_type> 
  struct is_armawrap_type<AWOp<T, op_type> > {
    static const bool value = true;
  }; 

  template<typename T1, typename T2>
  struct either_are_armawrap_type {
    static const bool value = is_armawrap_type<T1>::value || 
                              is_armawrap_type<T2>::value;
  };

  template<typename T1, typename T2>
  struct both_are_arma_type {
    static const bool value = arma::is_arma_type<T1>::value &&
                              arma::is_arma_type<T2>::value;
  };

  /*
   * The AWMatrix, AWColVector, and AWRowVector are deemed to be 'simple' 
   * 'simple' armawrap types, because .... These
   * macros are used by the call operator (see operator_call.hpp) to 
   * distinguish between 'simple' types and other types.
   */
  template<typename T>
  struct is_simple_armawrap_type {
    static const bool value = false;
  };

  template<typename elem_type>
  struct is_simple_armawrap_type<AWMatrix<elem_type> > {
    static const bool value = true;
  };

  template<typename elem_type>
  struct is_simple_armawrap_type<AWRowVector<elem_type> > {
    static const bool value = true;
  };

  template<typename elem_type>
  struct is_simple_armawrap_type<AWColVector<elem_type> > {
    static const bool value = true;
  };

  template<typename elem_type>
  struct is_simple_armawrap_type<AWDiagonalMatrix<elem_type> > {
    static const bool value = true;
  }; 

  /*
   * Template containers which map from armawrap::types 
   * to the corresponding arma::types.
   */
  template<typename T> struct armawrap_type_map;

  template<typename elem_type>
  struct armawrap_type_map<AWMatrix<elem_type> > {
    typedef arma::Mat<elem_type> type; 
  };

  template<typename elem_type>
  struct armawrap_type_map<AWRowVector<elem_type> > {
    typedef arma::Mat<elem_type> type; 
  };

  template<typename elem_type>
  struct armawrap_type_map<AWColVector<elem_type> > {
    typedef arma::Mat<elem_type> type; 
  };  

  template<typename elem_type>
  struct armawrap_type_map<AWDiagonalMatrix<elem_type> > {
    typedef arma::Mat<elem_type> type; 
  }; 

  template<typename elem_type>
  struct armawrap_type_map<AWIdentityMatrix<elem_type> > {
    typedef arma::Mat<elem_type> type; 
  }; 

  template<typename elem_type>
  struct armawrap_type_map<AWSymmetricMatrix<elem_type> > {
    typedef arma::Mat<elem_type> type; 
  }; 

  template<typename elem_type>
  struct armawrap_type_map<AWUpperTriangularMatrix<elem_type> > {
    typedef arma::Mat<elem_type> type; 
  }; 

  template<typename elem_type>
  struct armawrap_type_map<AWLowerTriangularMatrix<elem_type> > {
    typedef arma::Mat<elem_type> type; 
  };

  template<typename elem_type>
  struct armawrap_type_map<AWBandMatrix<elem_type> > {
    typedef arma::Mat<elem_type> type; 
  };

  template<typename elem_type>
  struct armawrap_type_map<AWLowerBandMatrix<elem_type> > {
    typedef arma::Mat<elem_type> type; 
  };
 
  template<typename elem_type>
  struct armawrap_type_map<AWUpperBandMatrix<elem_type> > {
    typedef arma::Mat<elem_type> type; 
  };

  template<typename elem_type>
  struct armawrap_type_map<AWSymmetricBandMatrix<elem_type> > {
    typedef arma::Mat<elem_type> type; 
  };

  template<typename elem_type>
  struct armawrap_type_map<AWInvertedCroutMatrix<elem_type> > {
    typedef arma::Mat<elem_type> type; 
  }; 

  template<typename elem_type, typename aw_type>
  struct armawrap_type_map<AWSubView<elem_type, aw_type> > {
    typedef arma::subview<elem_type> type;
  }; 

  template<typename T1, typename T2, typename glue_type>
  struct armawrap_type_map<AWGlue<T1, T2, glue_type> > {
    typedef arma::Glue<T1, T2, glue_type> type;
  }; 

  template<typename T1, typename T2, typename eglue_type>
  struct armawrap_type_map<AWEGlue<T1, T2, eglue_type> > {
    typedef arma::eGlue<T1, T2, eglue_type> type;
  }; 


  template<typename T, typename eop_type>
  struct armawrap_type_map<AWEOp<T, eop_type> > {
    typedef arma::eOp<T, eop_type> type;
  };


  template<typename T, typename op_type>
  struct armawrap_type_map<AWOp<T, op_type> > {
    typedef arma::Op<T, op_type> type;
  };

  /* 
   * Below are dummy template macros which map each of the major
   * arma::types to themselves, thus allowing all of said
   * arma::types to be passed to the armawrap_type_map template. 
   */

  template<typename elem_type>
  struct armawrap_type_map<arma::Mat<elem_type> > {
    typedef arma::Mat<elem_type> type; 
  }; 

  template<typename elem_type>
  struct armawrap_type_map<arma::Row<elem_type> > {
    typedef arma::Row<elem_type> type; 
  }; 

  template<typename elem_type>
  struct armawrap_type_map<arma::Col<elem_type> > {
    typedef arma::Col<elem_type> type; 
  }; 

  template<typename T1, typename T2, typename glue_type>
  struct armawrap_type_map<arma::Glue<T1, T2, glue_type> > {
    typedef arma::Glue<T1, T2, glue_type> type; 
  }; 

  template<typename T1, typename T2, typename eglue_type>
  struct armawrap_type_map<arma::eGlue<T1, T2, eglue_type> > {
    typedef arma::eGlue<T1, T2, eglue_type> type; 
  }; 

  template<typename T1, typename T2, typename mtglue_type, typename eT>
  struct armawrap_type_map<arma::mtGlue<eT, T1, T2, mtglue_type> > {
    typedef arma::mtGlue<eT, T1, T2, mtglue_type> type; 
  }; 

  template<typename T, typename op_type>
  struct armawrap_type_map<arma::Op<T, op_type> > {
    typedef arma::Op<T, op_type> type; 
  }; 

  template<typename T, typename eop_type>
  struct armawrap_type_map<arma::eOp<T, eop_type> > {
    typedef arma::eOp<T, eop_type> type; 
  }; 

  template<typename elem_type>
  struct armawrap_type_map<arma::subview<elem_type> > {
    typedef arma::subview<elem_type> type;
  }; 

  template<typename elem_type>
  struct armawrap_type_map<arma::subview_row<elem_type> > {
    typedef arma::subview_row<elem_type> type;
  }; 

  template<typename elem_type>
  struct armawrap_type_map<arma::subview_col<elem_type> > {
    typedef arma::subview_col<elem_type> type;
  };



  template<typename wT>
  struct armawrap_asrow_type_map {
    typedef arma::Op<typename armawrap_type_map<wT>::type, arma::op_vectorise_all> type;
  };

  template<typename wT>
  struct armawrap_ascol_type_map {
    typedef arma::Op<typename armawrap_type_map<wT>::type, arma::op_vectorise_all> type;
  };

  template<typename wT>
  struct armawrap_asdiag_type_map {
    typedef AWOp<
      arma::Op<typename armawrap_type_map<wT>::type,
               arma::op_vectorise_all>,
      arma::op_diagmat> type;
  };

  template<typename wT>
  struct armawrap_asmat_type_map {
    typedef AWOp<
      arma::Op<typename armawrap_type_map<wT>::type,
               arma::op_vectorise_all>,
      arma::op_reshape> type;
  };

  template<typename elem_type>
  struct armawrap_asrow_type_map<AWDiagonalMatrix<elem_type> > {
    typedef AWOp<arma::Col<elem_type>, arma::op_htrans> type;
  };

  template<typename elem_type>
  struct armawrap_ascol_type_map<AWDiagonalMatrix<elem_type> > {
    typedef arma::Col<elem_type> type;
  };

  template<typename elem_type>
  struct armawrap_asmat_type_map<AWDiagonalMatrix<elem_type> > {
    typedef AWOp<arma::Col<elem_type>, arma::op_reshape> type;
  }; 
}

#endif

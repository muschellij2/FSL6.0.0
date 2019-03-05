#ifndef __CLASS_DECS_HPP__
#define __CLASS_DECS_HPP__

/*
 * Forward declarations for all armawrap classes.
 */


namespace armawrap {
  template<typename elem_type, typename wT, typename aT>     class AWBase;
  template<typename elem_type>                               class AWMatrix;
  template<typename elem_type>                               class AWIdentityMatrix;
  template<typename elem_type>                               class AWDiagonalMatrix;
  template<typename elem_type>                               class AWSymmetricMatrix;
  template<typename elem_type>                               class AWUpperTriangularMatrix;
  template<typename elem_type>                               class AWLowerTriangularMatrix;
  template<typename elem_type>                               class AWBandMatrix;
  template<typename elem_type>                               class AWUpperBandMatrix;
  template<typename elem_type>                               class AWLowerBandMatrix;
  template<typename elem_type>                               class AWSymmetricBandMatrix;
  template<typename elem_type>                               class AWCroutMatrix;
  template<typename elem_type>                               class AWInvertedCroutMatrix;
  template<typename elem_type>                               class AWRowVector;
  template<typename elem_type>                               class AWColVector;
  template<typename elem_type, typename wT>                  class AWInsertManager;
  template<typename elem_type, typename wT>                  class AWCallManagerBase; 
  template<typename elem_type, typename wT, bool simpleType> class AWCallManager;
  template<typename elem_type, typename wT>                  class AWSubView;
  template<typename elem_type>                               class AWLogAndSign;

  template<typename T1, typename T2, typename eglue_type> class AWEGlue;
  template<typename T1, typename T2, typename glue_type>  class AWGlue;
  template<typename T,               typename eop_type>   class AWEOp;
  template<typename T,               typename op_type>    class AWOp;

  class AWMatrixBandWidth;
  class AWException;
}

#endif /* __CLASS_DECS_HPP__ */

#ifndef __OPERATOR_INSERT_HPP__
#define __OPERATOR_INSERT_HPP__

namespace armawrap {
  
  template<typename elem_type, typename wT>
  inline void doInsert(wT &target, unsigned int count, elem_type val) {
    target(count) = val;
  }

  template<typename elem_type>
  inline void doInsert(AWIdentityMatrix<elem_type> &target, unsigned int count, elem_type val) {
    target = val;
  }

  /**
   * Insert values via << into a armawrap::AWBase object.
   */
  template<typename elem_type, typename wT>
  class AWInsertManager {

    wT  & target;
    int   count;

  public:

    inline AWInsertManager(wT & target, elem_type firstVal) : 

      target(target), count(0) {
      target.zeros();
      doInsert(target, ++count, firstVal);

    }

    inline AWInsertManager & operator<<(elem_type val) {
      doInsert(target, ++count, val);
      return *this;
    }

    inline ~AWInsertManager() {

      if (count != target.Storage()) {
        throw AWException("A list of values was too short");
      }
    }
  };
}

#endif /* __OPERATOR_INSERT_HPP__ */

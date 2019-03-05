#ifndef __EXCEPTION_HPP__
#define __EXCEPTION_HPP__

/*
 * Exception class(es).
 */

#include <exception>

namespace armawrap {
  class AWException : public std::runtime_error {
  public:
    AWException(const char *msg) : std::runtime_error(msg) {};
  };
}

#endif /* __EXCEPTION_HPP__ */

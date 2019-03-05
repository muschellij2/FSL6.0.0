#ifndef __ARMAWRAP_HPP__
#define __ARMAWRAP_HPP__


#define __ARMAWRAP_VERSION__ 0.2.1


/*
 * armawrap. A library which looks like Newmat, but is actually armadillo.
 * Does that make any sense?
 */

#include "armadillo"

namespace armawrap {}

// Forward declarations of armawrap classes
#include "class_decs.hpp"

// Exception classes
#include "exception.hpp"

// Template magic for dealing with expressions containing
// both armawrap::types and arma::types.
#include "typemap.hpp"

// Forward declarations of armawrap functions
#include "function_decs.hpp"

// Basic data type definitions
#include "logandsign.hpp"
#include "base.hpp"
#include "matrix.hpp"
#include "vector.hpp"

// Other matrix types
#include "matrix_diag.hpp"
#include "matrix_tri.hpp"
#include "matrix_band.hpp"
#include "matrix_crout.hpp"

// Subviews
#include "subview.hpp"

// Expressions
#include "exprs.hpp"

// Unary and binary operators
#include "operator_call.hpp"
#include "operator_plus.hpp"
#include "operator_minus.hpp"
#include "operator_div.hpp"
#include "operator_times.hpp"
#include "operator_pipe.hpp"
#include "operator_ampersand.hpp"
#include "operator_equality.hpp"

// Data entry via the insertion
// operator, and stream output
#include "operator_insert.hpp"
#include "operator_stream.hpp"

// Miscellaneous functions
#include "function_schur.hpp"
#include "function_kronecker.hpp"
#include "function_dotproduct.hpp"
#include "function_svd.hpp"
#include "function_cholesky.hpp"
#include "function_qr.hpp"
#include "function_eigenvalues.hpp"
#include "function_fft.hpp"
#include "function_sort.hpp"

// Standalone functions on matrix/vector objects
#include "functions_base.hpp"
#include "functions_as.hpp"

#endif /* __ARMAWRAP_HPP__ */

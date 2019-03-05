//  This file was automatically generated on Thu Sep 18 23:32:18 2008
//  by libs/config/tools/generate.cpp
//  Copyright John Maddock 2002-4.
//  Use, modification and distribution are subject to the 
//  Boost Software License, Version 1.0. (See accompanying file 
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

//  See http://www.boost.org/libs/config for the most recent version.//
//  Revision $Id: no_typename_with_ctor_pass.cpp,v 1.1.1.1 2015/02/27 16:50:26 mwebster Exp $
//


// Test file for macro BOOST_NO_TYPENAME_WITH_CTOR
// This file should compile, if it does not then
// BOOST_NO_TYPENAME_WITH_CTOR should be defined.
// See file boost_no_typename_with_ctor.ipp for details

// Must not have BOOST_ASSERT_CONFIG set; it defeats
// the objective of this file:
#ifdef BOOST_ASSERT_CONFIG
#  undef BOOST_ASSERT_CONFIG
#endif

#include <boost/config.hpp>
#include "test.hpp"

#ifndef BOOST_NO_TYPENAME_WITH_CTOR
#include "boost_no_typename_with_ctor.ipp"
#else
namespace boost_no_typename_with_ctor = empty_boost;
#endif

int main( int, char *[] )
{
   return boost_no_typename_with_ctor::test();
}


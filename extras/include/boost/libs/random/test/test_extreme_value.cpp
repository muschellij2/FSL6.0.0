/* test_extreme_value.cpp
 *
 * Copyright Steven Watanabe 2010
 * Distributed under the Boost Software License, Version 1.0. (See
 * accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * $Id: test_extreme_value.cpp,v 1.1.1.1 2015/02/27 16:50:41 mwebster Exp $
 *
 */

#include <boost/random/extreme_value_distribution.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/math/distributions/extreme_value.hpp>

#define BOOST_RANDOM_DISTRIBUTION boost::random::extreme_value_distribution<>
#define BOOST_RANDOM_DISTRIBUTION_NAME extreme_value
#define BOOST_MATH_DISTRIBUTION boost::math::extreme_value
#define BOOST_RANDOM_ARG1_TYPE double
#define BOOST_RANDOM_ARG1_NAME a
#define BOOST_RANDOM_ARG1_DEFAULT 1000.0
#define BOOST_RANDOM_ARG1_DISTRIBUTION(n) boost::uniform_real<>(0.00001, n)
#define BOOST_RANDOM_ARG2_TYPE double
#define BOOST_RANDOM_ARG2_NAME b
#define BOOST_RANDOM_ARG2_DEFAULT 1000.0
#define BOOST_RANDOM_ARG2_DISTRIBUTION(n) boost::uniform_real<>(0.00001, n)

#include "test_real_distribution.ipp"

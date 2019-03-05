
#ifndef BOOST_FSM_EVENT_INCLUDED
#define BOOST_FSM_EVENT_INCLUDED

// Copyright Aleksey Gurtovoy 2002-2004
//
// Distributed under the Boost Software License, Version 1.0. 
// (See accompanying file LICENSE_1_0.txt or copy at 
// http://www.boost.org/LICENSE_1_0.txt)
//
// See http://www.boost.org/libs/mpl for documentation.

// $Id: event.hpp,v 1.1.1.2 2015/02/27 16:50:37 mwebster Exp $
// $Date: 2015/02/27 16:50:37 $
// $Revision: 1.1.1.2 $

#include "base_event.hpp"

namespace fsm { namespace aux {

template< typename Derived >
struct event
    : base_event
{
 public:
    typedef base_event base_t;

 private:
    virtual std::auto_ptr<base_event> do_clone() const
    {
        return std::auto_ptr<base_event>(
              new Derived(static_cast<Derived const&>(*this))
            );
    }
};

}}

#endif // BOOST_FSM_EVENT_INCLUDED

/***************************************************************************
                          GtestMissingError.h  -  description
                             -------------------
    begin                : Jul 2, 2017
    copyright            : (C) 2017 by Tomas Oberhuber et al.
    email                : tomas.oberhuber@fjfi.cvut.cz
 ***************************************************************************/

/* See Copyright Notice in tnl/Copyright */

#pragma once

#include <stdexcept>

struct GtestMissingError
   : public std::runtime_error
{
   GtestMissingError()
   : std::runtime_error( "The GTest library is needed to run the tests." )
   {}
};

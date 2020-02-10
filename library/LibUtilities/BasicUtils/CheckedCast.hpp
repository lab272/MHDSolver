///////////////////////////////////////////////////////////////////////////////
//
// File CheckedCast.hpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the "Software"),
// to deal in the Software without restriction, including without limitation
// the rights to use, copy, modify, merge, publish, distribute, sublicense,
// and/or sell copies of the Software, and to permit persons to whom the
// Software is furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
// THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
// DEALINGS IN THE SOFTWARE.
//
// Description: simple routines to check if the casting is narrowing
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIB_UTILITIES_CHECKEDCAST_H
#define NEKTAR_LIB_UTILITIES_CHECKEDCAST_H

#include <LibUtilities/LibUtilitiesDeclspec.h>
#include <LibUtilities/BasicUtils/ErrorUtil.hpp>
#include <limits>

namespace Nektar
{
namespace LibUtilities
{

/// checked cast for int types to float types
LIB_UTILITIES_EXPORT template <class To, class Ti>
inline To checked_cast(const Ti param)
{
    Ti min = std::numeric_limits<To>::min();
    Ti max = std::numeric_limits<To>::max();
    ASSERTL0(param >= min, "Casting would narrow (underflow).");
    ASSERTL0(param <= max, "Casting would narrow (overflow).");
    return static_cast<To>(param);
}

// do not allow for int types to int types conversion
LIB_UTILITIES_EXPORT template <class To>
To checked_cast(int param) = delete;

LIB_UTILITIES_EXPORT template <class To>
To checked_cast(long param) = delete;

LIB_UTILITIES_EXPORT template <class To>
To checked_cast(unsigned int param) = delete;

LIB_UTILITIES_EXPORT template <class To>
To checked_cast(unsigned long param) = delete;

}
}

#endif

#ifndef __CIFTI_VERSION_H__
#define __CIFTI_VERSION_H__

/*LICENSE_START*/ 
/*
 *  Copyright (c) 2014, Washington University School of Medicine
 *  All rights reserved.
 *
 *  Redistribution and use in source and binary forms, with or without modification,
 *  are permitted provided that the following conditions are met:
 *
 *  1. Redistributions of source code must retain the above copyright notice,
 *  this list of conditions and the following disclaimer.
 *
 *  2. Redistributions in binary form must reproduce the above copyright notice,
 *  this list of conditions and the following disclaimer in the documentation
 *  and/or other materials provided with the distribution.
 *
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 *  THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 *  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 *  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 *  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 *  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 *  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 *  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 *  EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include "AString.h"

#include "stdint.h"

namespace cifti
{
    class CiftiVersion
    {
        int16_t m_major, m_minor;
    public:
        int16_t getMajor() const { return m_major; }
        int16_t getMinor() const { return m_minor; }
        
        CiftiVersion();
        CiftiVersion(const int16_t& major, const int16_t& minor);
        CiftiVersion(const AString& versionString);
        AString toString() const;
        bool operator<(const CiftiVersion& rhs) const;
        bool operator>(const CiftiVersion& rhs) const;
        bool operator==(const CiftiVersion& rhs) const;
        bool operator!=(const CiftiVersion& rhs) const;
        bool operator<=(const CiftiVersion& rhs) const;
        bool operator>=(const CiftiVersion& rhs) const;
        ///quirk tests
        bool hasReversedFirstDims() const;
    };
}

#endif //__CIFTI_VERSION_H__

#ifndef __NIFTI_IO_H__
#define __NIFTI_IO_H__

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

#include "ByteSwapping.h"
#include "CiftiAssert.h"
#include "BinaryFile.h"
#include "CiftiException.h"
#include "NiftiHeader.h"

#include <cmath>
#include <limits>
#include <vector>

namespace cifti
{
    
    class NiftiIO
    {
        BinaryFile m_file;
        NiftiHeader m_header;
        std::vector<int64_t> m_dims;
        std::vector<char> m_scratch;//scratch memory for byteswapping, type conversion, etc
        int numBytesPerElem();//for resizing scratch
        template<typename TO, typename FROM>
        void convertRead(TO* out, FROM* in, const int64_t& count);//for reading from file
        template<typename TO, typename FROM>
        void convertWrite(TO* out, const FROM* in, const int64_t& count);//for writing to file
    public:
        void openRead(const AString& filename);
        void writeNew(const AString& filename, const NiftiHeader& header, const int& version = 1, const bool& withRead = false, const bool& swapEndian = false);
        AString getFilename() const { return m_file.getFilename(); }
        void overrideDimensions(const std::vector<int64_t>& newDims) { m_dims = newDims; }//HACK: deal with reading/writing CIFTI-1's broken headers
        void close();
        const NiftiHeader& getHeader() const { return m_header; }
        const std::vector<int64_t>& getDimensions() const { return m_dims; }
        int getNumComponents() const;
        //to read/write 1 frame of a standard volume file, call with fullDims = 3, indexSelect containing indexes for any of dims 4-7 that exist
        //NOTE: you need to provide storage for all components within the range, if getNumComponents() == 3 and fullDims == 0, you need 3 elements allocated
        template<typename T>
        void readData(T* dataOut, const int& fullDims, const std::vector<int64_t>& indexSelect, const bool& tolerateShortRead = false);
        template<typename T>
        void writeData(const T* dataIn, const int& fullDims, const std::vector<int64_t>& indexSelect);
    };
    
    template<typename T>
    void NiftiIO::readData(T* dataOut, const int& fullDims, const std::vector<int64_t>& indexSelect, const bool& tolerateShortRead)
    {
        CiftiAssert(fullDims >= 0 && fullDims <= (int)m_dims.size());
        CiftiAssert((size_t)fullDims + indexSelect.size() == m_dims.size());//could be >=, but should catch more stupid mistakes as ==
        int64_t numElems = getNumComponents();//for now, calculate read size on the fly, as the read call will be the slowest part
        int curDim;
        for (curDim = 0; curDim < fullDims; ++curDim)
        {
            numElems *= m_dims[curDim];
        }
        int64_t numDimSkip = numElems, numSkip = 0;
        for (; curDim < (int)m_dims.size(); ++curDim)
        {
            CiftiAssert(indexSelect[curDim - fullDims] >= 0 && indexSelect[curDim - fullDims] < m_dims[curDim]);
            numSkip += indexSelect[curDim - fullDims] * numDimSkip;
            numDimSkip *= m_dims[curDim];
        }
        m_scratch.resize(numElems * numBytesPerElem());
        m_file.seek(numSkip * numBytesPerElem() + m_header.getDataOffset());
        int64_t numRead = 0;
        m_file.read(m_scratch.data(), m_scratch.size(), &numRead);
        if ((numRead != (int64_t)m_scratch.size() && !tolerateShortRead) || numRead < 0)//for now, assume read giving -1 is always a problem
        {
            throw CiftiException("error while reading from file '" + m_file.getFilename() + "'");
        }
        switch (m_header.getDataType())
        {
            case NIFTI_TYPE_UINT8:
            case NIFTI_TYPE_RGB24://handled by components
                convertRead(dataOut, (uint8_t*)m_scratch.data(), numElems);
                break;
            case NIFTI_TYPE_INT8:
                convertRead(dataOut, (int8_t*)m_scratch.data(), numElems);
                break;
            case NIFTI_TYPE_UINT16:
                convertRead(dataOut, (uint16_t*)m_scratch.data(), numElems);
                break;
            case NIFTI_TYPE_INT16:
                convertRead(dataOut, (int16_t*)m_scratch.data(), numElems);
                break;
            case NIFTI_TYPE_UINT32:
                convertRead(dataOut, (uint32_t*)m_scratch.data(), numElems);
                break;
            case NIFTI_TYPE_INT32:
                convertRead(dataOut, (int32_t*)m_scratch.data(), numElems);
                break;
            case NIFTI_TYPE_UINT64:
                convertRead(dataOut, (uint64_t*)m_scratch.data(), numElems);
                break;
            case NIFTI_TYPE_INT64:
                convertRead(dataOut, (int64_t*)m_scratch.data(), numElems);
                break;
            case NIFTI_TYPE_FLOAT32:
            case NIFTI_TYPE_COMPLEX64://components
                convertRead(dataOut, (float*)m_scratch.data(), numElems);
                break;
            case NIFTI_TYPE_FLOAT64:
            case NIFTI_TYPE_COMPLEX128:
                convertRead(dataOut, (double*)m_scratch.data(), numElems);
                break;
            case NIFTI_TYPE_FLOAT128:
            case NIFTI_TYPE_COMPLEX256:
                convertRead(dataOut, (long double*)m_scratch.data(), numElems);
                break;
            default:
                CiftiAssert(0);
                throw CiftiException("internal error, tell the developers what you just tried to do");
        }
    }
    
    template<typename T>
    void NiftiIO::writeData(const T* dataIn, const int& fullDims, const std::vector<int64_t>& indexSelect)
    {
        CiftiAssert(fullDims >= 0 && fullDims <= (int)m_dims.size());
        CiftiAssert((size_t)fullDims + indexSelect.size() == m_dims.size());//could be >=, but should catch more stupid mistakes as ==
        int64_t numElems = getNumComponents();//for now, calculate read size on the fly, as the read call will be the slowest part
        int curDim;
        for (curDim = 0; curDim < fullDims; ++curDim)
        {
            numElems *= m_dims[curDim];
        }
        int64_t numDimSkip = numElems, numSkip = 0;
        for (; curDim < (int)m_dims.size(); ++curDim)
        {
            CiftiAssert(indexSelect[curDim - fullDims] >= 0 && indexSelect[curDim - fullDims] < m_dims[curDim]);
            numSkip += indexSelect[curDim - fullDims] * numDimSkip;
            numDimSkip *= m_dims[curDim];
        }
        m_scratch.resize(numElems * numBytesPerElem());
        m_file.seek(numSkip * numBytesPerElem() + m_header.getDataOffset());
        switch (m_header.getDataType())
        {
            case NIFTI_TYPE_UINT8:
            case NIFTI_TYPE_RGB24://handled by components
                convertWrite((uint8_t*)m_scratch.data(), dataIn, numElems);
                break;
            case NIFTI_TYPE_INT8:
                convertWrite((int8_t*)m_scratch.data(), dataIn, numElems);
                break;
            case NIFTI_TYPE_UINT16:
                convertWrite((uint16_t*)m_scratch.data(), dataIn, numElems);
                break;
            case NIFTI_TYPE_INT16:
                convertWrite((int16_t*)m_scratch.data(), dataIn, numElems);
                break;
            case NIFTI_TYPE_UINT32:
                convertWrite((uint32_t*)m_scratch.data(), dataIn, numElems);
                break;
            case NIFTI_TYPE_INT32:
                convertWrite((int32_t*)m_scratch.data(), dataIn, numElems);
                break;
            case NIFTI_TYPE_UINT64:
                convertWrite((uint64_t*)m_scratch.data(), dataIn, numElems);
                break;
            case NIFTI_TYPE_INT64:
                convertWrite((int64_t*)m_scratch.data(), dataIn, numElems);
                break;
            case NIFTI_TYPE_FLOAT32:
            case NIFTI_TYPE_COMPLEX64://components
                convertWrite((float*)m_scratch.data(), dataIn, numElems);
                break;
            case NIFTI_TYPE_FLOAT64:
            case NIFTI_TYPE_COMPLEX128:
                convertWrite((double*)m_scratch.data(), dataIn, numElems);
                break;
            case NIFTI_TYPE_FLOAT128:
            case NIFTI_TYPE_COMPLEX256:
                convertWrite((long double*)m_scratch.data(), dataIn, numElems);
                break;
            default:
                CiftiAssert(0);
                throw CiftiException("internal error, tell the developers what you just tried to do");
        }
        m_file.write(m_scratch.data(), m_scratch.size());
    }
    
    template<typename TO, typename FROM>
    void NiftiIO::convertRead(TO* out, FROM* in, const int64_t& count)
    {
        if (m_header.isSwapped())
        {
            ByteSwapping::swapArray(in, count);
        }
        double mult, offset;
        bool doScale = m_header.getDataScaling(mult, offset);
        if (std::numeric_limits<TO>::is_integer)//do round to nearest when integer output type
        {
            if (doScale)
            {
                for (int64_t i = 0; i < count; ++i)
                {
                    out[i] = (TO)floor(0.5 + offset + mult * (long double)in[i]);//we don't always need that much precision, but it will still be faster than hard drives
                }
            } else {
                for (int64_t i = 0; i < count; ++i)
                {
                    out[i] = (TO)floor(0.5 + in[i]);
                }
            }
        } else {
            if (doScale)
            {
                for (int64_t i = 0; i < count; ++i)
                {
                    out[i] = (TO)(offset + mult * (long double)in[i]);//we don't always need that much precision, but it will still be faster than hard drives
                }
            } else {
                for (int64_t i = 0; i < count; ++i)
                {
                    out[i] = (TO)in[i];//explicit cast to make sure the compiler doesn't squawk
                }
            }
        }
    }
    
    template<typename TO, typename FROM>
    void NiftiIO::convertWrite(TO* out, const FROM* in, const int64_t& count)
    {
        double mult, offset;
        bool doScale = m_header.getDataScaling(mult, offset);
        if (std::numeric_limits<TO>::is_integer)//do round to nearest when integer output type
        {
            if (doScale)
            {
                for (int64_t i = 0; i < count; ++i)
                {
                    out[i] = (TO)floor(0.5 + ((long double)in[i] - offset) / mult);//we don't always need that much precision, but it will still be faster than hard drives
                }
            } else {
                for (int64_t i = 0; i < count; ++i)
                {
                    out[i] = (TO)floor(0.5 + in[i]);
                }
            }
        } else {
            if (doScale)
            {
                for (int64_t i = 0; i < count; ++i)
                {
                    out[i] = (TO)(((long double)in[i] - offset) / mult);//we don't always need that much precision, but it will still be faster than hard drives
                }
            } else {
                for (int64_t i = 0; i < count; ++i)
                {
                    out[i] = (TO)in[i];//explicit cast to make sure the compiler doesn't squawk
                }
            }
        }
        if (m_header.isSwapped()) ByteSwapping::swapArray(out, count);
    }
    
}

#endif //__NIFTI_IO_H__

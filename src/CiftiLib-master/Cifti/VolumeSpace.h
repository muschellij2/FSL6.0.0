#ifndef __VOLUME_SPACE_H__
#define __VOLUME_SPACE_H__

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

#include "CiftiAssert.h"
#include "Vector3D.h"

#include "XmlAdapter.h"

#include "stdint.h"
#include <vector>

namespace cifti
{

    class VolumeSpace
    {
        int64_t m_dims[3];
        std::vector<std::vector<float> > m_sform, m_inverse;
        void computeInverse();
    public:
        enum OrientTypes
        {
            LEFT_TO_RIGHT = 0,
            RIGHT_TO_LEFT = 4,
            POSTERIOR_TO_ANTERIOR = 1,
            ANTERIOR_TO_POSTERIOR = 5,
            INFERIOR_TO_SUPERIOR = 2,
            SUPERIOR_TO_INFERIOR = 6
        };
        VolumeSpace();
        VolumeSpace(const int64_t dims[3], const std::vector<std::vector<float> >& sform);
        VolumeSpace(const int64_t dims[3], const float sform[12]);
        void setSpace(const int64_t dims[3], const std::vector<std::vector<float> >& sform);
        void setSpace(const int64_t dims[3], const float sform[12]);
        const int64_t* getDims() const { return m_dims; }
        const std::vector<std::vector<float> >& getSform() const { return m_sform; }
        void getSpacingVectors(Vector3D& iStep, Vector3D& jStep, Vector3D& kStep, Vector3D& origin) const;
        bool matchesVolumeSpace(const VolumeSpace& right) const;//allows slight mismatches
        bool operator==(const VolumeSpace& right) const;//requires that it be exact
        bool operator!=(const VolumeSpace& right) const { return !(*this == right); }

        ///returns true if volume space is not skew, and each axis and index is separate
        bool isPlumb() const;

        ///returns orientation, spacing, and center (spacing/center can be negative, spacing/center is LPI rearranged to ijk (first dimension uses first element), will assert false if isOblique is true)
        void getOrientAndSpacingForPlumb(OrientTypes* orientOut, float* spacingOut, float* originOut) const;
        
        ///get just orientation, even for non-plumb volumes
        void getOrientation(OrientTypes orientOut[3]) const;
        
        ///returns coordinate triplet of an index triplet
        template <typename T>
        inline void indexToSpace(const T* indexIn, float* coordOut) const
        { indexToSpace<T>(indexIn[0], indexIn[1], indexIn[2], coordOut[0], coordOut[1], coordOut[2]); }
        
        ///returns coordinate triplet of three indices
        template <typename T>
        inline void indexToSpace(const T& indexIn1, const T& indexIn2, const T& indexIn3, float* coordOut) const
        { indexToSpace<T>(indexIn1, indexIn2, indexIn3, coordOut[0], coordOut[1], coordOut[2]); }
        
        ///returns three coordinates of an index triplet
        template <typename T>
        inline void indexToSpace(const T* indexIn, float& coordOut1, float& coordOut2, float& coordOut3) const
        { indexToSpace<T>(indexIn[0], indexIn[1], indexIn[2], coordOut1, coordOut2, coordOut3); }
        
        ///returns three coordinates of three indices
        template <typename T>
        void indexToSpace(const T& indexIn1, const T& indexIn2, const T& indexIn3, float& coordOut1, float& coordOut2, float& coordOut3) const;

        ///returns floating point index triplet of a given coordinate triplet
        inline void spaceToIndex(const float* coordIn, float* indexOut) const { spaceToIndex(coordIn[0], coordIn[1], coordIn[2], indexOut[0], indexOut[1], indexOut[2]); }
        ///returns floating point index triplet of three given coordinates
        inline void spaceToIndex(const float& coordIn1, const float& coordIn2, const float& coordIn3, float* indexOut) const { spaceToIndex(coordIn1, coordIn2, coordIn3, indexOut[0], indexOut[1], indexOut[2]); }
        ///returns three floating point indexes of a given coordinate triplet
        inline void spaceToIndex(const float* coordIn, float& indexOut1, float& indexOut2, float& indexOut3) const { spaceToIndex(coordIn[0], coordIn[1], coordIn[2], indexOut1, indexOut2, indexOut3); }
        ///returns three floating point indexes of three given coordinates
        void spaceToIndex(const float& coordIn1, const float& coordIn2, const float& coordIn3, float& indexOut1, float& indexOut2, float& indexOut3) const;

        ///returns integer index triplet of voxel whose center is closest to the coordinate triplet
        inline void enclosingVoxel(const float* coordIn, int64_t* indexOut) const { enclosingVoxel(coordIn[0], coordIn[1], coordIn[2], indexOut[0], indexOut[1], indexOut[2]); }
        ///returns integer index triplet of voxel whose center is closest to the three coordinates
        inline void enclosingVoxel(const float& coordIn1, const float& coordIn2, const float& coordIn3, int64_t* indexOut) const { enclosingVoxel(coordIn1, coordIn2, coordIn3, indexOut[0], indexOut[1], indexOut[2]); }
        ///returns integer indexes of voxel whose center is closest to the coordinate triplet
        inline void enclosingVoxel(const float* coordIn, int64_t& indexOut1, int64_t& indexOut2, int64_t& indexOut3) const { enclosingVoxel(coordIn[0], coordIn[1], coordIn[2], indexOut1, indexOut2, indexOut3); }
        ///returns integer indexes of voxel whose center is closest to the three coordinates
        void enclosingVoxel(const float& coordIn1, const float& coordIn2, const float& coordIn3, int64_t& indexOut1, int64_t& indexOut2, int64_t& indexOut3) const;

        template <typename T>
        inline bool indexValid(const T* indexIn) const
        {
            return indexValid(indexIn[0], indexIn[1], indexIn[2]);//implicit cast to int64_t
        }

        ///checks if an index is within array dimensions
        inline bool indexValid(const int64_t& indexIn1, const int64_t& indexIn2, const int64_t& indexIn3) const
        {
            if (indexIn1 < 0 || indexIn1 >= m_dims[0]) return false;
            if (indexIn2 < 0 || indexIn2 >= m_dims[1]) return false;
            if (indexIn3 < 0 || indexIn3 >= m_dims[2]) return false;
            return true;
        }
        
        inline int64_t getIndex(const int64_t& indexIn1, const int64_t& indexIn2, const int64_t& indexIn3) const
        {
            CiftiAssert(indexValid(indexIn1, indexIn2, indexIn3));
            return indexIn1 + m_dims[0] * (indexIn2 + m_dims[1] * indexIn3);
        }
        
        template <typename T>
        inline int64_t getIndex(const T* indexIn) const
        {
            return getIndex(indexIn[0], indexIn[1], indexIn[2]);//implicit cast to int64_t
        }
        
        void readCiftiXML1(XmlReader& xml);//xml functions
        void readCiftiXML2(XmlReader& xml);
        void writeCiftiXML1(XmlWriter& xml) const;
        void writeCiftiXML2(XmlWriter& xml) const;
    };

    template <typename T>
    void VolumeSpace::indexToSpace(const T& indexIn1, const T& indexIn2, const T& indexIn3, float& coordOut1, float& coordOut2, float& coordOut3) const
    {
        coordOut1 = indexIn1 * m_sform[0][0] + indexIn2 * m_sform[0][1] + indexIn3 * m_sform[0][2] + m_sform[0][3];
        coordOut2 = indexIn1 * m_sform[1][0] + indexIn2 * m_sform[1][1] + indexIn3 * m_sform[1][2] + m_sform[1][3];
        coordOut3 = indexIn1 * m_sform[2][0] + indexIn2 * m_sform[2][1] + indexIn3 * m_sform[2][2] + m_sform[2][3];
    }

}

#endif //__VOLUME_SPACE_H__

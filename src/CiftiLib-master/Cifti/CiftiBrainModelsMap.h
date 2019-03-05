#ifndef __CIFTI_BRAIN_MODELS_MAP_H__
#define __CIFTI_BRAIN_MODELS_MAP_H__

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

#include "CiftiMappingType.h"

#include "Compact3DLookup.h"
#include "StructureEnum.h"
#include "VolumeSpace.h"

#include <map>
#include <utility>
#include <vector>

namespace cifti
{
    class CiftiBrainModelsMap : public CiftiMappingType
    {
    public:
        enum ModelType
        {
            SURFACE,
            VOXELS
        };
        struct SurfaceMap
        {
            int64_t m_ciftiIndex;
            int64_t m_surfaceNode;
        };
        struct VolumeMap
        {
            int64_t m_ciftiIndex;
            int64_t m_ijk[3];
        };
        struct ModelInfo
        {
            ModelType m_type;
            StructureEnum::Enum m_structure;
            int64_t m_indexStart, m_indexCount;//these are intended only for summary info, use getSurfaceMap, etc for the index to vertex/voxel mappings
        };
        struct IndexInfo
        {
            ModelType m_type;
            StructureEnum::Enum m_structure;
            int64_t m_surfaceNode;//only one of these two will be valid
            int64_t m_ijk[3];
        };
        bool hasVolumeData() const;
        bool hasVolumeData(const StructureEnum::Enum& structure) const;
        bool hasSurfaceData(const StructureEnum::Enum& structure) const;
        int64_t getIndexForNode(const int64_t& node, const StructureEnum::Enum& structure) const;
        int64_t getIndexForVoxel(const int64_t* ijk, StructureEnum::Enum* structureOut = NULL) const;
        int64_t getIndexForVoxel(const int64_t& i, const int64_t& j, const int64_t& k, StructureEnum::Enum* structureOut = NULL) const;
        IndexInfo getInfoForIndex(const int64_t index) const;
        std::vector<SurfaceMap> getSurfaceMap(const StructureEnum::Enum& structure) const;
        std::vector<VolumeMap> getFullVolumeMap() const;
        std::vector<VolumeMap> getVolumeStructureMap(const StructureEnum::Enum& structure) const;
        const VolumeSpace& getVolumeSpace() const;
        int64_t getSurfaceNumberOfNodes(const StructureEnum::Enum& structure) const;
        std::vector<StructureEnum::Enum> getSurfaceStructureList() const;
        std::vector<StructureEnum::Enum> getVolumeStructureList() const;
        const std::vector<int64_t>& getNodeList(const StructureEnum::Enum& structure) const;//useful for copying mappings to a new dense mapping
        const std::vector<int64_t>& getVoxelList(const StructureEnum::Enum& structure) const;
        std::vector<ModelInfo> getModelInfo() const;
        
        CiftiBrainModelsMap() { m_haveVolumeSpace = false; m_ignoreVolSpace = false; }
        void addSurfaceModel(const int64_t& numberOfNodes, const StructureEnum::Enum& structure, const float* roi = NULL);
        void addSurfaceModel(const int64_t& numberOfNodes, const StructureEnum::Enum& structure, const std::vector<int64_t>& nodeList);
        void addVolumeModel(const StructureEnum::Enum& structure, const std::vector<int64_t>& ijkList);
        void setVolumeSpace(const VolumeSpace& space);
        void clear();
        
        CiftiMappingType* clone() const { return new CiftiBrainModelsMap(*this); }
        MappingType getType() const { return BRAIN_MODELS; }
        int64_t getLength() const;
        bool operator==(const CiftiMappingType& rhs) const;
        bool approximateMatch(const CiftiMappingType& rhs) const;
        void readXML1(XmlReader& xml);
        void readXML2(XmlReader& xml);
        void writeXML1(XmlWriter& xml) const;
        void writeXML2(XmlWriter& xml) const;
    private:
        struct BrainModelPriv
        {
            ModelType m_type;
            StructureEnum::Enum m_brainStructure;
            int64_t m_surfaceNumberOfNodes;
            std::vector<int64_t> m_nodeIndices;
            std::vector<int64_t> m_voxelIndicesIJK;
            
            int64_t m_modelStart, m_modelEnd;//stuff only needed for optimization - models are kept in sorted order by their index ranges
            std::vector<int64_t> m_nodeToIndexLookup;
            bool operator==(const BrainModelPriv& rhs) const;
            bool operator!=(const BrainModelPriv& rhs) const { return !((*this) == rhs); }
            void setupSurface(const int64_t& start);
        };
        VolumeSpace m_volSpace;
        bool m_haveVolumeSpace, m_ignoreVolSpace;//second is needed for parsing cifti-1
        std::vector<BrainModelPriv> m_modelsInfo;
        std::map<StructureEnum::Enum, int> m_surfUsed, m_volUsed;
        Compact3DLookup<std::pair<int64_t, StructureEnum::Enum> > m_voxelToIndexLookup;//make one unified lookup rather than separate lookups per volume structure
        int64_t getNextStart() const;
        struct ParseHelperModel
        {//specifically to allow the parsed elements to be sorted before using addSurfaceModel/addVolumeModel
            ModelType m_type;
            StructureEnum::Enum m_brainStructure;
            int64_t m_surfaceNumberOfNodes;
            std::vector<int64_t> m_nodeIndices;
            std::vector<int64_t> m_voxelIndicesIJK;
            int64_t m_offset, m_count;
            bool operator<(const ParseHelperModel& rhs) const
            {
                if (m_offset < rhs.m_offset) return true;
                if (m_offset > rhs.m_offset) return false;//get the common cases first
                if (m_count < rhs.m_count) return true;//in case we have a zero-length model - this shouldn't happen, usually
                return false;
            }
            void parseBrainModel1(XmlReader& xml);
            void parseBrainModel2(XmlReader& xml);
            static std::vector<int64_t> readIndexArray(XmlReader& xml);
        };
    };
}

#endif //__CIFTI_BRAIN_MODELS_MAP_H__

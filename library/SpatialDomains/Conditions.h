////////////////////////////////////////////////////////////////////////////////
//
//  File:  Conditions.h
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description:
//
//
////////////////////////////////////////////////////////////////////////////////
#ifndef NEKTAR_SPATIALDOMAINS_BOUNDARYCONDITIONS_H
#define NEKTAR_SPATIALDOMAINS_BOUNDARYCONDITIONS_H

#include <string>
#include <map>
#include <iostream>
#include <sstream>

#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <SpatialDomains/MeshGraph.h>
#include <LibUtilities/BasicUtils/Equation.h>
#include <SpatialDomains/SpatialDomainsDeclspec.h>

class TiXmlElement;
class TiXmlDocument;

class SpatialDomainsDeclspec;

namespace Nektar
{
    namespace SpatialDomains
    {
        enum BoundaryConditionType
        {
            eDirichlet,
            eNeumann,
            eRobin,
            ePeriodic,
            eJunction,
            eBifurcation,
            eMerging
        };

        enum BndUserDefinedType
        {
            eI,
            eMG,
            eHigh,
            eWall,
            eWALL,
            eSymmetry,
            eRinglebFlow,
            eTimeDependent,
            eIsentropicVortex,
            eCalcBC,
            eQinflow,
            eNoUserDefined
        };

        struct BoundaryConditionBase
        {
            BoundaryConditionBase(
                    BoundaryConditionType type,
                    const std::string &userDefined = std::string("NoUserDefined"));

            virtual ~BoundaryConditionBase();

            SPATIAL_DOMAINS_EXPORT BoundaryConditionType GetBoundaryConditionType() const;
            SPATIAL_DOMAINS_EXPORT void SetUserDefined(BndUserDefinedType type);
            SPATIAL_DOMAINS_EXPORT BndUserDefinedType GetUserDefined() const;

        protected:
            BoundaryConditionType m_boundaryConditionType;
            BndUserDefinedType    m_userDefined;
        };


        struct DirichletBoundaryCondition : public BoundaryConditionBase
        {

             DirichletBoundaryCondition(
                    const std::string& eqn,
                    const std::string& userDefined = std::string("NoUserDefined"),
                    const std::string& filename=std::string("")):
                BoundaryConditionBase(eDirichlet, userDefined),
                m_dirichletCondition(eqn),
                m_filename(filename)
                {
                }

            LibUtilities::Equation m_dirichletCondition;
            std::string m_filename;
        };

        struct NeumannBoundaryCondition : public BoundaryConditionBase
        {
            NeumannBoundaryCondition(
                    const std::string& eqn,
                    const std::string& userDefined = std::string("NoUserDefined"),
                    const std::string& filename=std::string("")):
                BoundaryConditionBase(eNeumann, userDefined),
                m_neumannCondition(eqn),
                m_filename(filename)
            {
            }

            LibUtilities::Equation m_neumannCondition;
            std::string m_filename;
        };

        struct RobinBoundaryCondition : public BoundaryConditionBase
        {
            RobinBoundaryCondition(
                    const std::string &a,
                    const std::string &b,
                    const std::string &userDefined = std::string("NoUserDefined"),
                    const std::string& filename=std::string("")):
                BoundaryConditionBase(eRobin, userDefined),
                m_robinFunction(a),
                m_robinPrimitiveCoeff(b),
                m_filename(filename)
            {
            }
                // \frac{\partial {u}}{\partial{n}} +
                // m_robinPrimativeCoeff(x,y,z)*u = m_robinFunction(x,y,z)
            LibUtilities::Equation m_robinFunction;
            LibUtilities::Equation m_robinPrimitiveCoeff;
            std::string m_filename;
        };


        struct PeriodicBoundaryCondition : public BoundaryConditionBase
        {
            PeriodicBoundaryCondition(const unsigned int n):
                BoundaryConditionBase(ePeriodic),
                m_connectedBoundaryRegion(n)
            {
            }

            unsigned int m_connectedBoundaryRegion;
        };

        struct JunctionBoundaryCondition : public BoundaryConditionBase
        {
            JunctionBoundaryCondition(
                    const int& P,
                    const int& D1,
                    const std::string& userDefined = std::string("NoUserDefined"));

            SPATIAL_DOMAINS_EXPORT int GetParent() const;
            SPATIAL_DOMAINS_EXPORT int GetDaughter1() const;

            int m_parent;
            int m_daughter1;
        };

        struct BifurcationBoundaryCondition : public BoundaryConditionBase
        {
            BifurcationBoundaryCondition(
                    const int& P,
                    const int& D1,
                    const int& D2,
                    const std::string& userDefined = std::string("NoUserDefined"));

            SPATIAL_DOMAINS_EXPORT int GetParent() const;
            SPATIAL_DOMAINS_EXPORT int GetDaughter1() const;
            SPATIAL_DOMAINS_EXPORT int GetDaughter2() const;

            int m_parent;
            int m_daughter1;
            int m_daughter2;
        };

        struct MergingBoundaryCondition : public BoundaryConditionBase
        {
            MergingBoundaryCondition(
                    const int& P,
                    const int& D1,
                    const int& D2,
                    const std::string& userDefined = std::string("NoUserDefined"));

            SPATIAL_DOMAINS_EXPORT int GetParent() const;
            SPATIAL_DOMAINS_EXPORT int GetDaughter1() const;
            SPATIAL_DOMAINS_EXPORT int GetDaughter2() const;

            int m_parent;
            int m_daughter1;
            int m_daughter2;
        };


        typedef std::map<int, Composite> BoundaryRegion;
        typedef boost::shared_ptr<BoundaryRegion> BoundaryRegionShPtr;
        typedef boost::shared_ptr<const BoundaryRegion> ConstBoundaryRegionShPtr;
        typedef std::vector<BoundaryRegionShPtr> BoundaryRegionCollection;

        typedef boost::shared_ptr<BoundaryConditionBase> BoundaryConditionShPtr;
        typedef boost::shared_ptr<DirichletBoundaryCondition> DirichletBCShPtr;
        typedef boost::shared_ptr<NeumannBoundaryCondition>   NeumannBCShPtr;
        typedef boost::shared_ptr<RobinBoundaryCondition>     RobinBCShPtr;
        typedef boost::shared_ptr<JunctionBoundaryCondition>  JunctionBCShPtr;
        typedef boost::shared_ptr<BifurcationBoundaryCondition>  BifurcationBCShPtr;
        typedef boost::shared_ptr<MergingBoundaryCondition>   MergingBCShPtr;
        typedef std::map<std::string,BoundaryConditionShPtr>  BoundaryConditionMap;
        typedef boost::shared_ptr<BoundaryConditionMap>  BoundaryConditionMapShPtr;
        typedef std::map<int, BoundaryConditionMapShPtr> BoundaryConditionCollection;

        const static Array<OneD, BoundaryConditionShPtr> NullBoundaryConditionShPtrArray;

        class BoundaryConditions
        {
        public:
            SPATIAL_DOMAINS_EXPORT BoundaryConditions(
                    const LibUtilities::SessionReaderSharedPtr& pSession,
                    const MeshGraphSharedPtr&  meshGraph);
            SPATIAL_DOMAINS_EXPORT ~BoundaryConditions();

            const BoundaryRegionCollection &GetBoundaryRegions(void) const;
            const BoundaryConditionCollection &GetBoundaryConditions(void) const;
            const std::string GetVariable(unsigned int indx);

        protected:
            /// The mesh graph to use for referencing geometry info.
            MeshGraphSharedPtr                      m_meshGraph;
            LibUtilities::SessionReaderSharedPtr    m_session;

            BoundaryRegionCollection                m_boundaryRegions;
            BoundaryConditionCollection             m_boundaryConditions;

        private:
            BoundaryConditions();

            /// Read segments (and general MeshGraph) given TiXmlDocument.
            void Read(TiXmlElement *conditions);

            void ReadBoundaryRegions(TiXmlElement *regions);
            void ReadBoundaryConditions(TiXmlElement *conditions);
        };

        typedef boost::shared_ptr<BoundaryConditions> BoundaryConditionsSharedPtr;
    }
}

#endif //NEKTAR_SPATIALDOMAINS_BOUNDARYCONDITIONS_H


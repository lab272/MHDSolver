////////////////////////////////////////////////////////////////////////////////
//
//  File:  BoundaryConditions.h
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

#include <SpatialDomains/Equation.h>
#include <LibUtilities/BasicConst/NektarUnivTypeDefs.hpp>
#include <SpatialDomains/MeshGraph.h>

class TiXmlElement;
class TiXmlDocument;

class MeshGraph;

namespace Nektar
{
    namespace SpatialDomains
    {
        enum BoundaryConditionType
        {
            eDirichlet,
            eNeumann,
            eRobin,
            ePeriodic
        };

        struct BoundaryConditionBase
        {
            BoundaryConditionBase(BoundaryConditionType type):
                m_BoundaryConditionType(type)
            {
            }

            BoundaryConditionBase(BoundaryConditionType type, Equation userDefined):
                m_BoundaryConditionType(type), m_userDefined(userDefined)
            {
            }

            virtual ~BoundaryConditionBase(){};

            BoundaryConditionType GetBoundaryConditionType() const
            {
                return m_BoundaryConditionType;
            }

            void SetUserDefined(Equation equation){
                m_userDefined = equation;
            }
            
            Equation GetUserDefined() const{
                return m_userDefined;
            }


        protected:
            BoundaryConditionType m_BoundaryConditionType;
            Equation m_userDefined;
        };


        struct DirichletBoundaryCondition : public BoundaryConditionBase
        {

            DirichletBoundaryCondition(const std::string &eqn, const std::string &userDefined = std::string("NoUserDefined")):
                BoundaryConditionBase(eDirichlet, userDefined),
                m_DirichletCondition(eqn)
            {
            }
            
            Equation m_DirichletCondition;
        };

        struct NeumannBoundaryCondition : public BoundaryConditionBase
        {
            NeumannBoundaryCondition(const std::string &eqn, const std::string &userDefined = std::string("NoUserDefined")):
                BoundaryConditionBase(eNeumann, userDefined),
                m_NeumannCondition(eqn)
            {
            }

            Equation m_NeumannCondition;
        };

        struct RobinBoundaryCondition : public BoundaryConditionBase
        {
            RobinBoundaryCondition( const std::string &a, const std::string &b, const std::string &userDefined = std::string("NoUserDefined")):
                BoundaryConditionBase(eRobin, userDefined),
                m_a(a), m_b(b)
            {
            }

            // u = a(x,y,z) + b(x,y,z)*\frac{\partial{u}}{\partial{n}}
            Equation m_a;
            Equation m_b;
        };
        

        struct PeriodicBoundaryCondition : public BoundaryConditionBase
        {
            PeriodicBoundaryCondition(const unsigned int n):
                BoundaryConditionBase(ePeriodic),
                m_ConnectedBoundaryRegion(n)
            {
            }
             
            unsigned int m_ConnectedBoundaryRegion;
        };

        typedef std::map<std::string, NekDouble> ParamMap;
        typedef std::map<std::string, std::string> FunctionMap;
        typedef std::vector<std::string> Variable;

        typedef std::vector<Composite> BoundaryRegion;
        typedef boost::shared_ptr<BoundaryRegion> BoundaryRegionShPtr;
        typedef boost::shared_ptr<const BoundaryRegion> ConstBoundaryRegionShPtr;
        typedef std::vector<BoundaryRegionShPtr> BoundaryRegionCollection;

        typedef boost::shared_ptr<BoundaryConditionBase> BoundaryConditionShPtr;
        typedef boost::shared_ptr<DirichletBoundaryCondition> DirichletBCShPtr;
        typedef boost::shared_ptr<NeumannBoundaryCondition>   NeumannBCShPtr;
        typedef boost::shared_ptr<RobinBoundaryCondition>     RobinBCShPtr;
        typedef std::map<std::string,BoundaryConditionShPtr>  BoundaryConditionMap;
        typedef boost::shared_ptr<BoundaryConditionMap>  BoundaryConditionMapShPtr;
        typedef std::map<int, BoundaryConditionMapShPtr> BoundaryConditionCollection;

        const static Array<OneD, BoundaryConditionShPtr> NullBoundaryConditionShPtrArray;

        typedef Equation ForcingFunction;
        typedef boost::shared_ptr<ForcingFunction> ForcingFunctionShPtr;
        typedef boost::shared_ptr<const ForcingFunction> ConstForcingFunctionShPtr;
        typedef std::map<std::string, ForcingFunctionShPtr> ForcingFunctionsMap;

        typedef Equation ExactSolution;
        typedef boost::shared_ptr<ExactSolution> ExactSolutionShPtr;
        typedef boost::shared_ptr<const ExactSolution> ConstExactSolutionShPtr;
        typedef std::map<std::string, ExactSolutionShPtr> ExactSolutionMap;

        typedef Equation UserDefinedEqn;
        typedef boost::shared_ptr<UserDefinedEqn> UserDefinedEqnShPtr;
        typedef boost::shared_ptr<const UserDefinedEqn> ConstUserDefinedEqnShPtr;
        typedef std::map<std::string, UserDefinedEqnShPtr> UserDefinedEqnMap;

        typedef Equation InitialCondition;
        typedef boost::shared_ptr<InitialCondition> InitialConditionShPtr;
        typedef boost::shared_ptr<const InitialCondition> ConstInitialConditionShPtr;
        typedef std::map<std::string, InitialConditionShPtr> InitialConditionsMap;

        typedef std::map<std::string, std::string> SolverInfoMap;

        class BoundaryConditions
        {
        public:
            BoundaryConditions(const MeshGraph *meshGraph);
            ~BoundaryConditions();

            void Read(std::string &infilename);
            void Read(TiXmlDocument &doc);
            
            bool   CheckForParameter(const std::string &paramName);
            static NekDouble GetParameter(const std::string &parmName);

            BoundaryRegionCollection &GetBoundaryRegions(void) 
            {
                return m_BoundaryRegions;
            }

            BoundaryConditionCollection &GetBoundaryConditions(void)
            {
                return m_BoundaryConditions;
            }

            /// Get forcing function based on the index of the variable.
            /// The index is the order in which the variable was
            /// defined.
            ConstForcingFunctionShPtr GetForcingFunction(int indx) const;

            /// Get forcing function based on name of variable.
            ConstForcingFunctionShPtr GetForcingFunction(const std::string &var) const;

            bool ExactSolutionExists(int indx) const;
            ConstExactSolutionShPtr GetExactSolution(int indx) const;
            ConstExactSolutionShPtr GetExactSolution(const std::string &var) const;

            bool UserDefinedEqnExists(const std::string &var) const;
            ConstUserDefinedEqnShPtr GetUserDefinedEqn(int indx) const;
            ConstUserDefinedEqnShPtr GetUserDefinedEqn(const std::string &var) const;
            
            /// Get initial condition function based on the index of the variable.
            /// The index is the order in which the variable was
            /// defined.
            bool InitialConditionExists(int indx) const;
            ConstInitialConditionShPtr GetInitialCondition(int indx) const;

            /// Get initial condition function based on name of variable.
            ConstInitialConditionShPtr GetInitialCondition(const std::string &var) const;
            /// Check to see if initial condition exists in list. 
            bool FoundInitialCondition(const std::string &var);
            
            const std::string &GetVariable(unsigned int indx)
            {
                ASSERTL0(0 <= indx && indx < m_Variables.size(),"Variable index is out of range");
                return m_Variables[indx];
            }

            inline int GetNumVariables() const 
            {
                return m_Variables.size();
            }

            static const ParamMap &GetParameters(void)
            {
                return m_Parameters;
            }

            const std::string &GetSolverInfo(const std::string &lhs);
            bool SolverInfoExists(const std::string &property);

            const std::string &GetFunction(const std::string &lhs);
            Equation GetFunctionAsEquation(const std::string &lhs);

            /// Will look for the lhs equal to str and if found
            /// will return the function in str and return true.
            /// If not found it will return false and leave str
            /// as it was coming in.
            bool SubstituteFunction(std::string &str);

        protected:
            void ReadSolverInfo(TiXmlElement *functions);
            void ReadParameters(TiXmlElement *parameters);
            void ReadVariables(TiXmlElement *variables);
            void ReadFunctions(TiXmlElement *conditions);
            void ReadBoundaryRegions(TiXmlElement *regions);
            void ReadBoundaryConditions(TiXmlElement *conditions);
            void ReadForcingFunctions(TiXmlElement *functions);
            void ReadInitialConditions(TiXmlElement *conditions);
            void ReadExactSolution(TiXmlElement *solution);
            void ReadUserDefinedEqn(TiXmlElement *functions);

            // Containers to hold conditions and associated data
            static ParamMap m_Parameters;
            FunctionMap     m_Functions;
            Variable        m_Variables;
            BoundaryRegionCollection    m_BoundaryRegions;
            BoundaryConditionCollection m_BoundaryConditions;
            ForcingFunctionsMap         m_ForcingFunctions;
            InitialConditionsMap        m_InitialConditions;
            ExactSolutionMap            m_ExactSolution;
            UserDefinedEqnMap           m_UserDefinedEqn;

            SolverInfoMap               m_SolverInfo; //< Solver Information 

            /// The mesh graph to use for referencing geometry info.
            const MeshGraph *m_MeshGraph;

        private:
            BoundaryConditions();
        };
        
        typedef boost::shared_ptr<BoundaryConditions> BoundaryConditionsSharedPtr;
    }
}

#endif //NEKTAR_SPATIALDOMAINS_BOUNDARYCONDITIONS_H
    

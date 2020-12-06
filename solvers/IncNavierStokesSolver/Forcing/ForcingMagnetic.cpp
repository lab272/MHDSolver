///////////////////////////////////////////////////////////////////////////////
//
// File: ForcingMagnetic.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA), Alexander Proskurin (Rus).
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
// Description: Programmatic forcing
//
///////////////////////////////////////////////////////////////////////////////

#include <boost/core/ignore_unused.hpp>

#include <IncNavierStokesSolver/Forcing/ForcingMagnetic.h>

namespace Nektar
{
namespace SolverUtils
{

    std::string ForcingMagnetic::className = GetForcingFactory().
                                RegisterCreatorFunction("Magnetic",
                                                        ForcingMagnetic::create,
                                                        "Magnetic Forcing");

    ForcingMagnetic::ForcingMagnetic(
            const LibUtilities::SessionReaderSharedPtr &pSession,
            const std::weak_ptr<EquationSystem>      &pEquation)
        : Forcing(pSession, pEquation)
    {
    }

    Array<OneD, Array<OneD, NekDouble> >& ForcingMagnetic::UpdateForces()
    {
        return m_Forcing;
    }

    void ForcingMagnetic::v_InitObject(
            const Array<OneD, MultiRegions::ExpListSharedPtr>& pFields,
            const unsigned int& pNumForcingFields,
            const TiXmlElement* pForce)
    {
        int i;
        
        boost::ignore_unused(pForce);
                      
        m_adjoint = false;

        m_NumVariable = m_equ.lock()->GetSpaceDim();
        if(m_NumVariable == 1)
        {
            ASSERTL0(false,"ForcingMagnetic does not support 1D calculations.");
        }
        int phystot = pFields[0]->GetTotPoints();
        
        m_session->LoadParameter("Cndctvt", m_conductivity, 0.0);
 
//         Check MAGNETICFIELD function existence.
        const TiXmlElement* magfieldElmt = pForce->FirstChildElement("MAGNETICFIELD");
        ASSERTL0(magfieldElmt,
                "Requires MAGNETICFIELD tag, specifying the magnetic field function. "
                "the magnetic field function.");
        
        m_magneticFieldFuncName = magfieldElmt->GetText();
        ASSERTL0(m_session->DefinesFunction(m_magneticFieldFuncName),
                 "Function '" + m_magneticFieldFuncName + "' not defined.");
        
//         Check the magnetic field vector elements. They should be in general: "Hx", "Hy", "Hz".
//         If the flow completely 2D and use only "u" and "v" fields, "Hz" may be skipped.
//         In this case the "phi" field also may be skipped and the equation 
//         \Delta \varphi = \nabla(\boldsymbol{v} \times \boldsymbol{H}) is not solved.
        int magdim = 0;
        for (int i = 0; i < 3; ++i)
        {
            std::string s_FieldStr = kMagFieldStr[i]; 
            if(m_session->DefinesFunction(m_magneticFieldFuncName, s_FieldStr))
            {
                magdim++;
            }
            else if(s_FieldStr != kMagFieldStr[2])
            {
                ASSERTL0(false, "Variable '" + s_FieldStr + "' not defined.");
            }
        }
        
        if(m_NumVariable == 2 && magdim == 2)
        {
            m_2D = true;
        }
        else
        {
            m_2D = false;
        }

        m_magneticField = Array<OneD, Array<OneD, NekDouble> > (magdim);
        for (int i = 0; i < magdim; ++i)
        {
//             std::string s_FieldStr = kMagFieldStr[i];
            m_magneticField[i] = Array<OneD, NekDouble> (phystot, 0.0);
//             GetFunction(pFields, m_session, m_magneticFieldFuncName)->Evaluate(s_FieldStr, m_magneticField[i]);
        }
        updateMagneticField(pFields,0.0);
        
        if(!m_2D){
            int numfields = m_equ.lock()->UpdateFields().num_elements();
            bool findphi = false;
            for(i = 0; i < numfields; ++i)
            {
                std::string varName = m_equ.lock()->GetVariable(i);
                if(varName == "phi")
                {
                    m_phi = m_equ.lock()->UpdateFields()[i];
                    findphi = true;
                }
            }
            ASSERTL0(findphi,"Variable 'phi' must be defined.");
        }
        
        m_Forcing = Array<OneD, Array<OneD, NekDouble> > (m_NumVariable);
        for (int i = 0; i < m_NumVariable; ++i)
        {
            m_Forcing[i] = Array<OneD, NekDouble> (phystot, 0.0);
        }
        
        Array<OneD, Array<OneD, NekDouble> > tempfields(pFields.num_elements());
        for (int i = 0; i < m_NumVariable; ++i)
        {
            tempfields[i] = pFields[i]->GetPhys();
        }
        
        Update(pFields,tempfields,0.0);
        
    }
    
    void ForcingMagnetic::updateMagneticField(
            const Array< OneD, MultiRegions::ExpListSharedPtr > &pFields,
            const NekDouble &time)
    {
        for (int i = 0; i < m_magneticField.num_elements(); ++i)
        {
            std::string s_FieldStr = kMagFieldStr[i]; 
            GetFunction(pFields, m_session, m_magneticFieldFuncName)->Evaluate(s_FieldStr, m_magneticField[i],time);
        }
    }

    void ForcingMagnetic::v_Apply(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> > &inarray,
            Array<OneD, Array<OneD, NekDouble> > &outarray,
            const NekDouble &time)
    {
//         boost::ignore_unused(fields, inarray, time);
        
        Update(fields, inarray, time);

        for (int i = 0; i < m_NumVariable; i++)
        {
            Vmath::Vadd(outarray[i].num_elements(), outarray[i], 1,
                        m_Forcing[i], 1, outarray[i], 1);
        }
    }
    
    void ForcingMagnetic::Update(
            const Array< OneD, MultiRegions::ExpListSharedPtr > &pFields,
            const Array<OneD, Array<OneD, NekDouble> > &inarray,
            const NekDouble &time)
    {
        int phystot = pFields[0]->GetTotPoints();
        int i,j;
        
        updateMagneticField(pFields,time);
        
        if(m_2D)
        {
//             Calculates $\sigma(\boldsymbol{v} \times \boldsymbol{B}_0)\times \boldsymbol{B}_0$
//             in the case $\boldsymbol{B}_0 = (B_x,B_y,0)$ and $\boldsymbol{v} = (u,v,0)$
            Array<OneD, NekDouble> tempvalue = Array<OneD, NekDouble>(phystot, 0.0);
//             $F_x = - B_{y} \left(- B_{x} v + B_{y} u\right)$
            Vmath::Vmul(phystot,m_magneticField[0],1,inarray[1],1,tempvalue,1);
            Vmath::Vvtvm(phystot,m_magneticField[1],1,inarray[0],1,tempvalue,1,tempvalue,1);
            Vmath::Vmul(phystot,m_magneticField[1],1,tempvalue,1,tempvalue,1);
            Vmath::Smul(phystot,-1.0*m_conductivity,tempvalue,1,m_Forcing[0],1);
//             $F_y = B_{x} \left(- B_{x} v + B_{y} u\right)$
            Vmath::Vmul(phystot,m_magneticField[0],1,inarray[1],1,tempvalue,1);
            Vmath::Vvtvm(phystot,m_magneticField[1],1,inarray[0],1,tempvalue,1,tempvalue,1);
            Vmath::Vmul(phystot,m_magneticField[0],1,tempvalue,1,tempvalue,1);
            Vmath::Smul(phystot,m_conductivity,tempvalue,1,m_Forcing[1],1);

            return;
        }
        
        updatePhi(pFields,inarray,time);
        
        int dim = 3;
        Array<OneD, Array<OneD, NekDouble> > gradphi(dim),tempvalue(dim),tempvalue2;
        for(int i = 0; i < dim; i++)
        {
            gradphi[i] = Array<OneD, NekDouble>(phystot, 0.0);
            tempvalue[i] = Array<OneD, NekDouble>(phystot, 0.0);
        }
        
        VectorProduct(inarray, m_magneticField,tempvalue);
        
        NekDouble gradsign = -1.0;
        if(m_adjoint)
        {
            gradsign = 1.0;
        }
        for(int i = 0; i < dim; i++)
        {
            m_phi->PhysDeriv(i, m_phi->GetPhys(), gradphi[i]);
            Vmath::Svtvp(phystot,gradsign,gradphi[i],1,tempvalue[i],1,tempvalue[i],1);
        }
        
        tempvalue2 = gradphi;
        
        VectorProduct(tempvalue, m_magneticField,tempvalue2);
        
        for(int i = 0; i < dim; i++)
        {
            Vmath::Smul(phystot,m_conductivity,tempvalue[i],1,m_Forcing[i],1);
        }
        
        return;
        
    }
    
    
    void ForcingMagnetic::updatePhi(
            const Array< OneD, MultiRegions::ExpListSharedPtr > &pFields,
            const Array<OneD, Array<OneD, NekDouble> > &inarray,
            const NekDouble &time)
    {
        StdRegions::ConstFactorMap factors;
        // Setup coefficient for equation
        factors[StdRegions::eFactorLambda] = 0.0;

        int phystot = pFields[0]->GetTotPoints();
        int dim = 3;

        Array<OneD, Array<OneD, NekDouble> > productVB(dim);
        for(int i = 0; i < dim; i++){
           productVB[i] = Array<OneD, NekDouble>(phystot, 0.0);
        }

//      Calculate $\nabla \cdot \left( \boldsymbol{v} \times \boldsymbol{B}_0 \right)$
        VectorProduct(inarray, m_magneticField, productVB);
        Array<OneD, NekDouble> rhs(phystot, 0.0), div_i(phystot, 0.0);
        for(int i = 0; i < dim; i++){
           pFields[0]->PhysDeriv(i, productVB[i], div_i);
           Vmath::Vadd(phystot, rhs, 1, div_i, 1, rhs, 1);
        }
        
        if(m_adjoint)
        {
            Vmath::Neg(phystot, rhs, 1);
        }
                
        m_phi->HelmSolve(rhs, m_phi->UpdateCoeffs(), NullFlagList, factors);
        m_phi->BwdTrans(m_phi->GetCoeffs(), m_phi->UpdatePhys());
    }

// Calculate the vector product of inarray1 and inarray2  
    void ForcingMagnetic::VectorProduct(
        const Array<OneD,const Array<OneD, NekDouble> > &inarray1,
        const Array<OneD,const Array<OneD, NekDouble> > &inarray2,
        Array<OneD, Array<OneD, NekDouble> > &outarray)
    {

        unsigned int np = inarray1[0].num_elements();
        Array<OneD, NekDouble> temp(np, 0.0);

        Vmath::Vmul(np, inarray1[2], 1, inarray2[1], 1, temp, 1);
        Vmath::Vvtvm(np, inarray1[1], 1, inarray2[2], 1, temp, 1, outarray[0], 1);

        Vmath::Vmul(np, inarray1[0], 1, inarray2[2], 1, temp, 1);
        Vmath::Vvtvm(np, inarray1[2], 1, inarray2[0], 1, temp, 1, outarray[1], 1);

        Vmath::Vmul(np, inarray1[1], 1, inarray2[0], 1, temp, 1);
        Vmath::Vvtvm(np, inarray1[0], 1, inarray2[1], 1, temp, 1, outarray[2], 1);

    }
}
}

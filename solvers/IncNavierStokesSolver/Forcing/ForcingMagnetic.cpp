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
                      
        std::vector<std::string> fieldname{"Hx", "Hy", "Hz"};

        m_NumVariable = pNumForcingFields;
        int nq         = pFields[0]->GetTotPoints();
        
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
        
        m_session->LoadParameter("Cndctvt", m_conductivity, 0.0);
 
        const TiXmlElement* magfieldElmt = pForce->FirstChildElement("MAGNETICFIELD");
        ASSERTL0(magfieldElmt,
                "Requires MAGNETICFIELD tag, specifying the magnetic field function. "
                "the magnetic field function.");
        
        m_funcName = magfieldElmt->GetText();
        ASSERTL0(m_session->DefinesFunction(m_funcName),
                 "Function '" + m_funcName + "' not defined.");

        m_magneticField = Array<OneD, Array<OneD, NekDouble> > (fieldname.size());
        for (int i = 0; i < fieldname.size(); ++i)
        {
            std::string s_FieldStr = fieldname.at(i); 
            ASSERTL0(m_session->DefinesFunction(m_funcName, s_FieldStr),
                        "Variable '" + s_FieldStr + "' not defined.");
            m_magneticField[i] = Array<OneD, NekDouble> (nq, 0.0);
            GetFunction(pFields, m_session, m_funcName)->Evaluate(s_FieldStr, m_magneticField[i]);
        }
        
        m_Forcing = Array<OneD, Array<OneD, NekDouble> > (m_NumVariable);
        for (int i = 0; i < m_NumVariable; ++i)
        {
            m_Forcing[i] = Array<OneD, NekDouble> (nq, 0.0);
        }
        

        
    }

    void ForcingMagnetic::v_Apply(
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> > &inarray,
            Array<OneD, Array<OneD, NekDouble> > &outarray,
            const NekDouble &time)
    {
        boost::ignore_unused(fields, inarray, time);

        for (int i = 0; i < m_NumVariable; i++)
        {
            Vmath::Vadd(outarray[i].num_elements(), outarray[i], 1,
                        m_Forcing[i], 1, outarray[i], 1);
        }
    }

// Calculate the vector product of inarray1 and inarray2  
    void ForcingMagnetic::VectorProduct(
        const Array<OneD,
        const Array<OneD, NekDouble> > &inarray1,
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

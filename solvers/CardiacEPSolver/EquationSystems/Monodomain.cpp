///////////////////////////////////////////////////////////////////////////////
//
// File Monodomain.cpp
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
// Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
// Department of Aeronautics, Imperial College London (UK), and Scientific
// Computing and Imaging Institute, University of Utah (USA).
//
// License for the specific language governing rights and limitations under
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
// Description: Monodomain cardiac electrophysiology homogenised model.
//
///////////////////////////////////////////////////////////////////////////////

#include <iostream>

#include <CardiacEPSolver/EquationSystems/Monodomain.h>

namespace Nektar
{
    /**
     * @class Monodomain
     *
     * Base model of cardiac electrophysiology of the form
     * \f{align*}{
     *     \frac{\partial u}{\partial t} = \nabla^2 u + J_{ion},
     * \f}
     * where the reaction term, \f$J_{ion}\f$ is defined by a specific cell
     * model.
     *
     * This implementation, at present, treats the reaction terms explicitly
     * and the diffusive element implicitly.
     */

    /**
     * Registers the class with the Factory.
     */
    string Monodomain::className
            = GetEquationSystemFactory().RegisterCreatorFunction(
                "Monodomain",
                Monodomain::create,
                "Monodomain model of cardiac electrophysiology.");


    /**
     *
     */
    Monodomain::Monodomain(
            const LibUtilities::SessionReaderSharedPtr& pSession)
        : UnsteadySystem(pSession)
    {
    }

    void Monodomain::v_InitObject()
    {
        UnsteadySystem::v_InitObject();

        m_session->LoadParameter("Chi",        m_chi);
        m_session->LoadParameter("Cm",         m_capMembrane);

        std::string vCellModel;
        m_session->LoadSolverInfo("CELLMODEL", vCellModel, "");

        ASSERTL0(vCellModel != "", "Cell Model not specified.");

        m_cell = GetCellModelFactory().CreateInstance(vCellModel, m_session, m_fields[0]);

        m_intVariables.push_back(0);

        // Load variable coefficients
        /// @todo: move to an EvaluateFunction routine.
        std::string varCoeffs[3] = {"d00", "d11", "d22"};
        StdRegions::VarCoeffType varCoeffEnum[3] = {
                StdRegions::eVarCoeffD00,
                StdRegions::eVarCoeffD11,
                StdRegions::eVarCoeffD22
        };
        std::string varName = "intensity";

        for (int i = 0; i < 3; ++i)
        {
            if (m_session->DefinesFunction(varCoeffs[i]))
            {
                // Load from FLD file.
                if (m_session->GetFunctionType(varCoeffs[i]) == LibUtilities::eFunctionTypeFile)
                {
                    ImportFld(m_session->GetFunctionFilename(varCoeffs[i]),
                              m_fields[0],
                              varName);

                    m_fields[0]->BwdTrans(m_fields[0]->GetCoeffs(),
                                          m_fields[0]->UpdatePhys());

                    // Normalise and invert (assuming image intensity data)
                    int nq = m_fields[0]->GetNpoints();
                    NekDouble f_min = m_session->GetParameter("d_min");
                    NekDouble f_max = m_session->GetParameter("d_max");
                    NekDouble f_range = f_max - f_min;
                    NekDouble o_min = m_session->GetParameter("o_min");
                    NekDouble o_max = m_session->GetParameter("o_max");
                    Vmath::Sadd(nq, -f_min, m_fields[0]->GetPhys(), 1, m_fields[0]->UpdatePhys(), 1);
                    for (int j = 0; j < nq; ++j)
                    {
                        if (m_fields[0]->GetPhys()[j] < 0)
                        {
                            m_fields[0]->UpdatePhys()[j] = 0.0;
                        }
                        if (m_fields[0]->GetPhys()[j] > f_range)
                        {
                            m_fields[0]->UpdatePhys()[j] = f_range;
                        }
                    }
                    Vmath::Smul(nq, -1.0/f_range, m_fields[0]->GetPhys(), 1, m_fields[0]->UpdatePhys(), 1);
                    Vmath::Sadd(nq, 1.0, m_fields[0]->GetPhys(), 1, m_fields[0]->UpdatePhys(), 1);
                    Vmath::Smul(nq, o_max-o_min, m_fields[0]->GetPhys(), 1, m_fields[0]->UpdatePhys(), 1);
                    Vmath::Sadd(nq, o_min, m_fields[0]->GetPhys(), 1, m_fields[0]->UpdatePhys(), 1);

                    Array<OneD, NekDouble> tmp(nq);
                    Vmath::Vcopy(nq, m_fields[0]->GetPhys(), 1, tmp, 1);
                    m_vardiff[varCoeffEnum[i]] = tmp;
                }
                // Evaluate expression
                else
                {
                    int nq = m_fields[0]->GetNpoints();
                    Array<OneD,NekDouble> x0(nq);
                    Array<OneD,NekDouble> x1(nq);
                    Array<OneD,NekDouble> x2(nq);

                    // get the coordinates
                    m_fields[0]->GetCoords(x0,x1,x2);

                    Array<OneD, NekDouble> tmp(nq);

                    LibUtilities::EquationSharedPtr ifunc
                            = m_session->GetFunction(varCoeffs[i], varName);
                    ifunc->Evaluate(x0,x1,x2,0.0,tmp);
                    m_vardiff[varCoeffEnum[i]] = tmp;
                }

                // Dump actual variable coefficients for verification.
                m_fields[0]->FwdTrans_IterPerExp(m_fields[0]->GetPhys(),
                                                 m_fields[0]->UpdateCoeffs());
                std::stringstream filename;
                filename << varCoeffs[i];
                if (m_comm->GetSize() > 1)
                {
                    filename << "_P" << m_comm->GetRank();
                }
                filename << ".fld";
                WriteFld(filename.str());
            }
        }

        if (m_session->DefinesParameter("StimulusDuration"))
        {
            ASSERTL0(m_session->DefinesFunction("Stimulus", "u"),
                    "Stimulus function not defined.");
            m_session->LoadParameter("StimulusDuration", m_stimDuration);
        }
        else
        {
            m_stimDuration = 0;
        }

        if (!m_explicitDiffusion)
        {
            m_ode.DefineImplicitSolve (&Monodomain::DoImplicitSolve, this);
        }
        m_ode.DefineOdeRhs(&Monodomain::DoOdeRhs, this);
    }


    /**
     *
     */
    Monodomain::~Monodomain()
    {

    }


    /**
     * @param   inarray         Input array.
     * @param   outarray        Output array.
     * @param   time            Current simulation time.
     * @param   lambda          Timestep.
     */
    void Monodomain::DoImplicitSolve(
            const Array<OneD, const Array<OneD, NekDouble> >&inarray,
                  Array<OneD, Array<OneD, NekDouble> >&outarray,
            const NekDouble time,
            const NekDouble lambda)
    {
        int nvariables  = inarray.num_elements();
        int ncoeffs     = inarray[0].num_elements();
        int nq          = m_fields[0]->GetNpoints();
        StdRegions::ConstFactorMap factors;
        // lambda = \Delta t
        factors[StdRegions::eFactorLambda] = 1.0/lambda*m_chi*m_capMembrane;

        // We solve ( \nabla^2 - HHlambda ) Y[i] = rhs [i]
        // inarray = input: \hat{rhs} -> output: \hat{Y}
        // outarray = output: nabla^2 \hat{Y}
        // where \hat = modal coeffs
        for (int i = 0; i < nvariables; ++i)
        {
            // Multiply 1.0/timestep
            Vmath::Smul(nq, -factors[StdRegions::eFactorLambda], inarray[i], 1,
                                            m_fields[i]->UpdatePhys(), 1);

            // Solve a system of equations with Helmholtz solver and transform
            // back into physical space.
            m_fields[i]->HelmSolve(m_fields[i]->GetPhys(),
                                   m_fields[i]->UpdateCoeffs(), NullFlagList,
                                   factors, m_vardiff);

            m_fields[i]->BwdTrans( m_fields[i]->GetCoeffs(),
                                   m_fields[i]->UpdatePhys());
            m_fields[i]->SetPhysState(true);

            // Copy the solution vector (required as m_fields must be set).
            outarray[i] = m_fields[i]->GetPhys();
        }
    }


    void Monodomain::DoOdeRhs(
            const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                  Array<OneD,        Array<OneD, NekDouble> >&outarray,
            const NekDouble time)
    {
        int nq = m_fields[0]->GetNpoints();
        int nvar = inarray.num_elements();

        m_cell->TimeIntegrate(inarray, outarray, time);

        if (m_stimDuration > 0 && time < m_stimDuration)
        {
            Array<OneD,NekDouble> x0(nq);
            Array<OneD,NekDouble> x1(nq);
            Array<OneD,NekDouble> x2(nq);
            Array<OneD,NekDouble> result(nq);

            // get the coordinates
            m_fields[0]->GetCoords(x0,x1,x2);

            LibUtilities::EquationSharedPtr ifunc
                    = m_session->GetFunction("Stimulus", "u");
            ifunc->Evaluate(x0,x1,x2,time, result);

            Vmath::Vadd(nq, outarray[0], 1, result, 1, outarray[0], 1);
        }
        Vmath::Smul(nq, 1.0/m_capMembrane, outarray[0], 1, outarray[0], 1);
    }


    void Monodomain::v_SetInitialConditions(NekDouble initialtime,
                        bool dumpInitialConditions)
    {
        EquationSystem::v_SetInitialConditions(initialtime, dumpInitialConditions);
        m_cell->Initialise();
    }


    /**
     *
     */
    void Monodomain::v_PrintSummary(std::ostream &out)
    {
        UnsteadySystem::v_PrintSummary(out);
        if (m_session->DefinesFunction("d00") &&
            m_session->GetFunctionType("d00") == LibUtilities::eFunctionTypeExpression)
        {
            out << "\tDiffusivity-x   : "
                << m_session->GetFunction("d00", "intensity")->GetExpression()
                << endl;
        }
        if (m_session->DefinesFunction("d11") &&
            m_session->GetFunctionType("d11") == LibUtilities::eFunctionTypeExpression)
        {
            out << "\tDiffusivity-x   : "
                << m_session->GetFunction("d11", "intensity")->GetExpression()
                << endl;
        }
        if (m_session->DefinesFunction("d22") &&
            m_session->GetFunctionType("d22") == LibUtilities::eFunctionTypeExpression)
        {
            out << "\tDiffusivity-x   : "
                << m_session->GetFunction("d22", "intensity")->GetExpression()
                << endl;
        }
        m_cell->PrintSummary(out);
    }

}

///////////////////////////////////////////////////////////////////////////////
//
// File: VanDerWaalsEoS.cpp
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
// Description: Van der Waals equation of state
//
///////////////////////////////////////////////////////////////////////////////

#include "VanDerWaalsEoS.h"

using namespace std;

namespace Nektar
{

std::string VanDerWaalsEoS::className = GetEquationOfStateFactory().
    RegisterCreatorFunction("VanDerWaals",
                            VanDerWaalsEoS::create,
                            "Van der Waals equation of state.");

VanDerWaalsEoS::VanDerWaalsEoS(const LibUtilities::SessionReaderSharedPtr& pSession)
    : EquationOfState(pSession)
{
    NekDouble Tcrit, Pcrit;
    pSession->LoadParameter ("Tcrit", Tcrit);
    pSession->LoadParameter ("Pcrit", Pcrit);

    m_a = 27.0/64.0 * m_gasConstant * m_gasConstant * Tcrit * Tcrit / Pcrit;
    m_b = 1.0/8.0 * Tcrit / Pcrit;
}

NekDouble VanDerWaalsEoS::v_GetTemperature(
    const NekDouble &rho, const NekDouble &e)
{
    return (e + m_a * rho) * (m_gamma - 1) / m_gasConstant;
}

NekDouble VanDerWaalsEoS::v_GetPressure(
    const NekDouble &rho, const NekDouble &e)
{
    return (e + m_a * rho) * (m_gamma - 1) / (1.0/rho - m_b) - m_a * rho * rho;
}

NekDouble VanDerWaalsEoS::v_GetEntropy(
    const NekDouble &rho, const NekDouble &e)
{
    NekDouble T   = GetTemperature(rho,e);
    NekDouble sIg =
            m_gasConstant/(m_gamma-1) * log(T) - m_gasConstant * log(rho);

    return sIg + m_gasConstant * log(1-m_b*rho);
}

NekDouble VanDerWaalsEoS::v_GetDPDrho_e(
    const NekDouble &rho, const NekDouble &e)
{
    NekDouble result;

    result = (m_gamma - 1) * (e + 2*m_a*rho - m_a*m_b*rho*rho);
    result = result / ( (1 - m_b*rho) * (1 - m_b*rho) );
    result = result - 2*m_a*rho;

    return result;
}

NekDouble VanDerWaalsEoS::v_GetDPDe_rho(
    const NekDouble &rho, const NekDouble &e)
{
    return (m_gamma - 1) / (1.0/rho - m_b);
}

NekDouble VanDerWaalsEoS::v_GetEFromRhoP(
            const NekDouble &rho, const NekDouble &p)
{
    return (p + m_a*rho*rho) * (1.0/rho - m_b) / (m_gamma - 1) - m_a * rho;
}

}

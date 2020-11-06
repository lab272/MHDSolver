///////////////////////////////////////////////////////////////////////////////
//
// File VelocityCorrectionScheme.h
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
// Description: Velocity Correction Scheme header
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_VELOCITYCORRECTIONSCHEME_H
#define NEKTAR_SOLVERS_VELOCITYCORRECTIONSCHEME_H

#include <WeakMHDSolver/EquationSystems/IncNavierStokes.h>

namespace Nektar
{
    
enum NonLinDisturbance  {
    NonLinDisturbanceNone,
    NonLinDisturbanceRandom,
    NonLinDisturbanceSize
}; 

const std::string kNonLinDisturbanceStr[] = {
    "None",
    "Random"
};
    
class VelocityCorrectionScheme: public IncNavierStokes
{
public:

    /// Creates an instance of this class
    static SolverUtils::EquationSystemSharedPtr create (
        const LibUtilities::SessionReaderSharedPtr& pSession ) {
        SolverUtils::EquationSystemSharedPtr p =
            MemoryManager<VelocityCorrectionScheme>::
            AllocateSharedPtr ( pSession );
        p->InitObject();
        return p;
    }

    /// Name of class
    static std::string className;


    /// Constructor.
    VelocityCorrectionScheme ( const LibUtilities::SessionReaderSharedPtr& pSession );

    virtual ~VelocityCorrectionScheme();

    virtual void v_InitObject();
    
        void v_SetInitialConditions(
                NekDouble initialtime = 0.0,
                bool dumpInitialConditions = true,
                const int domain = 0);

    void SetUpPressureForcing (
        const Array<OneD, const Array<OneD, NekDouble> > &fields,
        Array<OneD, Array<OneD, NekDouble> > &Forcing,
        const NekDouble aii_Dt ) {
        v_SetUpPressureForcing ( fields, Forcing, aii_Dt );
    }

    void SolveElectricPotential ();

    void SpecialHVecMult ( const Array<OneD, const Array<OneD, NekDouble> > &inarray1, const Array<OneD, const Array<OneD, NekDouble> > &inarray2, Array<OneD, Array<OneD, NekDouble> > &outarray );

    void SetUpMagneticForcing( const Nektar::Array< Nektar::OneD, const Nektar::Array< Nektar::OneD, Nektar::NekDouble > >& inarray, Nektar::Array< Nektar::OneD, Nektar::Array< Nektar::OneD, Nektar::NekDouble > >& outarray );

    
    
    //     void SetUpMagneticForcing (
//         const Nektar::Array< Nektar::OneD, const Nektar::Array< Nektar::OneD, Nektar::NekDouble > >& inarray, Nektar::Array< Nektar::OneD, Nektar::Array< Nektar::OneD, Nektar::NekDouble > >& Forcing );

//     void SetUpMagneticForcingNS (
//         const Nektar::Array< Nektar::OneD, const Nektar::Array< Nektar::OneD, Nektar::NekDouble > >& inarray, Nektar::Array< Nektar::OneD, Nektar::Array< Nektar::OneD, Nektar::NekDouble > >& Forcing );

    void SetUpViscousForcing (
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
        Array<OneD, Array<OneD, NekDouble> > &Forcing,
        const NekDouble aii_Dt ) {
        v_SetUpViscousForcing ( inarray, Forcing, aii_Dt );
    }

    void SolvePressure ( const Array<OneD, NekDouble>  &Forcing ) {
        v_SolvePressure ( Forcing );
    }

    void SolveViscous (
        const Array<OneD, const Array<OneD, NekDouble> > &Forcing,
        Array<OneD, Array<OneD, NekDouble> > &outarray,
        const NekDouble aii_Dt ) {
        v_SolveViscous ( Forcing, outarray, aii_Dt );
    }

    void SolveUnsteadyStokesSystem (
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
        Array<OneD, Array<OneD, NekDouble> > &outarray,
        const NekDouble time,
        const NekDouble a_iixDt );

    void EvaluateAdvection_SetPressureBCs (
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
        Array<OneD, Array<OneD, NekDouble> > &outarray,
        const NekDouble time ) {
        v_EvaluateAdvection_SetPressureBCs ( inarray, outarray, time );
    }

protected:
    /// Set up initial conditions for nonlinear disturbance.
    NonLinDisturbance m_nonlinDist; 
    ///Initial nonlinear disturbance amplotude.
    NekDouble m_initialDisturbanceAmplitude;
    /// bool to set up $\phi$ field to zero in 2D case.
    bool m_PhiZero;
    /// bool to set up masgnetic force to zero.
    bool m_MagForceZero;
    /// bool to identify if spectral vanishing viscosity is active.
    bool m_useHomo1DSpecVanVisc;
    /// bool to identify if spectral vanishing viscosity is active.
    bool m_useSpecVanVisc;
    /// cutt off ratio from which to start decayhing modes
    NekDouble m_sVVCutoffRatio;
    /// Diffusion coefficient of SVV modes
    NekDouble m_sVVDiffCoeff;

    /// Save aiiDt value to use as a trip to reset global matrix setup
    Array<OneD, NekDouble> m_saved_aii_Dt;

    /// Variable Coefficient map for the Laplacian which can be activated as part of SVV or otherwise
    StdRegions::VarCoeffMap m_varCoeffLap;

    // Virtual functions
    virtual void v_GenerateSummary ( SolverUtils::SummaryList& s );

    virtual void v_TransCoeffToPhys ( void );

    virtual void v_TransPhysToCoeff ( void );

    virtual void v_DoInitialise ( void );

    virtual Array<OneD, bool> v_GetSystemSingularChecks();

    virtual int v_GetForceDimension();

    virtual void v_SetUpPressureForcing (
        const Array<OneD, const Array<OneD, NekDouble> > &fields,
        Array<OneD, Array<OneD, NekDouble> > &Forcing,
        const NekDouble aii_Dt );

    virtual void v_SetUpViscousForcing (
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
        Array<OneD, Array<OneD, NekDouble> > &Forcing,
        const NekDouble aii_Dt );

    virtual void v_SolvePressure ( const Array<OneD, NekDouble>  &Forcing );

    virtual void v_SolveViscous (
        const Array<OneD, const Array<OneD, NekDouble> > &Forcing,
        Array<OneD, Array<OneD, NekDouble> > &outarray,
        const NekDouble aii_Dt );

    virtual void v_EvaluateAdvection_SetPressureBCs (
        const Array<OneD, const Array<OneD, NekDouble> > &inarray,
        Array<OneD, Array<OneD, NekDouble> > &outarray,
        const NekDouble time );

private:

};

typedef boost::shared_ptr<VelocityCorrectionScheme>
VelocityCorrectionSchemeSharedPtr;

} //end of namespace


#endif //VELOCITY_CORRECTION_SCHEME_H

///////////////////////////////////////////////////////////////////////////////
//
// File Stimulus.h
//
// For more information, please see: http://www.nektar.info
//
// The MIT License
//
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
// Description: Stimulus class.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SOLVERS_CARDIACEPSOLVER_STIMULI_STIMULUSRECT
#define NEKTAR_SOLVERS_CARDIACEPSOLVER_STIMULI_STIMULUSRECT

#include <LibUtilities/BasicUtils/NekFactory.hpp>
#include <LibUtilities/BasicUtils/SessionReader.h>
#include <LibUtilities/BasicUtils/SharedArray.hpp>
//#include <SpatialDomains/SpatialData.h>
#include <MultiRegions/ExpList.h>
#include <StdRegions/StdNodalTriExp.h>
#include <StdRegions/StdNodalTetExp.h>
#include <CardiacEPSolver/Stimuli/Stimulus.h>

namespace Nektar
{
    // Forward declaration
    class StimulusRect;
    
    
    
    /// Protocol base class.
    class StimulusRect: public Stimulus
    {
    public:
        Stimulus(const LibUtilities::SessionReaderSharedPtr& pSession,
                 const MultiRegions::ExpListSharedPtr& pField);
        
        virtual ~Stimulus() {}
        
        /// Initialise the cell model storage and set initial conditions
        void Initialise();
        
    protected:
        NekDouble m_px1;
        NekDouble m_py1;
        NekDouble m_pz1;
        NekDouble m_px2;
        NekDouble m_py2;
        NekDouble m_pz2;
        
        
        
        virtual void v_Update(
                              const Array<OneD, const  Array<OneD, NekDouble> >&inarray,
                              Array<OneD,        Array<OneD, NekDouble> >&outarray,
                              const NekDouble time) = 0;
        
        virtual void v_PrintSummary(std::ostream &out) = 0;
        
    };
    
}

#endif 

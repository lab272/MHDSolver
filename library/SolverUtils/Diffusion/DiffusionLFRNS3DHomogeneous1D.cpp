///////////////////////////////////////////////////////////////////////////////
//
// File: DiffusionLFRNS3DHomogeneous1D.cpp
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
///////////////////////////////////////////////////////////////////////////////
//
// File: DiffusionLFR.cpp
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
// Description: LFR diffusion 3DHomogeneous1D class.
//
///////////////////////////////////////////////////////////////////////////////

#include <SolverUtils/Diffusion/DiffusionLFRNS3DHomogeneous1D.h>
#include <LibUtilities/Polylib/Polylib.h>
#include <boost/math/special_functions/gamma.hpp>
#include <iostream>
#include <iomanip>

namespace Nektar
{
    namespace SolverUtils
    {
        std::string DiffusionLFRNS3DHomogeneous1D::type[] = {
            GetDiffusionFactory().RegisterCreatorFunction(
                "LFRDGNS3DHomogeneous1D", DiffusionLFRNS3DHomogeneous1D::create),
            GetDiffusionFactory().RegisterCreatorFunction(
                "LFRSDNS3DHomogeneous1D", DiffusionLFRNS3DHomogeneous1D::create),
            GetDiffusionFactory().RegisterCreatorFunction(
               "LFRHUNS3DHomogeneous1D", DiffusionLFRNS3DHomogeneous1D::create),
            GetDiffusionFactory().RegisterCreatorFunction(
                "LFRcminNS3DHomogeneous1D", DiffusionLFRNS3DHomogeneous1D::create),
            GetDiffusionFactory().RegisterCreatorFunction(
                "LFRcinfNS3DHomogeneous1D", DiffusionLFRNS3DHomogeneous1D::create)};
        
        /**
         * @brief DiffusionLFRNS3DHomogeneous1D uses the 2D Flux Reconstruction (FR) 
         * approach to compute the diffusion term looping on the planes in the z
         * direction and adding the flux in z direction at the end.
         * The implementation is only for segments,
         * quadrilaterals and hexahedra at the moment.
         *
         * \todo Extension to triangles, tetrahedra and other shapes.
         * (Long term objective)
         */
        DiffusionLFRNS3DHomogeneous1D::DiffusionLFRNS3DHomogeneous1D
                                    (std::string diffType):m_diffType(diffType)
        {
            string diffName;
            
            if (diffType == "LFRDGNS3DHomogeneous1D")
            {
                diffName = "LFRDGNS";
            }
            else if (diffType == "LFRSDNS3DHomogeneous1D")
            {
                diffName = "LFRSDNS";
            }
            else if (diffType == "LFRHUNS3DHomogeneous1D")
            {
                diffName = "LFRHUNS";
            }
            else if (diffType == "LFRcminNS3DHomogeneous1D")
            {
                diffName = "LFRcminNS";
            }
            else if (diffType == "LFRcinfNS3DHomogeneous1D")
            {
                diffName = "LFRcinfNS";
            }
            
            m_planeDiff = GetDiffusionFactory().CreateInstance(diffName, diffName);
        }
        
        /**
         * @brief Initiliase DiffusionLFRNS3DHomogeneous1D objects and store them 
         * before starting the time-stepping.
         *
         * This routine calls the virtual functions #v_SetupMetrics,
         * #v_SetupCFunctions and #v_SetupInterpolationMatrices to
         * initialise the objects needed by DiffusionLFR.
         *
         * @param pSession  Pointer to session reader.
         * @param pFields   Pointer to fields.
         */
        void DiffusionLFRNS3DHomogeneous1D::v_InitObject(
            LibUtilities::SessionReaderSharedPtr        pSession,
            Array<OneD, MultiRegions::ExpListSharedPtr> pFields)
        {
            int nConvectiveFields = pFields.num_elements();
            nVariables = nConvectiveFields - 1;
            
            Array<OneD, MultiRegions::ExpListSharedPtr>
            pFields_plane0(nConvectiveFields);
            
            for (int i = 0; i < nConvectiveFields; ++i)
            {
                pFields_plane0[i] = pFields[i]->GetPlane(0);
            }
            
            m_planeDiff->InitObject(pSession, pFields_plane0);
            
            spaceDim = 3;
            nPointsTot      = pFields[0]->GetTotPoints();
            nCoeffs         = pFields[0]->GetNcoeffs();
            
            planes = pFields[0]->GetZIDs();
            num_planes = planes.num_elements();
            
            nPointsTot_plane = nPointsTot/num_planes;
            nCoeffs_plane = nCoeffs/num_planes;
            
            diffLFRNS =
                boost::dynamic_pointer_cast<DiffusionLFRNS>(m_planeDiff);
            
            fluxvector =
            Array<OneD, Array<OneD, NekDouble> > (nConvectiveFields);
            for (j = 0; j < nConvectiveFields; j ++)
            {
                fluxvector[j] = Array<OneD, NekDouble>(nPointsTot);
            }
            
            derivatives01_homo =
            Array<OneD, Array<OneD, NekDouble> > (nVariables);
            for (i = 0; i < nVariables; ++i)
            {
                derivatives01_homo[i] = Array<OneD, NekDouble>(nPointsTot, 0.0);
            }
            
            fields_plane =
                    Array <OneD, Array<OneD, MultiRegions::ExpListSharedPtr> >
                                                                  (num_planes);
            inarray_plane =
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (num_planes);
            outarray_plane =
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (num_planes);
            derivatives01_plane =
                Array<OneD, Array<OneD, Array<OneD, NekDouble> > > (num_planes);
            fluxvector_homo =
                Array<OneD, Array<OneD, Array<OneD, Array<OneD, NekDouble> > > >
                                                                   (num_planes);
            
            for (i = 0; i < num_planes; ++i)
            {
                fields_plane[i] = Array<OneD, MultiRegions::ExpListSharedPtr>
                                                            (nConvectiveFields);
                inarray_plane[i] = Array<OneD, Array<OneD, NekDouble> >
                                                            (nVariables);
                derivatives01_plane[i] = Array<OneD, Array<OneD, NekDouble> >
                                                            (nVariables);
                outarray_plane[i] = Array<OneD, Array<OneD, NekDouble> >
                                                            (nConvectiveFields);
                fluxvector_homo[i] =
                    Array<OneD, Array<OneD, Array<OneD, NekDouble> > >(spaceDim);
                
                for (j = 0; j < nConvectiveFields; j ++)
                {
                    outarray_plane[i][j] = Array<OneD, NekDouble>
                    (nPointsTot_plane, 0.0);
                }
                
                for (j = 0; j < nVariables; j ++)
                {
                    inarray_plane[i][j] = Array<OneD, NekDouble>
                                                        (nPointsTot_plane, 0.0);
                    derivatives01_plane[i][j] = Array<OneD, NekDouble>
                                                        (nPointsTot_plane, 0.0);
                }
                
                for (j = 0; j < spaceDim ; ++j)
                {
                    fluxvector_homo[i][j] = Array<OneD, Array<OneD, NekDouble> >
                                                            (nConvectiveFields);
                    
                    for (int k = 0; k < nConvectiveFields ; ++k)
                    {
                        fluxvector_homo[i][j][k] =
                            Array<OneD, NekDouble>(nPointsTot_plane, 0.0);
                    }
                }
            }
            
            outarray_homo = Array<OneD, Array<OneD, NekDouble> >
                                                        (nConvectiveFields);
            for (j = 0; j < nConvectiveFields; j ++)
            {
                outarray_homo[j] = Array<OneD, NekDouble>(nPointsTot, 0.0);
            }
            
            m_planeDiff->SetRiemannSolver(m_riemann);
            m_planeDiff->SetFluxVectorVecNS(m_fluxVectorNS);
        }
        
        /**
         * @brief Calculate FR Diffusion for the Navier-Stokes (NS) equations
         * using an LDG interface flux and the the flux in the third direction.
         *
         * The equations that need a diffusion operator are those related
         * with the velocities and with the energy.
         */
        void DiffusionLFRNS3DHomogeneous1D::v_Diffuse(
            const int                                         nConvectiveFields,
            const Array<OneD, MultiRegions::ExpListSharedPtr> &fields,
            const Array<OneD, Array<OneD, NekDouble> >        &inarray,
                  Array<OneD, Array<OneD, NekDouble> >        &outarray)
        {
            for (i = 0; i < nVariables; ++i)
            {
                fields[0]->PhysDeriv(2, inarray[i], derivatives01_homo[i]);
            }
            
            for (i = 0; i < num_planes; ++i)
            {
                for (j = 0; j < nConvectiveFields; j ++)
                {
                    fields_plane[i][j]= fields[j]->GetPlane(i);
                }
                
                for (j = 0; j < nVariables; j ++)
                {
                    Vmath::Vcopy(nPointsTot_plane,
                                 &inarray[j][i * nPointsTot_plane], 1,
                                 &inarray_plane[i][j][0], 1);
                    
                    Vmath::Vcopy(nPointsTot_plane,
                                 &derivatives01_homo[j][i * nPointsTot_plane], 1,
                                 &derivatives01_plane[i][j][0], 1);
                }
                
                // Set the first order derivatives in 3rd direction
                diffLFRNS->SetHomoDerivs(derivatives01_plane[i]);
                
                m_planeDiff->Diffuse(nConvectiveFields,
                                     fields_plane[i],
                                     inarray_plane[i],
                                     outarray_plane[i]);
                
                for (j = 0; j < nConvectiveFields; j ++)
                {
                    Vmath::Vcopy(nPointsTot_plane,
                                 &outarray_plane[i][j][0], 1,
                                 &outarray[j][i * nPointsTot_plane], 1);
                }
                
                m_planeDiff->FluxVec(fluxvector_homo[i]);
                
                for ( j = 0; j < nConvectiveFields; ++j)
                {
                    Vmath::Vcopy(nPointsTot_plane,
                                 &fluxvector_homo[i][2][j][0], 1,
                                 &fluxvector[j][i * nPointsTot_plane], 1);
                }
            }
            
            for (i = 0; i < nConvectiveFields; ++i)
            {
                fields[0]->PhysDeriv(2, fluxvector[i], outarray_homo[i]);
                
                Vmath::Vadd(nPointsTot,
                            outarray[i], 1,
                            outarray_homo[i], 1,
                            outarray[i], 1);
            }
        }
    }// close namespace SolverUtils
}// close namespace nektar++

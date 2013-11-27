///////////////////////////////////////////////////////////////////////////////
//
// File: GeomFactors3D.cpp
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
// Description: Implementation of 3D geometric factors.
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_SPATIALDOMAINS_GEOMFACTORS3D_H
#define NEKTAR_SPATIALDOMAINS_GEOMFACTORS3D_H

#include <SpatialDomains/GeomFactors.h>
#include <SpatialDomains/SpatialDomainsDeclspec.h>
#include <StdRegions/StdExpansion3D.h>

namespace Nektar
{
    namespace SpatialDomains
    {
        /// Geometric factors relating to the coordinate mapping
        /// \f$\chi:\Omega_{st}\rightarrow \Omega_e \f$ for three-dimensional
        /// elements
        class GeomFactors3D : public GeomFactors
        {
        public:
            /// Constructor.
            SPATIAL_DOMAINS_EXPORT GeomFactors3D(
                const GeomType gtype,
                const int coordim,
                const Array<OneD, const StdRegions
                    ::StdExpansion3DSharedPtr> &Coords,
                const Array<OneD, const LibUtilities::BasisSharedPtr>
                    &tbasis);
            
        private:
            // Set up 3D geometric factors
            void SetUpJacGmat3D(
                const Array<OneD, Array<OneD, NekDouble> > d1,
                const Array<OneD, Array<OneD, NekDouble> > d2,
                const Array<OneD, Array<OneD, NekDouble> > d3);

            /// Tests if the element is valid and not self-intersecting.
            void CheckIfValid();

            /// Performs 3D interpolation between two sets of point
            /// distributions.
            virtual void v_Interp(
                        const PointsKeyArray &map_points,
                        const Array<OneD, const NekDouble> &src,
                        const PointsKeyArray &tpoints,
                        Array<OneD, NekDouble> &tgt) const;

            virtual void v_Adjoint(
                        const Array<TwoD, const NekDouble>& src,
                        Array<TwoD, NekDouble>& tgt) const;

        };

        /// Shared pointer to GeomFactors3D object.
        typedef boost::shared_ptr<GeomFactors3D>      GeomFactors3DSharedPtr;
        /// Vector of shared pointers to GeomFactors3D objects.
        typedef std::vector< GeomFactors3DSharedPtr > GeomFactors3DVector;
        /// Iterator for the vector of shared pointers.
        typedef GeomFactors3DVector::iterator GeomFactors3DVectorIter;
    }
}

#endif

///////////////////////////////////////////////////////////////////////////////
//
// File $Source: /usr/sci/projects/Nektar/cvs/Nektar++/library/LocalRegions/PyrExp.h,v $
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
// Description: Header file for PyrExp routines
//
///////////////////////////////////////////////////////////////////////////////

#ifndef PYREXP_H
#define PYREXP_H

#include <LocalRegions/LocalRegions.hpp>

#include <StdRegions/StdPyrExp.h>
#include <SpatialDomains/PyrGeom.h>

#include <SpatialDomains/GeomFactors.h>
#include <LocalRegions/MatrixKey.h>




namespace Nektar
{
  namespace LocalRegions
  {

    class PyrExp: public StdRegions::StdPyrExp
    {

    public:

      /** \brief Constructor using BasisKey class for quadrature points and order
      definition */
        PyrExp(const LibUtilities::BasisKey &Ba,
               const LibUtilities::BasisKey &Bb,
	       const LibUtilities::BasisKey &Bc,
               const SpatialDomains::PyrGeomSharedPtr &geom);

        PyrExp(const LibUtilities::BasisKey &Ba,
	       const LibUtilities::BasisKey &Bb,
	       const LibUtilities::BasisKey &Bc);

	    
      /// Copy Constructor
      PyrExp(const PyrExp &T);

      /// Destructor
      ~PyrExp();

       void GetCoords(Array<OneD,NekDouble> &coords_0,
		       Array<OneD,NekDouble> &coords_1,
          	       Array<OneD,NekDouble> &coords_2);

	void GetCoord(const ConstArray<OneD,NekDouble> &Lcoords, 
                      Array<OneD,NekDouble> &coords);

	 //----------------------------
        // Integration Methods
        //----------------------------

        /// \brief Integrate the physical point list \a inarray over region
       NekDouble Integral(const ConstArray<OneD,NekDouble> &inarray);    

       void IProductWRTBase(const ConstArray<OneD,NekDouble> &inarray, Array<OneD,NekDouble> &outarray);
       void FwdTrans(const ConstArray<OneD,NekDouble> & inarray,Array<OneD,NekDouble> &outarray);

       NekDouble PhysEvaluate(const ConstArray<OneD,NekDouble> &coord);

      //-----------------------------
      // Differentiation Methods
      //-----------------------------
        void PhysDeriv(const ConstArray<OneD, NekDouble> &inarray, 
                       Array<OneD, NekDouble> &out_d0,
                       Array<OneD, NekDouble> &out_d1,
                       Array<OneD, NekDouble> &out_d2);

      /// Return Shape of region, using  ShapeType enum list. i.e. Pyramid
      StdRegions::ShapeType DetShapeType() const
      {
    	return StdRegions::ePyramid;
      }

       SpatialDomains::PyrGeomSharedPtr GetGeom()
      {
          return m_geom;
      }

      void WriteToFile(FILE *outfile);
      void WriteToFile(std::ofstream &outfile, const int dumpVar);


    protected:
      //SpatialDomains::PyrGeom * m_geom;
	void GenMetricInfo();

	void IProductWRTBase(const ConstArray<OneD,NekDouble> &base0, 
			     const ConstArray<OneD,NekDouble> &base1, 
			     const ConstArray<OneD,NekDouble> &base2, 
			     const ConstArray<OneD,NekDouble> &inarray,
			     Array<OneD,NekDouble> &outarray);

        DNekMatSharedPtr CreateStdMatrix(const StdRegions::StdMatrixKey &mkey);
        DNekScalMatSharedPtr  CreateMatrix(const MatrixKey &mkey);
        DNekScalBlkMatSharedPtr  CreateStaticCondMatrix(const MatrixKey &mkey);

	SpatialDomains::PyrGeomSharedPtr m_geom;
        SpatialDomains::GeomFactorsSharedPtr  m_metricinfo;

	LibUtilities::NekManager<MatrixKey, DNekScalMat, MatrixKey::opLess> m_matrixManager;
        LibUtilities::NekManager<MatrixKey, DNekScalBlkMat, MatrixKey::opLess> m_staticCondMatrixManager;


    private:
      PyrExp();

      virtual StdRegions::ShapeType v_DetShapeType() const
      {
        return DetShapeType();
      }

      virtual SpatialDomains::GeomFactorsSharedPtr v_GetMetricInfo() const
        {
            return m_metricinfo;
        }

         virtual void v_GetCoords(Array<OneD, NekDouble> &coords_0,
                                 Array<OneD, NekDouble> &coords_1 = NullNekDouble1DArray,
                                 Array<OneD, NekDouble> &coords_2 = NullNekDouble1DArray)
         {
             GetCoords(coords_0, coords_1, coords_2);
         }

        virtual void v_GetCoord(const ConstArray<OneD, NekDouble> &lcoord, 
                                Array<OneD, NekDouble> &coord)
        {
            GetCoord(lcoord, coord);
        }

        virtual int v_GetCoordim()
        {
            return m_geom->GetCoordim();
        }

        /// Virtual call to SegExp::PhysDeriv
        virtual void v_StdPhysDeriv(const ConstArray<OneD, NekDouble> &inarray, 
                                    Array<OneD, NekDouble> &out_d0,
                                    Array<OneD, NekDouble> &out_d1,
                                    Array<OneD, NekDouble> &out_d2)
        {
            StdPyrExp::PhysDeriv(inarray, out_d0, out_d1, out_d2);
        }

        virtual void v_PhysDeriv(const ConstArray<OneD, NekDouble> &inarray,
                                 Array<OneD, NekDouble> &out_d0,
                                 Array<OneD, NekDouble> &out_d1,
                                 Array<OneD, NekDouble> &out_d2)
        {
            PhysDeriv(inarray, out_d0, out_d1, out_d2);
        }


        virtual void v_WriteToFile(FILE *outfile)
        {
            WriteToFile(outfile);
        }
        
        virtual void v_WriteToFile(std::ofstream &outfile, const int dumpVar)
        {
            WriteToFile(outfile,dumpVar);
        }

        /** \brief Virtual call to integrate the physical point list \a inarray
        over region (see SegExp::Integral) */
        virtual NekDouble v_Integral(const ConstArray<OneD, NekDouble> &inarray )
        {
            return Integral(inarray);
        }

        /** \brief Virtual call to TriExp::IProduct_WRT_B */
        virtual void v_IProductWRTBase(const ConstArray<OneD, NekDouble> &inarray,
                                       Array<OneD, NekDouble> &outarray)
        {
            IProductWRTBase(inarray,outarray);
        }

        virtual void v_IProductWRTDerivBase (const int dir,
                                             const ConstArray<OneD,NekDouble> &inarray,
                                             Array<OneD, NekDouble> &outarray)
        {
            IProductWRTDerivBase(dir,inarray,outarray);
        }

        /// Virtual call to SegExp::FwdTrans
        virtual void v_FwdTrans(const ConstArray<OneD, NekDouble> &inarray,
                                Array<OneD, NekDouble> &outarray)
        {
            FwdTrans(inarray, outarray);
        }
        
        /// Virtual call to PyrExp::Evaluate
        virtual NekDouble v_PhysEvaluate(const ConstArray<OneD, NekDouble> &coords)
        {
            return PhysEvaluate(coords);
        }

        virtual NekDouble v_Linf(const ConstArray<OneD, NekDouble> &sol)
        {
            return Linf(sol);
        }

        virtual NekDouble v_Linf()
        {
            return Linf();
        }

        virtual NekDouble v_L2(const ConstArray<OneD, NekDouble> &sol)
        {
            return StdExpansion::L2(sol);
        }

        virtual NekDouble v_L2()
        {
            return StdExpansion::L2();
        }

        virtual DNekMatSharedPtr v_CreateStdMatrix(const StdRegions::StdMatrixKey &mkey)
        {
            return CreateStdMatrix(mkey);
        }

        virtual DNekScalMatSharedPtr& v_GetLocMatrix(const MatrixKey &mkey)
        {
            return m_matrixManager[mkey];
        }

        virtual DNekScalBlkMatSharedPtr& v_GetLocStaticCondMatrix(const MatrixKey &mkey)
        {
            return m_staticCondMatrixManager[mkey];
        }




    };

    // type defines for use of PyrExp in a boost vector
    typedef boost::shared_ptr<PyrExp> PyrExpSharedPtr;
    typedef std::vector< PyrExpSharedPtr > PyrExpVector;
    typedef std::vector< PyrExpSharedPtr >::iterator PyrExpVectorIter;

  } //end of namespace
} //end of namespace

#endif  //PYREXP_H

/**
 *    $Log: PyrExp.h,v $
 *    Revision 1.11  2008/03/17 10:35:15  pvos
 *    Clean up of the code
 *
 *    Revision 1.10  2008/03/12 15:24:29  pvos
 *    Clean up of the code
 *
 *    Revision 1.9  2008/02/16 05:51:36  ehan
 *    Added PhysDeriv and virtual functions.
 *
 *    Revision 1.8  2008/02/05 00:40:19  ehan
 *    Added initial pyramidic expansion.
 *
 *    Revision 1.7  2007/07/22 23:04:18  bnelson
 *    Backed out Nektar::ptr.
 *
 *    Revision 1.6  2007/07/20 00:45:51  bnelson
 *    Replaced boost::shared_ptr with Nektar::ptr
 *
 *    Revision 1.5  2007/01/15 21:12:26  sherwin
 *    First definition
 *
 *    Revision 1.4  2006/06/13 18:05:01  sherwin
 *    Modifications to make MultiRegions demo ProjectLoc2D execute properly.
 *
 *    Revision 1.3  2006/06/05 00:09:22  bnelson
 *    Fixed a gcc 4.1.0 compilation warning with the #endif at the end of the file.
 *
 *    Revision 1.2  2006/05/30 14:00:03  sherwin
 *    Updates to make MultiRegions and its Demos work
 *
 *    Revision 1.1  2006/05/04 18:58:46  kirby
 *    *** empty log message ***
 *
 *    Revision 1.9  2006/03/12 21:59:48  sherwin
 *
 *    compiling version of LocalRegions
 *
 *    Revision 1.8  2006/03/12 07:43:32  sherwin
 *
 *    First revision to meet coding standard. Needs to be compiled
 *
 **/

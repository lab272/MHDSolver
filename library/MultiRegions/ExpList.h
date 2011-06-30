///////////////////////////////////////////////////////////////////////////////
//
// File ExpList.h
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
// Description: Expansion list top class definition
//
///////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_LIBS_MULTIREGIONS_EXPLIST_H
#define NEKTAR_LIBS_MULTIREGIONS_EXPLIST_H
#include <MultiRegions/MultiRegionsDeclspec.h>
#include <LibUtilities/Communication/Comm.h>
#include <MultiRegions/MultiRegions.hpp>
#include <StdRegions/StdExpansion.h>
#include <MultiRegions/LocalToGlobalBaseMap.h>
#include <MultiRegions/GlobalMatrix.h>
#include <MultiRegions/GlobalMatrixKey.h>
#include <MultiRegions/GlobalLinSysKey.h>

#include <LocalRegions/MatrixKey.h>
#include <SpatialDomains/SegGeom.h>
#include <LocalRegions/PointExp.h>

#include <SpatialDomains/MeshGraph.h>
#include <SpatialDomains/MeshGraph2D.h>
#include <MultiRegions/GlobalOptimizationParameters.h>
#include <boost/enable_shared_from_this.hpp>

namespace Nektar
{
    namespace MultiRegions
    {
        // Forward declarations
        class GlobalLinSys;
        class LocalToGlobalC0ContMap;
        class LocalToGlobalBaseMap;
        class LocalToGlobalDGMap;
        class ExpList0D;
        class ExpList1D;
	class ExpList2D;
	
	enum Direction
	{
	    eX,
	    eY,
	    eZ,
	    eS
	};
		

        /// A map between global matrix keys and their associated block
        /// matrices.
        typedef map<GlobalMatrixKey,DNekScalBlkMatSharedPtr> BlockMatrixMap;
        /// A shared pointer to a BlockMatrixMap.
        typedef boost::shared_ptr<BlockMatrixMap> BlockMatrixMapShPtr;
		
		enum FourierSpaceType
		{
			ePhys,
			eCoef,
			eNotDef
		};
		
        /// Base class for all multi-elemental spectral/hp expansions.
        class ExpList: public boost::enable_shared_from_this<ExpList>
        {
        public:
            /// The default constructor.
            MULTI_REGIONS_EXPORT ExpList();

            /// The default constructor.
            ExpList(LibUtilities::CommSharedPtr &pComm);

            /// The copy constructor.
            MULTI_REGIONS_EXPORT ExpList(const ExpList &in, bool DeclareCoeffPhysArrays = true);

            /// The default destructor.
            MULTI_REGIONS_EXPORT virtual ~ExpList();

            /// Copy coefficients from concatenated list to expansion list.
            MULTI_REGIONS_EXPORT void PutCoeffsInToElmtExp(void);

            /// Copy coefficients from expansion list to concatenated list.
            MULTI_REGIONS_EXPORT void PutElmtExpInToCoeffs(void);

            /// Copy one elements coefficients from the concatenated list
            /// to the expansion list.
            MULTI_REGIONS_EXPORT void PutCoeffsInToElmtExp(int eid);

            /// Copy one elements coefficients from the expansion list to
            /// the concatenated list.
            MULTI_REGIONS_EXPORT void PutElmtExpInToCoeffs(int eid);

            /// Copy physical data from \a m_phys to expansion list.
            MULTI_REGIONS_EXPORT void PutPhysInToElmtExp(void);

            /// Copy physical data from given array to expansion list.
            MULTI_REGIONS_EXPORT void PutPhysInToElmtExp(Array<OneD, const NekDouble> &in);

            /// Copy expansion list physical data to given array.
            MULTI_REGIONS_EXPORT void PutElmtExpInToPhys(Array<OneD,NekDouble> &out);

            /// Copy expansion list physical data from one element to array.
            MULTI_REGIONS_EXPORT void PutElmtExpInToPhys(int eid, Array<OneD,NekDouble> &out);

            /// Returns the total number of local degrees of freedom
            /// \f$N_{\mathrm{eof}}=\sum_{e=1}^{{N_{\mathrm{el}}}}N^{e}_m\f$.
            inline int GetNcoeffs(void) const;

            // Returns the total number of local degrees of freedom
            // for element eid
            MULTI_REGIONS_EXPORT int GetNcoeffs(const int eid) const;
			
			inline int GetContNcoeffs(void) const;

            /// Evaulates the maximum number of modes in the elemental basis
            /// order over all elements
            inline int EvalBasisNumModesMax(void) const;

            /// Returns the vector of the number of modes in the elemental
            /// basis order over all elements.
            MULTI_REGIONS_EXPORT const Array<OneD,int> EvalBasisNumModesMaxPerExp(void) const;

            /// Returns the total number of quadrature points #m_npoints
            /// \f$=Q_{\mathrm{tot}}\f$.
            inline int GetTotPoints(void) const;

            /// Returns the total number of quadrature points for eid's element
            /// \f$=Q_{\mathrm{tot}}\f$.
            inline int GetTotPoints(const int eid) const;

            /// Returns the total number of quadrature points #m_npoints
            /// \f$=Q_{\mathrm{tot}}\f$.
            inline int GetNpoints(void) const;

            /// Sets the transformed state #m_transState of the coefficient
            /// arrays.
            inline void SetTransState(const TransState transState);

            /// This function returns the transformed state #m_transState of
            /// the coefficient arrays.
            inline TransState GetTransState(void) const;
			
			/// Sets the Fourier Space to the one of the possible configuration
			/// ePhys, eCoef or eNotSet
            inline void SetFourierSpace(const FourierSpaceType fourierspace);
			
            /// This function returns the Fourie Space condition, which is stored
			/// in the variable m_FourierSpace.
            inline FourierSpaceType GetFourierSpace(void) const;

            /// Fills the array #m_phys
            inline void SetPhys(const Array<OneD, const NekDouble> &inarray);

            /// Sets the array #m_phys
            inline void SetPhysArray(Array<OneD, NekDouble> &inarray);

            /// This function manually sets whether the array of physical
            /// values \f$\boldsymbol{u}_l\f$ (implemented as #m_phys) is
            /// filled or not.
            inline void SetPhysState(const bool physState);

            /// This function indicates whether the array of physical values
            /// \f$\boldsymbol{u}_l\f$ (implemented as #m_phys) is filled or
            /// not.
            inline bool GetPhysState(void) const;

            /// This function integrates a function \f$f(\boldsymbol{x})\f$
            /// over the domain consisting of all the elements of the expansion.
            MULTI_REGIONS_EXPORT NekDouble PhysIntegral (void);

            /// This function integrates a function \f$f(\boldsymbol{x})\f$
            /// over the domain consisting of all the elements of the expansion.
            MULTI_REGIONS_EXPORT NekDouble PhysIntegral(const Array<OneD, const NekDouble> &inarray);

            /// This function calculates the inner product of a function
            /// \f$f(\boldsymbol{x})\f$ with respect to all \emph{local}
            /// expansion modes \f$\phi_n^e(\boldsymbol{x})\f$.
            MULTI_REGIONS_EXPORT void   IProductWRTBase_IterPerExp(
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,       NekDouble> &outarray);

            ///
            inline void IProductWRTBase(
                                const Array<OneD, const NekDouble> &inarray,
                                Array<OneD,       NekDouble> &outarray,
                                bool  UseContCoeffs = false);

            /// This function calculates the inner product of a function
            /// \f$f(\boldsymbol{x})\f$ with respect to the derivative (in
            /// direction \param dir) of all \emph{local} expansion modes
            /// \f$\phi_n^e(\boldsymbol{x})\f$.
            MULTI_REGIONS_EXPORT void   IProductWRTDerivBase(const int dir,
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,       NekDouble> &outarray);

            /// This function elementally evaluates the forward transformation
            /// of a function \f$u(\boldsymbol{x})\f$ onto the global
            /// spectral/hp expansion.
            MULTI_REGIONS_EXPORT void   FwdTrans_IterPerExp (
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,       NekDouble> &outarray);

            ///
            inline void FwdTrans(
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,       NekDouble> &outarray,
                                bool  UseContCoeffs = false);

            /// This function elementally mulplies the coefficient space of
            /// Sin my the elemental inverse of the mass matrix.
            MULTI_REGIONS_EXPORT void  MultiplyByElmtInvMass (
                                const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,       NekDouble> &outarray);

            ///
            inline void MultiplyByInvMassMatrix(
                                const Array<OneD,const NekDouble> &inarray,
                                      Array<OneD,      NekDouble> &outarray,
                                bool  UseContCoeffs = false);

            /// Solve helmholtz problem
            inline void HelmSolve(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                          NekDouble lambda = 0.0,
                    const Array<OneD, const NekDouble> &varLambda
                                                = NullNekDouble1DArray,
                    const Array<OneD, const Array<OneD, NekDouble> > &varCoeff
                                                = NullNekDoubleArrayofArray);

            /// Solve helmholtz problem (continuous case parameters).
            inline void HelmSolve(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                          NekDouble lambda,
                          bool UseContCoeffs,
                    const Array<OneD, const NekDouble> &dirForcing
                                                = NullNekDouble1DArray,
                    const Array<OneD, const NekDouble> &varLambda
                                                = NullNekDouble1DArray,
                    const Array<OneD, const Array<OneD, NekDouble> > &varCoeff
                                                = NullNekDoubleArrayofArray);

            /// Solve helmholtz problem (discontinuous case parameters).
            inline void HelmSolve(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray,
                          NekDouble lambda,
                          NekDouble tau,
                    const Array<OneD, const NekDouble> &varLambda
                                                = NullNekDouble1DArray,
                    const Array<OneD, const Array<OneD, NekDouble> > &varCoeff
                                                = NullNekDoubleArrayofArray);

            /// Solve Advection Diffusion Reaction
            inline void LinearAdvectionDiffusionReactionSolve(
                       const Array<OneD, Array<OneD, NekDouble> > &velocity,
                       const Array<OneD, const NekDouble> &inarray,
                       Array<OneD, NekDouble> &outarray,
                       const NekDouble lambda,
                       bool  UseContCoeffs = false,
                       const Array<OneD, const NekDouble>&
                       dirForcing = NullNekDouble1DArray);


            /// Solve Advection Diffusion Reaction
            inline void LinearAdvectionReactionSolve(
                       const Array<OneD, Array<OneD, NekDouble> > &velocity,
                       const Array<OneD, const NekDouble> &inarray,
                       Array<OneD, NekDouble> &outarray,
                       const NekDouble lambda,
                       bool  UseContCoeffs = false,
                       const Array<OneD, const NekDouble>&
                       dirForcing = NullNekDouble1DArray);

            ///
            MULTI_REGIONS_EXPORT void FwdTrans_BndConstrained(
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray);


            /// This function elementally evaluates the backward transformation
            /// of the global spectral/hp element expansion.
            MULTI_REGIONS_EXPORT void BwdTrans_IterPerExp (
                    const Array<OneD, const NekDouble> &inarray,
                          Array<OneD,       NekDouble> &outarray);

            ///
            inline void BwdTrans (const Array<OneD, const NekDouble> &inarray,
                                  Array<OneD,NekDouble> &outarray,
                                  bool  UseContCoeffs = false);

            /// This function calculates the coordinates of all the elemental
            /// quadrature points \f$\boldsymbol{x}_i\f$.
            inline void GetCoords(Array<OneD, NekDouble> &coord_0,
                                  Array<OneD, NekDouble> &coord_1 = NullNekDouble1DArray,
                                  Array<OneD, NekDouble> &coord_2 = NullNekDouble1DArray);
			
			/// This function calculates the coordinates of all the elemental
            /// quadrature points \f$\boldsymbol{x}_i\f$.
            inline void GetCoords(NekDouble &x, NekDouble &y, NekDouble &z);
			
			/// This function calculates the coordinates of all the elemental
            /// quadrature points \f$\boldsymbol{x}_i\f$.
            inline void GetCoord(Array<OneD, NekDouble> &coords);
			
			// Homogeneous transforms
			inline void HomogeneousFwdTrans(const Array<OneD, const NekDouble> &inarray, 
											  Array<OneD, NekDouble> &outarray, 
											  bool UseContCoeffs = false);
			
			inline void HomogeneousBwdTrans(const Array<OneD, const NekDouble> &inarray, 
											  Array<OneD, NekDouble> &outarray, 
											  bool UseContCoeffs = false);
			
			inline void GetBCValues(Array<OneD, NekDouble> &BndVals, 
									const Array<OneD, NekDouble> &TotField, 
									int BndID);
			
			inline void NormVectorIProductWRTBase(Array<OneD, const NekDouble> &V1,
												  Array<OneD, const NekDouble> &V2,
												  Array<OneD, NekDouble> &outarray,
												  int BndID);
			
            /// This function calculates Surface Normal vector of a smooth
            /// manifold.
            MULTI_REGIONS_EXPORT void GetSurfaceNormal(Array<OneD,NekDouble> &SurfaceNormal,
                                  const int k);

            /// Populate tangents vector with tangents from each element.
            MULTI_REGIONS_EXPORT void GetTangents(
                             Array<OneD, Array<OneD, Array<OneD, NekDouble> > > &tangents);

            /// Apply geometry information to each expansion.
            MULTI_REGIONS_EXPORT void ApplyGeomInfo(SpatialDomains::MeshGraph &graph);

            /// This function writes the spectral/hp element solution to the
            /// file \a out.
            MULTI_REGIONS_EXPORT void WriteToFile(std::ofstream &out,
                             OutputFormat format = eTecplot,
                             std::string var = "v");

            void WriteTecplotHeader(std::ofstream &outfile,
                                    std::string var = "v")
            {
                v_WriteTecplotHeader(outfile,var);
            }

            void WriteTecplotZone(std::ofstream &outfile, int expansion)
            {
                v_WriteTecplotZone(outfile,expansion);
            }

            void WriteTecplotField(std::ofstream &outfile, int expansion)
            {
                v_WriteTecplotField(outfile,expansion);
            }

            MULTI_REGIONS_EXPORT void WriteVtkHeader(std::ofstream &outfile);
            MULTI_REGIONS_EXPORT void WriteVtkFooter(std::ofstream &outfile);

            void WriteVtkPieceHeader(std::ofstream &outfile, int expansion)
            {
                v_WriteVtkPieceHeader(outfile, expansion);
            }

            MULTI_REGIONS_EXPORT void WriteVtkPieceFooter(std::ofstream &outfile, int expansion);

            void WriteVtkPieceData  (std::ofstream &outfile, int expansion,
                                     std::string var = "v")
            {
                v_WriteVtkPieceData(outfile, expansion, var);
            }

            MULTI_REGIONS_EXPORT void ReadFromFile(std::ifstream &in,
                              OutputFormat format = eTecplot);

            /// This function returns the dimension of the coordinates of the
            /// element \a eid.
            // inline
            MULTI_REGIONS_EXPORT int GetCoordim(int eid);

            /// Set the \a i th coefficiient in \a m_coeffs to value \a val
            inline void SetCoeff(int i, NekDouble val);
			
			/// Set the coefficiient in \a m_coeffs to value \a val (0D Exapnsion)
            inline void SetCoeff(NekDouble val);
			
			/// Set the physical value in \a m_coeffs to value \a val (0D Exapnsion)
            inline void SetPhys(NekDouble val);
			
			inline const SpatialDomains::VertexComponentSharedPtr &GetGeom(void) const;
			
			inline const SpatialDomains::VertexComponentSharedPtr &GetVertex(void) const;

            /// Set the \a i th coefficiient in  #m_coeffs to value \a val
            inline void SetCoeffs(int i, NekDouble val);

            /// Set the  #m_coeffs array to inarray
            inline void SetCoeffsArray(Array<OneD, NekDouble> &inarray);

            /// Set the  #m_contCoeffs array to inarray
            inline void SetContCoeffsArray(Array<OneD, NekDouble> &inarray);

            /// This function returns (a reference to) the array
            /// \f$\boldsymbol{\hat{u}}_l\f$ (implemented as #m_coeffs)
            /// containing all local expansion coefficients.
            inline const Array<OneD, const NekDouble> &GetCoeffs() const;

            /// Returns (a reference to) the array \f$\boldsymbol{\hat{u}}_g\f$
            /// (implemented as #m_contCoeffs) containing all global expansion
            /// coefficients.
            inline const Array<OneD, const NekDouble> &GetContCoeffs() const;

            /// Put the coefficients from m_coeffs into m_contCoeffs using
            /// local to global mapping
            inline void LocalToGlobal(void);

            /// Put the coefficients from m_contCoeffs into m_contCoeffs using
            /// local to global mapping
            inline void GlobalToLocal(void);

            /// Get the \a i th value  (coefficient) of #m_coeffs
            inline NekDouble GetCoeff(int i);

            /// Get the \a i th value  (coefficient) of #m_coeffs
            inline NekDouble GetCoeffs(int i);

            /// This function returns (a reference to) the array
            /// \f$\boldsymbol{u}_l\f$ (implemented as #m_phys) containing the
            /// function \f$u^{\delta}(\boldsymbol{x})\f$ evaluated at the
            /// quadrature points.
            // inline
            MULTI_REGIONS_EXPORT const Array<OneD, const NekDouble> &GetPhys()  const;

            /// This function calculates the \f$L_\infty\f$ error of the global
            /// spectral/hp element approximation.
            MULTI_REGIONS_EXPORT NekDouble Linf (const Array<OneD, const NekDouble> &soln);

            /// This function calculates the \f$L_2\f$ error with
            /// respect to soln of the global
            /// spectral/hp element approximation.
            NekDouble L2 (const Array<OneD, const NekDouble> &soln)
            {
                return v_L2(soln);
            }

            /// This function calculates the \f$L_2\f$ measure of the global
            /// spectral/hp element approximation.
            NekDouble L2 (void)
            {
                return v_L2();
            }

            /// Calculates the \f$H^1\f$ error of the global spectral/hp
            /// element approximation.
            MULTI_REGIONS_EXPORT NekDouble H1 (const Array<OneD, const NekDouble> &soln);

            /// This function returns the number of elements in the expansion.
            inline int GetExpSize(void);

            /// This function returns the number of elements in the
            /// expansion which may be different for a homogeoenous extended
            /// expansionp.
            inline int GetNumElmts(void)
            {
                return v_GetNumElmts();
            }

            /// This function returns the vector of elements in the expansion.
            inline const boost::shared_ptr<StdRegions::StdExpansionVector> GetExp() const;

            /// This function returns (a shared pointer to) the local elemental
            /// expansion of the \f$n^{\mathrm{th}}\f$ element.
            inline StdRegions::StdExpansionSharedPtr& GetExp(int n) const;

            /// This function returns (a shared pointer to) the local elemental
            /// expansion containing the arbitrary point given by \a gloCoord.
            MULTI_REGIONS_EXPORT StdRegions::StdExpansionSharedPtr& GetExp(
                                                      const Array<OneD, const NekDouble> &gloCoord);

            /// This function returns the index of the local elemental
            /// expansion containing the arbitrary point given by \a gloCoord.
            MULTI_REGIONS_EXPORT int GetExpIndex(const Array<OneD, const NekDouble> &gloCoord);

            /// Get the start offset position for a global list of #m_coeffs
            /// correspoinding to element n.
            inline const int GetCoeff_Offset(int n) const;

            /// Get the start offset position for a global list of m_phys
            /// correspoinding to element n.
            inline const int GetPhys_Offset(int n) const;


            /// Get the element id associated with the n th
            /// consecutive block of data in  #m_phys and #m_coeffs
            inline const int GetOffset_Elmt_Id(int n) const;

            /// This function returns (a reference to) the array
            /// \f$\boldsymbol{\hat{u}}_l\f$ (implemented as #m_coeffs)
            /// containing all local expansion coefficients.
            inline Array<OneD, NekDouble> &UpdateCoeffs();

            /// Returns (a reference to) the array \f$\boldsymbol{\hat{u}}_g\f$
            /// (implemented as #m_contCoeffs) containing all global expansion
            /// coefficients in ContField
            inline Array<OneD, NekDouble> &UpdateContCoeffs()
            {
                return v_UpdateContCoeffs();
            }

            /// This function returns (a reference to) the array
            /// \f$\boldsymbol{u}_l\f$ (implemented as #m_phys) containing the
            /// function \f$u^{\delta}(\boldsymbol{x})\f$ evaluated at the
            /// quadrature points.
            inline Array<OneD, NekDouble> &UpdatePhys();

	    inline void PhysDeriv(Direction edir, 
								  const Array<OneD, const NekDouble> &inarray,
								  Array<OneD, NekDouble> &out_d, bool UseContCoeffs = false);	
	    
            /// This function discretely evaluates the derivative of a function
            /// \f$f(\boldsymbol{x})\f$ on the domain consisting of all
            /// elements of the expansion.
			inline void PhysDeriv(const Array<OneD, const NekDouble> &inarray,
								  Array<OneD, NekDouble> &out_d0,
								  Array<OneD, NekDouble> &out_d1 = NullNekDouble1DArray,
								  Array<OneD, NekDouble> &out_d2 = NullNekDouble1DArray, bool UseContCoeffs = false);
			
            inline void PhysDeriv(const int dir,
								  const Array<OneD, const NekDouble> &inarray,
								  Array<OneD, NekDouble> &out_d, bool UseContCoeffs = false);


            // functions associated with DisContField
            inline const Array<OneD, const  boost::shared_ptr<ExpList> > &GetBndCondExpansions();

            inline boost::shared_ptr<ExpList> &UpdateBndCondExpansion(int i);
			
            inline boost::shared_ptr<ExpList1D> &GetTrace();

            inline boost::shared_ptr<LocalToGlobalDGMap> &GetTraceMap(void);

            inline void AddTraceIntegral(
                                         const Array<OneD, const NekDouble> &Fx,
                                         const Array<OneD, const NekDouble> &Fy,
                                         Array<OneD, NekDouble> &outarray);

            inline void AddTraceIntegral(
                                         const Array<OneD, const NekDouble> &Fn,
                                         Array<OneD, NekDouble> &outarray);

            inline void AddTraceBiIntegral(
                                           const Array<OneD, const NekDouble> &Fwd,
                                           const Array<OneD, const NekDouble> &Bwd,
                                           Array<OneD, NekDouble> &outarray);

            inline void GetFwdBwdTracePhys( Array<OneD,NekDouble> &Fwd,
                                            Array<OneD,NekDouble> &Bwd);

            inline void GetFwdBwdTracePhys(
                                           const Array<OneD,const NekDouble> &field,
                                           Array<OneD,NekDouble> &Fwd,
                                           Array<OneD,NekDouble> &Bwd);

            inline void ExtractTracePhys(
                                         Array<OneD,NekDouble> &outarray);

            inline void ExtractTracePhys(
                                         const Array<OneD, const NekDouble> &inarray,
                                         Array<OneD,NekDouble> &outarray);

            inline const Array<OneD, const SpatialDomains
                ::BoundaryConditionShPtr>& GetBndConditions();

            inline Array<OneD, SpatialDomains::BoundaryConditionShPtr>& UpdateBndConditions();

            inline void EvaluateBoundaryConditions(const NekDouble time = 0.0, const NekDouble = NekConstants::kNekUnsetDouble, const NekDouble = NekConstants::kNekUnsetDouble);


            // Routines for continous matrix solution
            /// This function calculates the result of the multiplication of a
            /// matrix of type specified by \a mkey with a vector given by \a
            /// inarray.
            inline void GeneralMatrixOp(const GlobalMatrixKey             &gkey,
                                        const Array<OneD,const NekDouble> &inarray,
                                        Array<OneD,      NekDouble> &outarray,
                                        bool  UseContCoeffs = false);

            void GeneralMatrixOp_IterPerExp( const GlobalMatrixKey      &gkey,
                                             const Array<OneD,const NekDouble> &inarray,
                                             Array<OneD,      NekDouble> &outarray);

            inline void SetUpPhysNormals(const StdRegions::StdExpansionVector &locexp);

 	    inline void SetUpPhysTangents(const StdRegions::StdExpansionVector &locexp);
 	                

            inline void SetUpTangents();

            inline void GetBoundaryToElmtMap(Array<OneD, int> &ElmtID,
                                             Array<OneD,int> &EdgeID);

            MULTI_REGIONS_EXPORT void  GeneralGetFieldDefinitions(std::vector<SpatialDomains::FieldDefinitionsSharedPtr> &fielddef, int NumHomoDir = 0, Array<OneD, LibUtilities::BasisSharedPtr> &HomoBasis = LibUtilities::NullBasisSharedPtr1DArray, std::vector<NekDouble> &HomoLen = SpatialDomains::NullNekDoubleVector);

            /// load global optimisation parameters
            void ReadGlobalOptimizationParameters(
                                                  const std::string &infilename)
            {
                v_ReadGlobalOptimizationParameters(infilename);
            }

            const NekOptimize::GlobalOptParamSharedPtr &GetGlobalOptParam(void)
            {
                return m_globalOptParam;
            }

            map<int, RobinBCInfoSharedPtr> GetRobinBCInfo()
            {
                return v_GetRobinBCInfo();
            }

            void GetPeriodicEdges(SpatialDomains::MeshGraph2D &graph2D,
                                  SpatialDomains::BoundaryConditions &bcs,
                                  const std::string variable,
                                  vector<map<int,int> > & periodicVertices,
                                  map<int,int>& periodicEdges)
            {
                v_GetPeriodicEdges(graph2D,bcs,variable,periodicVertices,periodicEdges);
            }

            std::vector<SpatialDomains::FieldDefinitionsSharedPtr>
                GetFieldDefinitions()
            {
                return v_GetFieldDefinitions();
            }


            void GetFieldDefinitions(std::vector<SpatialDomains::FieldDefinitionsSharedPtr> &fielddef)
            {
                v_GetFieldDefinitions(fielddef);
            }



            /// Append the element data listed in elements
            /// fielddef->m_ElementIDs onto fielddata
            void AppendFieldData(
                                 SpatialDomains::FieldDefinitionsSharedPtr &fielddef,
                                 std::vector<NekDouble> &fielddata)
            {
                v_AppendFieldData(fielddef,fielddata);
            }

            
            /// Append the data in coeffs listed in elements
            /// fielddef->m_ElementIDs onto fielddata
            void AppendFieldData(
                                 SpatialDomains::FieldDefinitionsSharedPtr &fielddef,
                                 std::vector<NekDouble> &fielddata,
                                 Array<OneD, NekDouble> &coeffs)
            {
                v_AppendFieldData(fielddef,fielddata,coeffs);
            }

            /// Extract the data in fielddata into the m_coeff list
            void ExtractDataToCoeffs(
                                     SpatialDomains::FieldDefinitionsSharedPtr &fielddef,
                                     std::vector<NekDouble> &fielddata,
                                     std::string &field)
            {
                v_ExtractDataToCoeffs(fielddef,fielddata,field);
            }


            /// Extract the data in fielddata into the coeffs
            void ExtractDataToCoeffs(
                                     SpatialDomains::FieldDefinitionsSharedPtr &fielddef,
                                     std::vector<NekDouble> &fielddata,
                                     std::string &field,
                                     Array<OneD, NekDouble> &coeffs)
            {
                v_ExtractDataToCoeffs(fielddef,fielddata,field,coeffs);
            }

            /// Returns a shared pointer to the current object.
            boost::shared_ptr<ExpList> GetSharedThisPtr()
            {
                return shared_from_this();
            }

            /// Returns the comm object
            boost::shared_ptr<LibUtilities::Comm> GetComm()
            {
                return m_comm;
            }

            // Wrapper functions for Homogeneous Expansions
            inline LibUtilities::BasisSharedPtr  GetHomogeneousBasis(void)
            {
                return GetHomogeneousBasis();
            }

        protected:
            boost::shared_ptr<DNekMat> GenGlobalMatrixFull(
                                                           const GlobalLinSysKey &mkey,
                                                           const boost::shared_ptr<LocalToGlobalC0ContMap> &locToGloMap);

            /// Communicator
            LibUtilities::CommSharedPtr m_comm;

            /// The total number of local degrees of freedom. #m_ncoeffs
            /// \f$=N_{\mathrm{eof}}=\sum_{e=1}^{{N_{\mathrm{el}}}}N^{e}_l\f$
            int m_ncoeffs;

            /// The total number of quadrature points. #m_npoints
            /// \f$=Q_{\mathrm{tot}}=\sum_{e=1}^{{N_{\mathrm{el}}}}N^{e}_Q\f$
            int m_npoints;

            /**
             * \brief Concatenation of all local expansion coefficients.
             *
             * The array of length #m_ncoeffs\f$=N_{\mathrm{eof}}\f$ which is
             * the concatenation of the local expansion coefficients
             * \f$\hat{u}_n^e\f$ over all \f$N_{\mathrm{el}}\f$ elements
             * \f[\mathrm{\texttt{m\_coeffs}}=\boldsymbol{\hat{u}}_{l} =
             * \underline{\boldsymbol{\hat{u}}}^e = \left [ \begin{array}{c}
             * \boldsymbol{\hat{u}}^{1} \       \
             * \boldsymbol{\hat{u}}^{2} \       \
             * \vdots \                                                 \
             * \boldsymbol{\hat{u}}^{{{N_{\mathrm{el}}}}} \end{array} \right ],
             * \quad
             * \mathrm{where}\quad \boldsymbol{\hat{u}}^{e}[n]=\hat{u}_n^{e}\f]
             */
            Array<OneD, NekDouble> m_coeffs;

            /**
             * \brief The global expansion evaluated at the quadrature points
             *
             * The array of length #m_npoints\f$=Q_{\mathrm{tot}}\f$ containing
             * the evaluation of \f$u^{\delta}(\boldsymbol{x})\f$ at the
             * quadrature points \f$\boldsymbol{x}_i\f$.
             * \f[\mathrm{\texttt{m\_phys}}=\boldsymbol{u}_{l} =
             * \underline{\boldsymbol{u}}^e = \left [ \begin{array}{c}
             * \boldsymbol{u}^{1} \             \
             * \boldsymbol{u}^{2} \             \
             * \vdots \                                                 \
             * \boldsymbol{u}^{{{N_{\mathrm{el}}}}} \end{array} \right ],\quad
             * \mathrm{where}\quad
             * \boldsymbol{u}^{e}[i]=u^{\delta}(\boldsymbol{x}_i)\f]
             */
            Array<OneD, NekDouble> m_phys;

            /**
             * \brief The transformed state of the array of coefficients of the
             * expansion.
             *
             * #TransState is the enumeration which can attain the following
             * values:
             * - <em>eNotSet</em>: The coefficients are not set.
             * - <em>eLocal</em>: The array #m_coeffs is filled with the proper
             *   local coefficients.
             * - <em>eContinuous</em>: The array #m_contCoeffs is filled with
             *   the proper global coefficients.
             * - <em>eLocalCont</em>: Both the arrays #m_coeffs and
             *   #m_contCoeffs are filled with the proper coefficients.
             */
            TransState m_transState;

            /**
             * \brief The state of the array #m_phys.
             *
             * Indicates whether the array #m_phys, created to contain the
             * evaluation of \f$u^{\delta}(\boldsymbol{x})\f$ at the quadrature
             * points \f$\boldsymbol{x}_i\f$, is filled with these values.
             */
            bool       m_physState;

            /**
             * \brief The list of local expansions.
             *
             * The (shared pointer to the) vector containing (shared pointers
             * to) all local expansions. The fact that the local expansions are
             * all stored as a (pointer to a) #StdExpansion, the abstract base
             * class for all local expansions, allows a general implementation
             * where most of the routines for the derived classes are defined
             * in the #ExpList base class.
             */
            boost::shared_ptr<StdRegions::StdExpansionVector> m_exp;

            /// Offset of elemental data into the array #m_coeffs
            Array<OneD, int>  m_coeff_offset;

            /// Offset of elemental data into the array #m_phys
            Array<OneD, int>  m_phys_offset;

            /// Array containing the element id #m_offset_elmt_id[n]
            /// that the n^th consecutive block of data in #m_coeffs
            /// and #m_phys is associated, i.e. for an array of
            /// constant expansion size and single shape elements
            /// m_phys[n*m_npoints] is the data related to
            /// m_exp[m_offset_elmt_id[n]];
            Array<OneD, int>  m_offset_elmt_id;

            NekOptimize::GlobalOptParamSharedPtr m_globalOptParam;

            BlockMatrixMapShPtr  m_blockMat;
			
			FourierSpaceType m_FourierSpace;

            /// This function assembles the block diagonal matrix of local
            /// matrices of the type \a mtype.
            const DNekScalBlkMatSharedPtr GenBlockMatrix(
                                                         const GlobalMatrixKey &gkey);

            const DNekScalBlkMatSharedPtr& GetBlockMatrix(
                                                          const GlobalMatrixKey &gkey);

            void MultiplyByBlockMatrix(
                                       const GlobalMatrixKey             &gkey,
                                       const Array<OneD,const NekDouble> &inarray,
                                       Array<OneD,      NekDouble> &outarray);

            /// Generates a global matrix from the given key and map.
            boost::shared_ptr<GlobalMatrix>  GenGlobalMatrix(
                                                             const GlobalMatrixKey &mkey,
                                                             const boost::shared_ptr<LocalToGlobalC0ContMap> &locToGloMap);


            void GlobalEigenSystem(const boost::shared_ptr<DNekMat> &Gmat,
                                   Array<OneD, NekDouble> &EigValsReal,
                                   Array<OneD, NekDouble> &EigValsImag,
                                   Array<OneD, NekDouble> &EigVecs
                                   = NullNekDouble1DArray);


            /// This operation constructs the global linear system of type \a
            /// mkey.
            boost::shared_ptr<GlobalLinSys>  GenGlobalLinSys(
                                                             const GlobalLinSysKey &mkey,
                                                             const boost::shared_ptr<LocalToGlobalC0ContMap> &locToGloMap);

            /// Generate a GlobalLinSys from information provided by the key
            /// "mkey" and the mapping provided in LocToGloBaseMap.
            boost::shared_ptr<GlobalLinSys> GenGlobalBndLinSys(
                                                               const GlobalLinSysKey     &mkey,
                                                               const LocalToGlobalBaseMapSharedPtr &locToGloMap);

            // Virtual prototypes

            virtual int v_GetNumElmts(void)
            {
                return (*m_exp).size();
            }

            virtual const Array<OneD,const boost::shared_ptr<ExpList> > &v_GetBndCondExpansions(void);

            virtual boost::shared_ptr<ExpList> &v_UpdateBndCondExpansion(int i);
			
            virtual boost::shared_ptr<ExpList1D> &v_GetTrace();

            virtual boost::shared_ptr<LocalToGlobalDGMap> &v_GetTraceMap();

            virtual void v_AddTraceIntegral(
                                            const Array<OneD, const NekDouble> &Fx,
                                            const Array<OneD, const NekDouble> &Fy,
                                            Array<OneD, NekDouble> &outarray);

            virtual void v_AddTraceIntegral(
                                            const Array<OneD, const NekDouble> &Fn,
                                            Array<OneD, NekDouble> &outarray);

            virtual void v_AddTraceBiIntegral(
                                              const Array<OneD, const NekDouble> &Fwd,
                                              const Array<OneD, const NekDouble> &Bwd,
                                              Array<OneD, NekDouble> &outarray);

            virtual void v_GetFwdBwdTracePhys(
                                              Array<OneD,NekDouble> &Fwd,
                                              Array<OneD,NekDouble> &Bwd);

            virtual void v_GetFwdBwdTracePhys(
                                              const Array<OneD,const NekDouble>  &field,
                                              Array<OneD,NekDouble> &Fwd,
                                              Array<OneD,NekDouble> &Bwd);

            virtual void v_ExtractTracePhys(
                                            Array<OneD,NekDouble> &outarray);

            virtual void v_ExtractTracePhys(
                                            const Array<OneD, const NekDouble> &inarray,
                                            Array<OneD,NekDouble> &outarray);

            virtual void v_MultiplyByInvMassMatrix(
                                                   const Array<OneD,const NekDouble> &inarray,
                                                   Array<OneD,      NekDouble> &outarray,
                                                   bool  UseContCoeffs);

            virtual void v_HelmSolve(
                                     const Array<OneD, const NekDouble> &inarray,
                                     Array<OneD,       NekDouble> &outarray,
                                     NekDouble lambda,
                                     const Array<OneD, const NekDouble> &Sigma,
                                     const Array<OneD, const Array<OneD, NekDouble> > &varcoeff);

            virtual void v_HelmSolveCG(
                                       const Array<OneD, const NekDouble> &inarray,
                                       Array<OneD,       NekDouble> &outarray,
                                       NekDouble lambda,
                                       const Array<OneD, const NekDouble> &Sigma,
                                       const Array<OneD, const Array<OneD, NekDouble> > &varcoeff,
                                       bool UseContCoeffs,
                                       const Array<OneD, const NekDouble> &dirForcing);

            virtual void v_HelmSolveDG(
                                       const Array<OneD, const NekDouble> &inarray,
                                       Array<OneD,       NekDouble> &outarray,
                                       NekDouble lambda,
                                       const Array<OneD, const NekDouble> &Sigma,
                                       const Array<OneD, const Array<OneD, NekDouble> > &varcoeff,
                                       NekDouble tau);

            virtual void v_HelmSolve(
                                     const Array<OneD, const NekDouble> &inarray,
                                     Array<OneD,       NekDouble> &outarray,
                                     const Array<OneD, const Array<OneD, NekDouble> > &varcoeffs,
                                     const Array<OneD, NekDouble> &lambda,
                                     NekDouble tau);

            virtual void v_LinearAdvectionDiffusionReactionSolve(
                                                                 const Array<OneD, Array<OneD, NekDouble> > &velocity,
                                                                 const Array<OneD, const NekDouble> &inarray,
                                                                 Array<OneD, NekDouble> &outarray,
                                                                 const NekDouble lambda,
                                                                 bool  UseContCoeffs = false,
                                                                 const Array<OneD, const NekDouble>&
                                                                 dirForcing = NullNekDouble1DArray);

            virtual void v_LinearAdvectionReactionSolve(
                                                        const Array<OneD, Array<OneD, NekDouble> > &velocity,
                                                        const Array<OneD, const NekDouble> &inarray,
                                                        Array<OneD, NekDouble> &outarray,
                                                        const NekDouble lambda,
                                                        bool  UseContCoeffs = false,
                                                        const Array<OneD, const NekDouble>&
                                                        dirForcing = NullNekDouble1DArray);

            // wrapper functions about virtual functions
            virtual Array<OneD, NekDouble> &v_UpdateContCoeffs();

            virtual const
                Array<OneD, const NekDouble> &v_GetContCoeffs() const;

            virtual void v_LocalToGlobal();

            virtual void v_GlobalToLocal();

            virtual void v_BwdTrans(
                                    const Array<OneD,const NekDouble> &inarray,
                                    Array<OneD,      NekDouble> &outarray,
                                    bool  UseContCoeffs);

  	    virtual void v_SetUpPhysTangents(
		                	    const StdRegions::StdExpansionVector &locexp);
	    
            virtual void v_FwdTrans(
                                    const Array<OneD,const NekDouble> &inarray,
                                    Array<OneD,      NekDouble> &outarray,
                                    bool  UseContCoeffs);


            virtual void v_IProductWRTBase(
                                           const Array<OneD,const NekDouble> &inarray,
                                           Array<OneD,      NekDouble> &outarray,
                                           bool  UseContCoeffs);

            virtual void v_GeneralMatrixOp(
                                           const GlobalMatrixKey             &gkey,
                                           const Array<OneD,const NekDouble> &inarray,
                                           Array<OneD,      NekDouble> &outarray,
                                           bool  UseContCoeffs);

            virtual void v_GetCoords(Array<OneD, NekDouble> &coord_0,
                                     Array<OneD, NekDouble> &coord_1,
                                     Array<OneD, NekDouble> &coord_2 = NullNekDouble1DArray);
			
			virtual void v_GetCoords(NekDouble &x,NekDouble &y,NekDouble &z);
			
			virtual void v_GetCoord(Array<OneD, NekDouble> &coords);

			virtual void v_SetCoeff(NekDouble val);
			
			virtual void v_SetPhys(NekDouble val);
			
			virtual const SpatialDomains::VertexComponentSharedPtr &v_GetGeom(void) const;
			
			virtual const SpatialDomains::VertexComponentSharedPtr &v_GetVertex(void) const;
			
			virtual void v_PhysDeriv(const Array<OneD, const NekDouble> &inarray,
									 Array<OneD, NekDouble> &out_d0,
									 Array<OneD, NekDouble> &out_d1, 
									 Array<OneD, NekDouble> &out_d2, bool UseContCoeffs = false);
			
			virtual void v_PhysDeriv(Direction edir, 
									 const Array<OneD, const NekDouble> &inarray,
									 Array<OneD, NekDouble> &out_d, bool UseContCoeffs = false);			
			virtual void v_HomogeneousFwdTrans(const Array<OneD, const NekDouble> &inarray, 
												 Array<OneD, NekDouble> &outarray, 
												 bool UseContCoeffs = false);
			
			virtual void v_HomogeneousBwdTrans(const Array<OneD, const NekDouble> &inarray, 
												 Array<OneD, NekDouble> &outarray, 
												 bool UseContCoeffs = false);
			
			virtual void v_GetBCValues(Array<OneD, NekDouble> &BndVals, 
									   const Array<OneD, NekDouble> &TotField, 
									   int BndID);
			
			virtual void v_NormVectorIProductWRTBase(Array<OneD, const NekDouble> &V1,
												     Array<OneD, const NekDouble> &V2,
												     Array<OneD, NekDouble> &outarray,
												     int BndID);
			
            virtual void v_SetUpPhysNormals(const StdRegions::StdExpansionVector &locexp);

            virtual void v_SetUpTangents();

            virtual void v_GetBoundaryToElmtMap(Array<OneD, int> &ElmtID,
                                                Array<OneD,int> &EdgeID);

            virtual void v_ReadGlobalOptimizationParameters(
                                                            const std::string &infilename);

            virtual std::vector<SpatialDomains::FieldDefinitionsSharedPtr>
                v_GetFieldDefinitions(void);

            virtual void  v_GetFieldDefinitions(std::vector<SpatialDomains::FieldDefinitionsSharedPtr> &fielddef);


            virtual void v_AppendFieldData(SpatialDomains::FieldDefinitionsSharedPtr &fielddef, std::vector<NekDouble> &fielddata);

            virtual void v_AppendFieldData(SpatialDomains::FieldDefinitionsSharedPtr &fielddef, std::vector<NekDouble> &fielddata, Array<OneD, NekDouble> &coeffs);

            virtual void v_ExtractDataToCoeffs(SpatialDomains::FieldDefinitionsSharedPtr &fielddef, std::vector<NekDouble> &fielddata, std::string &field);

            virtual void v_ExtractDataToCoeffs(SpatialDomains::FieldDefinitionsSharedPtr &fielddef, std::vector<NekDouble> &fielddata, std::string &field,
                                               Array<OneD, NekDouble> &coeffs);

            virtual void v_WriteTecplotHeader(std::ofstream &outfile,
                                            std::string var = "v");
            virtual void v_WriteTecplotZone(std::ofstream &outfile,
                                            int expansion);
            virtual void v_WriteTecplotField(std::ofstream &outfile,
                                             int expansion);

            virtual void v_WriteVtkPieceHeader(std::ofstream &outfile, int expansion);
            virtual void v_WriteVtkPieceData(std::ofstream &outfile, int expansion,
                                        std::string var);

            virtual NekDouble v_L2(void);
            virtual NekDouble v_L2(const Array<OneD, const NekDouble> &soln);

            virtual void v_SetContCoeffsArray(Array<OneD, NekDouble> &inarray);

        private:
            virtual int v_GetContNcoeffs() const;

            virtual const Array<OneD,const SpatialDomains::BoundaryConditionShPtr> &v_GetBndConditions();

            virtual Array<OneD, SpatialDomains::BoundaryConditionShPtr> &v_UpdateBndConditions();

            virtual void v_EvaluateBoundaryConditions(const NekDouble time = 0.0, 
													  const NekDouble x2_in = NekConstants::kNekUnsetDouble,
													  const NekDouble x3_in = NekConstants::kNekUnsetDouble);

            virtual map<int, RobinBCInfoSharedPtr> v_GetRobinBCInfo(void);


            virtual void v_GetPeriodicEdges(SpatialDomains::MeshGraph2D &graph2D,
                                            SpatialDomains::BoundaryConditions &bcs,
                                            const std::string variable,
                                            vector<map<int,int> > & periodicVertices,
                                            map<int,int>& periodicEdges);

            // Homogeneous direction wrapper functions. 
            virtual LibUtilities::BasisSharedPtr  v_GetHomogeneousBasis(void)
            {
                ASSERTL0(false,
                         "This method is not defined or valid for this class type");
                return LibUtilities::NullBasisSharedPtr; 
            }

        };


        /// Shared pointer to an ExpList object.
        typedef boost::shared_ptr<ExpList>      ExpListSharedPtr;
        /// An empty ExpList object.
        static ExpList NullExpList;


        // Inline routines follow.

        /**
         * Returns the total number of local degrees of freedom
         * \f$N_{\mathrm{eof}}=\sum_{e=1}^{{N_{\mathrm{el}}}}N^{e}_m\f$.
         */
        inline int ExpList::GetNcoeffs() const
        {
            return m_ncoeffs;
        }

        inline int ExpList::GetNcoeffs(const int eid) const
        {
            return (*m_exp)[eid]->GetNcoeffs();
        }

        inline int ExpList::GetContNcoeffs() const
        {
            return v_GetContNcoeffs();
        }

        /**
         * Evaulates the maximum number of modes in the elemental basis
         * order over all elements
         */
        inline int ExpList::EvalBasisNumModesMax() const
        {
            int i;
            int returnval = 0;

            for(i= 0; i < (*m_exp).size(); ++i)
            {
                returnval = max(returnval,
                                (*m_exp)[i]->EvalBasisNumModesMax());
            }

            return returnval;
        }

        /**
         *
         */
        inline const Array<OneD,int> ExpList::EvalBasisNumModesMaxPerExp(void)
            const
        {
            int i;
            Array<OneD,int> returnval((*m_exp).size(),0);

            for(i= 0; i < (*m_exp).size(); ++i)
            {
                returnval[i]
                    = max(returnval[i],(*m_exp)[i]->EvalBasisNumModesMax());
            }

            return returnval;
        }


        /**
         *
         */
        inline int ExpList::GetTotPoints() const
        {
            return m_npoints;
        }

        inline int ExpList::GetTotPoints(const int eid) const
        {
            return (*m_exp)[eid]->GetTotPoints();
        }

        /**
         *
         */
        inline int ExpList::GetNpoints() const
        {
            return m_npoints;
        }

        /**
         *
         */
        inline void ExpList::SetTransState(const TransState transState)
        {
            m_transState = transState;
        }

        /**
         *
         */
        inline TransState ExpList::GetTransState() const
        {
            return m_transState;
        }
		
		/**
         *
         */
        inline void ExpList::SetFourierSpace(const FourierSpaceType fourierspace)
        {
            m_FourierSpace = fourierspace;
        }
		
        /**
         *
         */
        inline FourierSpaceType ExpList::GetFourierSpace() const
        {
            return m_FourierSpace;
        }

        /**
         * This function fills the array \f$\boldsymbol{u}_l\f$, the evaluation
         * of the expansion at the quadrature points (implemented as #m_phys),
         * with the values of the array \a inarray.
         *
         * @param   inarray         The array containing the values where
         *                          #m_phys should be filled with.
         */
        inline void ExpList::SetPhys(
                                     const Array<OneD, const NekDouble> &inarray)
        {
            ASSERTL0(inarray.num_elements() == m_npoints,
                     "Input array does not have correct number of elements.");

            Vmath::Vcopy(m_npoints,&inarray[0],1,&m_phys[0],1);
            m_physState = true;
        }


        inline void ExpList::SetPhysArray(Array<OneD, NekDouble> &inarray)
        {
            m_phys = inarray;
        }


        /**
         * @param   physState       \a true (=filled) or \a false (=not filled).
         */
        inline void ExpList::SetPhysState(const bool physState)
        {
            m_physState = physState;
        }


        /**
         * @return  physState       \a true (=filled) or \a false (=not filled).
         */
        inline bool ExpList::GetPhysState() const
        {
            return m_physState;
        }

        /**
         *
         */
        inline void ExpList::IProductWRTBase(
                                             const Array<OneD, const NekDouble> &inarray,
                                             Array<OneD,       NekDouble> &outarray,
                                             bool  UseContCoeffs)
        {
            v_IProductWRTBase(inarray,outarray,UseContCoeffs);
        }


        /**
         *
         */
        inline void ExpList::FwdTrans(
                                      const Array<OneD, const NekDouble> &inarray,
                                      Array<OneD,       NekDouble> &outarray,
                                      bool  UseContCoeffs)
        {
            v_FwdTrans(inarray,outarray,UseContCoeffs);
        }


        /**
         *
         */
        inline void ExpList::MultiplyByInvMassMatrix(
                                                     const Array<OneD,const NekDouble> &inarray,
                                                     Array<OneD,      NekDouble> &outarray,
                                                     bool  UseContCoeffs)
        {
            v_MultiplyByInvMassMatrix(inarray,outarray,UseContCoeffs);
        }

        /**
         *
         */
        inline void ExpList::HelmSolve(
                                       const Array<OneD, const NekDouble> &inarray,
                                       Array<OneD,       NekDouble> &outarray,
                                       NekDouble lambda,
                                       const Array<OneD, const NekDouble> &Sigma,
                                       const Array<OneD, const Array<OneD, NekDouble> > &varcoeff)
        {
            // HelmSolve(inarray, outarray);
            // HelmSolve(inarray, outarray, lambda);
            // HelmSolve(inarray, outarray, lambda, Lambda);
            // HelmSolve(inarray, outarray, lambda, Lambda, varcoeff);
            v_HelmSolve(inarray, outarray, lambda, Sigma, varcoeff);
            // v_HelmSolve -> v_HelmSolveCG or v_HelmSolveDG in derived classes
        }


        /**
         *
         */
        inline void ExpList::HelmSolve(
                                       const Array<OneD, const NekDouble> &inarray,
                                       Array<OneD,       NekDouble> &outarray,
                                       NekDouble lambda,
                                       bool UseContCoeffs,
                                       const Array<OneD, const NekDouble> &dirForcing,
                                       const Array<OneD, const NekDouble> &Sigma,
                                       const Array<OneD, const Array<OneD, NekDouble> > &varcoeff)
        {
            // HelmSolve(inarray, outarray, lambda, useContCoeff);
            // HelmSolve(inarray, outarray, lambda, useContCoeff, dirForcing);
            // HelmSolve(inarray, outarray, lambda, useContCoeff, dirForcing,
            //                                              Lambda);
            // HelmSolve(inarray, outarray, lambda, useContCoeff, dirForcing,
            //                                              Lambda, varcoeff);
            v_HelmSolveCG(inarray, outarray, lambda, Sigma, varcoeff,
                          UseContCoeffs, dirForcing);
        }


        /**
         *
         */
        inline void ExpList::HelmSolve(
                                       const Array<OneD, const NekDouble> &inarray,
                                       Array<OneD,       NekDouble> &outarray,
                                       NekDouble lambda,
                                       NekDouble tau,
                                       const Array<OneD, const NekDouble> &Sigma,
                                       const Array<OneD, const Array<OneD, NekDouble> > &varcoeff)
        {
            // HelmSolve(inarray, outarray, lambda, tau);
            // HelmSolve(inarray, outarray, lambda, tau, Lambda);
            // HelmSolve(inarray, outarray, lambda, tau, Lambda, varcoeff);
            v_HelmSolveDG(inarray, outarray, lambda, Sigma, varcoeff, tau);
        }

        /**
         *
         */
        inline void ExpList::LinearAdvectionDiffusionReactionSolve(
                                                                   const Array<OneD, Array<OneD, NekDouble> > &velocity,
                                                                   const Array<OneD, const NekDouble> &inarray,
                                                                   Array<OneD, NekDouble> &outarray,
                                                                   const NekDouble lambda,
                                                                   bool  UseContCoeffs,
                                                                   const Array<OneD, const NekDouble>&  dirForcing)
        {
            v_LinearAdvectionDiffusionReactionSolve(velocity,inarray, outarray, lambda, UseContCoeffs,dirForcing);
        }

        inline void ExpList::LinearAdvectionReactionSolve(
                                                          const Array<OneD, Array<OneD, NekDouble> > &velocity,
                                                          const Array<OneD, const NekDouble> &inarray,
                                                          Array<OneD, NekDouble> &outarray,
                                                          const NekDouble lambda,
                                                          bool  UseContCoeffs,
                                                          const Array<OneD, const NekDouble>&  dirForcing)
        {
            v_LinearAdvectionReactionSolve(velocity,inarray, outarray, lambda, UseContCoeffs,dirForcing);
        }

        /**
         *
         */
        inline void ExpList::BwdTrans (
                                       const Array<OneD, const NekDouble> &inarray,
                                       Array<OneD,       NekDouble> &outarray,
                                       bool  UseContCoeffs)
        {
            v_BwdTrans(inarray,outarray,UseContCoeffs);
        }


        /**
         *
         */
        inline void ExpList::GetCoords(Array<OneD, NekDouble> &coord_0,
                                       Array<OneD, NekDouble> &coord_1,
                                       Array<OneD, NekDouble> &coord_2)

        {
            v_GetCoords(coord_0,coord_1,coord_2);
        }
		
		/**
         *
         */
        inline void ExpList::SetCoeff(NekDouble val)
		
        {
            v_SetCoeff(val);
        }
		
		/**
         *
         */
		inline const SpatialDomains::VertexComponentSharedPtr &ExpList::GetGeom(void) const
		{
			return v_GetGeom();
		}
		
		/**
         *
         */
		inline const SpatialDomains::VertexComponentSharedPtr &ExpList::GetVertex(void) const
		{
			return v_GetVertex();
		}
		
		
		/**
         *
         */
        inline void ExpList::SetPhys(NekDouble val)
		
        {
            v_SetPhys(val);
        }
		
		/**
         *
         */
        inline void ExpList::GetCoords(NekDouble &x,NekDouble &y,NekDouble &z)
        {
            v_GetCoords(x,y,z);
        }
		
		inline void ExpList::GetCoord(Array<OneD, NekDouble> &coords)
        {
            v_GetCoord(coords);
        }

		/**
		 *
		 */
		inline void ExpList::PhysDeriv(const Array<OneD, const NekDouble> &inarray,
									   Array<OneD, NekDouble> &out_d0,
									   Array<OneD, NekDouble> &out_d1, 
									   Array<OneD, NekDouble> &out_d2, bool UseContCoeffs)
		{
			v_PhysDeriv(inarray,out_d0,out_d1,out_d2,UseContCoeffs);
		}

	inline void ExpList::PhysDeriv(Direction edir,
							   const Array<OneD, const NekDouble> &inarray,
							   Array<OneD, NekDouble> &out_d, bool UseContCoeffs)
	{
	    v_PhysDeriv(edir, inarray,out_d, UseContCoeffs);
	}		
		
		/**
		 *
		 */
		inline void ExpList::PhysDeriv(const int dir,
									   const Array<OneD, const NekDouble> &inarray,
									   Array<OneD, NekDouble> &out_d, bool UseContCoeffs)
		{
			MultiRegions::Direction edir =(MultiRegions::Direction)dir;
			v_PhysDeriv(edir,inarray, out_d,UseContCoeffs);
		}
		
		/**
		 *
		 */
		inline void ExpList::HomogeneousFwdTrans(const Array<OneD, const NekDouble> &inarray, 
												   Array<OneD, NekDouble> &outarray, 
												   bool UseContCoeffs)
		{
			v_HomogeneousFwdTrans(inarray,outarray,UseContCoeffs);
		}
		
		/**
		 *
		 */
		inline void ExpList::HomogeneousBwdTrans(const Array<OneD, const NekDouble> &inarray, 
												   Array<OneD, NekDouble> &outarray, 
												   bool UseContCoeffs)
		{
			v_HomogeneousBwdTrans(inarray,outarray,UseContCoeffs);
		}
		
		/**
		 *
		 */
		inline void ExpList::GetBCValues(Array<OneD, NekDouble> &BndVals, 
						                 const Array<OneD, NekDouble> &TotField, 
						                 int BndID)
		{
			v_GetBCValues(BndVals,TotField,BndID);
		}
		
		/**
		 *
		 */
		inline void ExpList::NormVectorIProductWRTBase(Array<OneD, const NekDouble> &V1,
									                   Array<OneD, const NekDouble> &V2,
									                   Array<OneD, NekDouble> &outarray,
									                   int BndID)
		{
			v_NormVectorIProductWRTBase(V1,V2,outarray,BndID);
		}
		
        /**
         * @param   eid         The index of the element to be checked.
         * @return  The dimension of the coordinates of the specific element.
         */
        inline int ExpList::GetCoordim(int eid)
        {
            ASSERTL2(eid <= (*m_exp).size(),
                     "eid is larger than number of elements");
            return (*m_exp)[eid]->GetCoordim();
        }

        /**
         * @param   i           The index of m_coeffs to be set
         * @param   val         The value which m_coeffs[i] is to be set to.
         */
        inline void ExpList::SetCoeff(int i, NekDouble val)
        {
            m_coeffs[i] = val;
        }


        /**
         * @param   i           The index of #m_coeffs to be set.
         * @param   val         The value which #m_coeffs[i] is to be set to.
         */
        inline void ExpList::SetCoeffs(int i, NekDouble val)
        {
            m_coeffs[i] = val;
        }


        inline void ExpList::SetCoeffsArray(Array<OneD, NekDouble> &inarray)
        {
            m_coeffs = inarray;
        }

        inline void ExpList::SetContCoeffsArray(Array<OneD, NekDouble> &inarray)
        {
            v_SetContCoeffsArray(inarray);
        }

        /**
         * As the function returns a constant reference to a
         * <em>const Array</em>, it is not possible to modify the
         * underlying data of the array #m_coeffs. In order to do
         * so, use the function #UpdateCoeffs instead.
         *
         * @return  (A constant reference to) the array #m_coeffs.
         */
        inline const Array<OneD, const NekDouble> &ExpList::GetCoeffs() const
        {
            return m_coeffs;
        }

        inline const Array<OneD, const NekDouble> &ExpList::GetContCoeffs()
                                                                        const
        {
            return v_GetContCoeffs();
        }

        inline void ExpList::LocalToGlobal()
        {
            v_LocalToGlobal();
        }

        inline void ExpList::GlobalToLocal()
        {
            v_GlobalToLocal();
        }


        /**
         * @param   i               The index of #m_coeffs to be returned
         * @return  The NekDouble held in #m_coeffs[i].
         */
        inline NekDouble ExpList::GetCoeff(int i)
        {
            return m_coeffs[i];
        }

        /**
         * @param   i               The index of #m_coeffs to be returned
         * @return  The NekDouble held in #m_coeffs[i].
         */
        inline NekDouble ExpList::GetCoeffs(int i)
        {
            return m_coeffs[i];
        }

        /**
         * As the function returns a constant reference to a
         * <em>const Array</em> it is not possible to modify the
         * underlying data of the array #m_phys. In order to do
         * so, use the function #UpdatePhys instead.
         *
         * @return  (A constant reference to) the array #m_phys.
         */
        inline const Array<OneD, const NekDouble> &ExpList::GetPhys()  const
        {
            return m_phys;
        }

        /**
         * @return  \f$N_{\mathrm{el}}\f$, the number of elements in the
         *          expansion.
         */
        inline int ExpList::GetExpSize(void)
        {
            return (*m_exp).size();
        }


        /**
         * @param   n               The index of the element concerned.
         *
         * @return  (A shared pointer to) the local expansion of the
         *          \f$n^{\mathrm{th}}\f$ element.
         */
        inline StdRegions::StdExpansionSharedPtr& ExpList::GetExp(int n) const
        {
            return (*m_exp)[n];
        }

        /**
         * @return  (A const shared pointer to) the local expansion vector #m_exp
         */
        inline const boost::shared_ptr<StdRegions::StdExpansionVector> ExpList::GetExp(void) const
        {
            return m_exp;
        }


        /**
         *
         */
        inline const int ExpList::GetCoeff_Offset(int n) const
        {
            return m_coeff_offset[n];
        }

        /**
         *
         */
        inline const int ExpList::GetPhys_Offset(int n) const
        {
            return m_phys_offset[n];
        }

        /**
         *
         */
        inline const int ExpList::GetOffset_Elmt_Id(int n) const
        {
            return m_offset_elmt_id[n];
        }

        /**
         * If one wants to get hold of the underlying data without modifying
         * them, rather use the function #GetCoeffs insead.
         *
         * @return  (A reference to) the array #m_coeffs.
         */
        inline Array<OneD, NekDouble> &ExpList::UpdateCoeffs()
        {
            m_transState = eLocal;
            return m_coeffs;
        }

        /**
         * If one wants to get hold of the underlying data without modifying
         * them,  rather use the function #GetPhys instead.
         *
         * @return  (A reference to) the array #m_phys.
         */
        inline Array<OneD, NekDouble> &ExpList::UpdatePhys()
        {
            m_physState = true;
            return m_phys;
        }


        // functions associated with DisContField
        inline const Array<OneD, const  boost::shared_ptr<ExpList> > &ExpList::GetBndCondExpansions()
        {
            return v_GetBndCondExpansions();
        }

        inline boost::shared_ptr<ExpList>  &ExpList::UpdateBndCondExpansion(int i)
        {
            return v_UpdateBndCondExpansion(i);
        }
		
        inline boost::shared_ptr<ExpList1D> &ExpList::GetTrace()
        {
            return v_GetTrace();
        }

        inline boost::shared_ptr<LocalToGlobalDGMap> &ExpList::GetTraceMap()
        {
            return v_GetTraceMap();
        }

        inline void ExpList::AddTraceIntegral(
                                const Array<OneD, const NekDouble> &Fx,
                                const Array<OneD, const NekDouble> &Fy,
                                      Array<OneD, NekDouble> &outarray)
        {
            v_AddTraceIntegral(Fx,Fy,outarray);
        }

        inline void ExpList::AddTraceIntegral(
                                const Array<OneD, const NekDouble> &Fn,
                                      Array<OneD, NekDouble> &outarray)
        {
            v_AddTraceIntegral(Fn,outarray);
        }

        inline void ExpList::AddTraceBiIntegral(
                                    const Array<OneD, const NekDouble> &Fwd,
                                    const Array<OneD, const NekDouble> &Bwd,
                                          Array<OneD, NekDouble> &outarray)
        {
            v_AddTraceBiIntegral(Fwd,Bwd,outarray);
        }

        inline void ExpList::GetFwdBwdTracePhys(Array<OneD,NekDouble> &Fwd,
                                         Array<OneD,NekDouble> &Bwd)
        {
            v_GetFwdBwdTracePhys(Fwd,Bwd);
        }

        inline void ExpList::GetFwdBwdTracePhys(
                                const Array<OneD,const NekDouble>  &field,
                                      Array<OneD,NekDouble> &Fwd,
                                      Array<OneD,NekDouble> &Bwd)
        {
            v_GetFwdBwdTracePhys(field,Fwd,Bwd);
        }

        inline void ExpList::ExtractTracePhys(Array<OneD,NekDouble> &outarray)
        {
            v_ExtractTracePhys(outarray);
        }


        inline void ExpList::ExtractTracePhys(
                            const Array<OneD, const NekDouble> &inarray,
                                  Array<OneD,NekDouble> &outarray)
        {
            v_ExtractTracePhys(inarray,outarray);
        }

        inline const Array<OneD,const SpatialDomains::BoundaryConditionShPtr>
                            &ExpList::GetBndConditions()
        {
            return v_GetBndConditions();
        }


        inline Array<OneD, SpatialDomains::BoundaryConditionShPtr>
            &ExpList::UpdateBndConditions()
        {
            return v_UpdateBndConditions();
        }

        inline void ExpList::EvaluateBoundaryConditions(const NekDouble time,
                                                        const NekDouble x2_in,
														const NekDouble x3_in)
        {
            v_EvaluateBoundaryConditions(time,x2_in,x3_in);
        }

        // Routines for continous matrix solution
        /**
         * This operation is equivalent to the evaluation of
         * \f$\underline{\boldsymbol{M}}^e\boldsymbol{\hat{u}}_l\f$, that is,
         * \f[ \left[
         * \begin{array}{cccc}
         * \boldsymbol{M}^1 & 0 & \hspace{3mm}0 \hspace{3mm}& 0 \\
         * 0 & \boldsymbol{M}^2 & 0 & 0 \\
         * 0 &  0 & \ddots &  0 \\
         * 0 &  0 & 0 & \boldsymbol{M}^{N_{\mathrm{el}}} \end{array} \right]
         *\left [ \begin{array}{c}
         * \boldsymbol{\hat{u}}^{1} \\
         * \boldsymbol{\hat{u}}^{2} \\
         * \vdots \\
         * \boldsymbol{\hat{u}}^{{{N_{\mathrm{el}}}}} \end{array} \right ]\f]
         * where \f$\boldsymbol{M}^e\f$ are the local matrices of type
         * specified by the key \a mkey. The decoupling of the local matrices
         * allows for a local evaluation of the operation. However, rather than
         * a local matrix-vector multiplication, the local operations are
         * evaluated as implemented in the function
         * StdRegions#StdExpansion#GeneralMatrixOp.
         *
         * @param   mkey            This key uniquely defines the type matrix
         *                          required for the operation.
         * @param   inarray         The vector \f$\boldsymbol{\hat{u}}_l\f$ of
         *                          size \f$N_{\mathrm{eof}}\f$.
         * @param   outarray        The resulting vector of size
         *                          \f$N_{\mathrm{eof}}\f$.
         */
        inline void ExpList::GeneralMatrixOp(
                                const GlobalMatrixKey             &gkey,
                                const Array<OneD,const NekDouble> &inarray,
                                      Array<OneD,      NekDouble> &outarray,
                                bool  UseContCoeffs)
        {
            v_GeneralMatrixOp(gkey,inarray,outarray,UseContCoeffs);
        }


        inline void ExpList::SetUpPhysNormals(
                                const StdRegions::StdExpansionVector &locexp)
        {
            v_SetUpPhysNormals(locexp);
        }

        inline void ExpList::SetUpPhysTangents(
                                const StdRegions::StdExpansionVector &locexp)
        {
            v_SetUpPhysTangents(locexp);
        }
        
        inline void ExpList::SetUpTangents()
        {
            v_SetUpTangents();
        }

        inline void ExpList::GetBoundaryToElmtMap( Array<OneD, int> &ElmtID,
                                            Array<OneD,int> &EdgeID)
        {
            v_GetBoundaryToElmtMap(ElmtID,EdgeID);
        }

        const static Array<OneD, ExpListSharedPtr> NullExpListSharedPtrArray;
        
  } //end of namespace
} //end of namespace

#endif // EXPLIST_H

/**
* $Log: ExpList.h,v $
* Revision 1.89  2010/03/02 23:50:23  sherwin
* Updates related to making IncNavierStokesSolver able to use ContCoeffs
*
* Revision 1.88  2010/03/01 17:57:28  cantwell
* Fixed 3D global matrix operations.
* Fixed ProjectCont{1,2,3}D demos.
* Fixed incorrectly placed ASSERT in boundary conditions.
* Updated TimingGeneralMatrixOp3D to use contfield3d rather than explist3d.
*
* Revision 1.87  2010/01/27 13:19:13  cantwell
* Added functions to write history/probe data during timestepping.
*
* Revision 1.86  2010/01/20 18:05:09  cantwell
* Added utility for probing a line of points in a FLD file.
*
* Revision 1.85  2010/01/06 13:24:48  cantwell
* Corrected naming of parameters in HelmSolve routines.
* Documentation of DisContField1D.
*
* Revision 1.84  2010/01/03 19:39:09  cantwell
* Added FldToVtk converter.
* Added XmlToVtk converter.
*
* Revision 1.83  2009/12/15 18:09:02  cantwell
* Split GeomFactors into 1D, 2D and 3D
* Added generation of tangential basis into GeomFactors
* Updated ADR2DManifold solver to use GeomFactors for tangents
* Added <GEOMINFO> XML session section support in MeshGraph
* Fixed const-correctness in VmathArray
* Cleaned up LocalRegions code to generate GeomFactors
* Removed GenSegExp
* Temporary fix to SubStructuredGraph
* Documentation for GlobalLinSys and GlobalMatrix classes
*
* Revision 1.82  2009/12/14 18:01:08  cbiotto
* Adding functions for printing out tecplot file
*
* Revision 1.81  2009/12/08 15:10:50  sehunchun
* HelmholtzSolver with additional variables are added
*
* Revision 1.80  2009/11/19 14:06:00  sehunchun
* *** empty log message ***
*
* Revision 1.79  2009/11/19 13:48:06  sehunchun
* *** empty log message ***
*
* Revision 1.78  2009/11/19 11:41:07  pvos
* Fixed various bugs
*
* Revision 1.77  2009/11/10 19:05:34  sehunchun
* *** empty log message ***
*
* Revision 1.76  2009/11/07 17:15:17  sehunchun
* Add GetTotPoints(idx)
*
* Revision 1.75  2009/11/06 21:51:18  sherwin
* Added L2_DGDeriv method
*
* Revision 1.74  2009/11/04 20:30:15  cantwell
* Added documentation to ExpList and ExpList1D and tidied up code.
*
* Revision 1.73  2009/11/02 19:15:43  cantwell
* Moved ContField1D to inherit from DisContField1D.
* Moved ContField3D to inherit from DisContField3D.
* Incorporated GenExpList1D functionality into ExpList1D.
* Tidied up and added documentation to various classes.
* Moved Namespace documentation and introductions to separate files along with
* doxygen configuration.
* Added option to use system ZLIB library instead of libboost_zlib on UNIX.
* Added extra search paths to FindMetis.cmake and FindNektar++.cmake.
* Updated Linux compiling instructions.
* Updated regDemo to use Helmholtz2D-g when built as debug.
*
* Revision 1.72  2009/10/30 14:02:55  pvos
* Multi-level static condensation updates
*
* Revision 1.71  2009/10/25 18:54:38  sherwin
* Added H1 norm for error evaluation
*
* Revision 1.70  2009/10/22 16:46:15  cbiotto
* Adding function EvalBasisNumModesMaxPerExp
*
* Revision 1.69  2009/10/22 16:40:35  cbiotto
* *** empty log message ***
*
* Revision 1.68  2009/09/06 22:28:45  sherwin
* Updates for Navier-Stokes solver
*
* Revision 1.67  2009/07/08 17:22:47  sehunchun
* Deleting GetTanBasis
*
* Revision 1.66  2009/07/08 11:13:54  sehunchun
* Adding GetSurfaceNormal function
*
* Revision 1.65  2009/07/07 16:36:45  sehunchun
* Adding AddTraceBiIntegral...
*
* Revision 1.64  2009/07/03 15:38:25  sehunchun
* Adding GetTanBasis function
*
* Revision 1.63  2009/05/14 14:26:41  pvos
* Updates to apply the dirichlet boundary condition forcing inside the static condensation algorithm
*
* Revision 1.62  2009/04/27 21:34:58  sherwin
* Modified WriteToField Method
*
* Revision 1.61  2009/04/27 15:02:04  pvos
* From h-to-p efficiently updates
*
* Revision 1.60  2009/04/22 22:32:10  sherwin
* Added in method to read dat file
*
* Revision 1.59  2009/04/20 16:14:06  sherwin
* Updates for optimising bandwidth of DG solver and allowing write import on explist
*
* Revision 1.58  2009/04/03 20:33:57  sherwin
* Update for Eigenfunction evaluation
*
* Revision 1.57  2009/03/23 11:52:15  pvos
* NekMatrix updates
*
* Revision 1.56  2009/03/23 10:51:52  pvos
* Added BlockMatrix support
*
* Revision 1.55  2009/03/04 14:17:38  pvos
* Removed all methods that take and Expansion as argument
*
* Revision 1.54  2009/03/04 05:58:49  bnelson
* Fixed visual studio compile errors.
*
* Revision 1.53  2009/02/08 09:11:49  sherwin
* General updates to introduce multiple matrix definitions based on different boundary types
*
* Revision 1.52  2009/02/03 14:33:08  pvos
* Modifications for solvers with time-dependent dirichlet BC's
*
* Revision 1.51  2009/02/02 16:43:26  claes
* Added virtual functions for solver access
*
* Revision 1.50  2009/01/06 21:05:57  sherwin
* Added virtual function calls for BwdTrans, FwdTrans and IProductWRTBase from the class ExpList. Introduced _IterPerExp versions of these methods in ExpList.cpp§
*
* Revision 1.49  2008/11/19 10:52:55  pvos
* Changed MultiplyByInvMassMatrix + added some virtual functions
*
* Revision 1.48  2008/11/01 22:06:45  bnelson
* Fixed Visual Studio compile error.
*
* Revision 1.47  2008/10/29 22:46:35  sherwin
* Updates for const correctness and a few other bits
*
* Revision 1.46  2008/10/19 15:57:52  sherwin
* Added method EvalBasisNumModesMax
*
* Revision 1.45  2008/10/16 10:21:42  sherwin
* Updates to make methods consisten with AdvectionDiffusionReactionsSolver. Modified MultiplyByInvMassMatrix to take local or global arrays
*
* Revision 1.44  2008/10/04 20:04:26  sherwin
* Modifications for solver access
*
* Revision 1.43  2008/09/16 13:36:06  pvos
* Restructured the LocalToGlobalMap classes
*
* Revision 1.42  2008/08/14 22:15:51  sherwin
* Added LocalToglobalMap and DGMap and depracted LocalToGlobalBndryMap1D,2D. Made DisContField classes compatible with updated ContField formats
*
* Revision 1.41  2008/07/29 22:27:33  sherwin
* Updates for DG solvers, including using GenSegExp, fixed forcing function on UDG HelmSolve and started to tidy up the mapping arrays to be 1D rather than 2D
*
* Revision 1.40  2008/07/15 13:00:04  sherwin
* Updates for DG advection solver - not yet debugged
*
* Revision 1.39  2008/07/12 19:08:29  sherwin
* Modifications for DG advection routines
*
* Revision 1.38  2008/07/12 17:31:39  sherwin
* Added m_phys_offset and rename m_exp_offset to m_coeff_offset
*
* Revision 1.37  2008/07/11 15:48:32  pvos
* Added Advection classes
*
* Revision 1.36  2008/06/06 23:27:20  ehan
* Added doxygen documentation
*
* Revision 1.35  2008/06/05 15:06:58  pvos
* Added documentation
*
* Revision 1.34  2008/05/29 21:35:03  pvos
* Added WriteToFile routines for Gmsh output format + modification of BndCond implementation in MultiRegions
*
* Revision 1.33  2008/05/10 18:27:33  sherwin
* Modifications necessary for QuadExp Unified DG Solver
*
* Revision 1.32  2008/04/06 06:00:07  bnelson
* Changed ConstArray to Array<const>
*
* Revision 1.31  2008/03/12 15:25:45  pvos
* Clean up of the code
*
* Revision 1.30  2008/01/23 21:50:52  sherwin
* Update from EdgeComponents to SegGeoms
*
* Revision 1.29  2007/12/17 13:05:04  sherwin
* Made files compatible with modifications in StdMatrixKey which now holds constants
*
* Revision 1.28  2007/12/06 22:52:30  pvos
* 2D Helmholtz solver updates
*
* Revision 1.27  2007/11/20 16:27:16  sherwin
* Zero Dirichlet version of UDG Helmholtz solver
*
* Revision 1.26  2007/10/03 11:37:50  sherwin
* Updates relating to static condensation implementation
*
* Revision 1.25  2007/09/03 19:58:31  jfrazier
* Formatting.
*
* Revision 1.24  2007/07/27 03:10:49  bnelson
* Fixed g++ compile error.
*
* Revision 1.23  2007/07/26 08:40:49  sherwin
* Update to use generalised i/o hooks in Helmholtz1D
*
* Revision 1.22  2007/07/22 23:04:20  bnelson
* Backed out Nektar::ptr.
*
* Revision 1.21  2007/07/20 02:04:12  bnelson
* Replaced boost::shared_ptr with Nektar::ptr
*
* Revision 1.20  2007/07/17 07:11:05  sherwin
* Chaned definition of NullExpList
*
* Revision 1.19  2007/07/16 18:28:43  sherwin
* Modification to introduce non-zero Dirichlet boundary conditions into the Helmholtz1D Demo
*
* Revision 1.18  2007/07/13 09:02:24  sherwin
* Mods for Helmholtz solver
*
* Revision 1.17  2007/06/07 15:54:19  pvos
* Modificications to make Demos/MultiRegions/ProjectCont2D work correctly.
* Also made corrections to various ASSERTL2 calls
*
* Revision 1.16  2007/06/05 16:36:55  pvos
* Updated Explist2D ContExpList2D and corresponding demo-codes
*
**/

namespace Nektar {
/**
* \page pageCoupledSolver Direct solver (coupled approach)
 
 Here, we will explain how to use the incompressible Navier-Stokes solver of Nektar++ using the direct solver. 
 Its source code is located at <code>Nektar++/solvers/IncNavierStokesSolver/EquationSystems</code>.
 We will start by showing the formulation and demonstrating how to run a first example with the incompressible Navier-Stokes 
 solver and briefly explain the options of the solver specified in the input file. 
 After these introductory explanations, we will present some numerical results 
 in order to demonstrate the capabilities and the accuracy of the solver.
 - @ref sectionFormulation
 - @ref sectionRunning1example
 - @ref sectionInputOptions
 - @ref sectionResult
 - @ref sectionRef
 
 \section sectionFormulation Formulation
 
 We consider the weak form of the Stokes problem for the velocity field \f$\boldsymbol{v}=[u,v]^{T}\f$ and the pressure field \f$p\f$:
 \f[
 (\nabla \phi,\nu \nabla \boldsymbol{v}) - (\nabla\cdot\phi,p)=(\phi,\boldsymbol{f})
 \f]
 \f[
 (q,\nabla \cdot \boldsymbol{v}) = 0
 \f]
 
 where \f$\boldsymbol{v},\phi \in \boldsymbol{V}\f$, \f$p,q \in \boldsymbol{W}\f$ and \f$\boldsymbol{V},\boldsymbol{W}\f$ are appropriate spaces for the velocity
 and the pressure system to satisfy the inf-sup condition.
 Using a matrix notation the discrete system may be written as:
 \f{displaymath}
 \left[ \begin{array}{ccc}
 A & D_b^T & B\\
 D_b & 0 & D_i^T\\
 B^T & D_i & C
 \end{array}\right]
 \left[ \begin{array}{c}
 \boldsymbol{v_b}\\
 p\\
 \boldsymbol{v_i}
 \end{array}\right] =
 \left[ \begin{array}{c}
 \boldsymbol{f_b}\\
 0\\
 \boldsymbol{f_i}
 \end{array}\right]
 \f}
 
 where the components of \f$A,B\f$ and \f$C\f$ are \f$(\nabla\phi_b,\nu\nabla\boldsymbol{v_b})\f$, \f$(\nabla\phi_b,\nu\nabla\boldsymbol{v_i})\f$ and
 \f$(\nabla\phi_i,\nu\nabla\boldsymbol{v_i})\f$ and the components \f$D_b\f$ and \f$D_i\f$ are \f$(q,\nabla\boldsymbol{v_b})\f$ and \f$(q,\nabla\boldsymbol{v_i})\f$.
 The indices \f$b\f$ and \f$i\f$ refer to the degrees of freedom on the elemental boundary and interior, respectively. In constructing the system we have lumped the 
 contributions form each component of the velocity field into matrices \f$A,B\f$ and \f$C\f$. However we note that for a Newtonian fluid the contribution from
 each field is decoupled. Since the inetrior degrees of freedom of the velocity field do not overlap, the matrix  \f$C\f$ is block diagonal and to take advantage
 of this structure we can statically condense out the \f$C\f$ matrix to obtain the system:
 \f{displaymath}
 \left[ \begin{array}{ccc}
 A-BC^{-1}B^T & D_b^T-BC^{-1}D_i & 0\\
 D_b-D_i^TC^{-1}B^T & -D_i^TC^{-1}D_i & 0\\
 B^T & D_i & C
 \end{array}\right]
 \left[ \begin{array}{c}
 \boldsymbol{v_b}\\
 p\\
 \boldsymbol{v_i}
 \end{array}\right] =
 \left[ \begin{array}{c}
 \boldsymbol{f_b} - BC^{-1}\boldsymbol{f_i}\\
 -D_i^TC^{-1}\boldsymbol{f_i}\\
 \boldsymbol{f_i}
 \end{array}\right]
 \f}
 
 To extend the aboce Stokes solver to an unsteady Navier-Stokes solver we first introduce the unsteady term, \f$\boldsymbol{v_t}\f$, into the Stokes problem.
 This has the principal effect of modifying the weak Laplacian operator \f$(\nabla\phi,\nu\nabla\boldsymbol{v})\f$ into a weak Helmholts operator
 \f$(\nabla\phi,\nu\nabla\boldsymbol{v})-\lambda(\phi,\boldsymbol{v})\f$ where \f$\lambda\f$ depends on the time integration scheme. The second modification requires
 the explicit discretisation of the non-linear terms in a similar manner to the splitting scheme and this term is then introduced as the forcing term \f$\boldsymbol{f}\f$.
 
\section sectionRunning1example Running a first example
 
 The folder <code>Nektar++/regressionTests/Solvers/IncNavierStokesSolver/InputFiles</code> contains several <code>*.xml</code> files.
 These are input files for the Navier-Stokes solver specifying the geometry (i.e. the mesh and 
 the spectal/hp expansion), the parameters and boundary conditions. Further details on the structure 
 of the input file can be found in @ref pageXML.
 
 Now, lets try to run one of them with the coupled solver.
 
 - Copy the input file, <code>Test_ChanFlow_LinNS_m8.xml</code> to the directory where the solver is compiled, e.g. 
 <code>Nektar++/solvers/IncNavierStokesSolver/build/IncNavierStokesSolver</code>.
 - Run the code by typing
 @code
 ./IncNavierStokesSolver Test_ChanFlow_LinNS_m8.xml
 @endcode
 
 The solution should now have been written to the file <code>Test_ChanFlow_LinNS_m8.fld</code>. 
 This file is formatted in the Nektar++ output format.
 To visualise the solution, we can convert the fld-file into TecPlot, Gmsh or Tecplot file formats using 
 the Post-processing tools in <code>Nektar++/utilities/builds/PostProcessing/</code>. 
 Here, we will demonstrate how to convert the <code>fld</code>-file into TecPlot-file format. 
 
 We convert the <code>fld</code>-file into Tecplot-file format by typing
 @code
 ../../../utilities/builds/PostProcessing/FldToTecplot Test_ChanFlow_LinNS_m8.xml Test_ChanFlow_LinNS_m8.fld
 @endcode
 
 It will create <code>Test_ChanFlow_LinNS_m8.dat</code> which can be loaed in TecPlot.
 
 \image html CoupChanFlow.png "Channel Flow (u-velocity component)"
 
 

\section sectionInputOptions Input Options
 
\section sectionResult Numerical Results
 
\section sectionRef References
 
*/	
	
}

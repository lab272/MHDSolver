#include <cstdio>
#include <cstdlib>

#include <MultiRegions/MultiRegions.hpp>
#include <MultiRegions/ContExpList3D.h>

using namespace Nektar;

// This routine projects a polynomial which has energy in all mdoes of
// the expansions and report an error.

int main(int argc, char *argv[])
{ 
    MultiRegions::ContExpList3DSharedPtr Exp,Fce;
    int     i, j, nq,  coordim;
    Array<OneD,NekDouble>  fce; 
    Array<OneD,NekDouble>  xc0,xc1,xc2; 
    NekDouble  lambda;
    
    if(argc != 2)
    {
        fprintf(stderr,"Usage: ProjectCont3D  meshfile \n");
        exit(1);
    }

    //----------------------------------------------
    // Read in mesh from input file
    string meshfile(argv[1]);
    SpatialDomains::MeshGraph3D graph3D; 
    graph3D.ReadGeometry(meshfile);
    graph3D.ReadExpansions(meshfile);
    //----------------------------------------------

    //----------------------------------------------
    // Print summary of solution details
    const SpatialDomains::ExpansionVector &expansions = graph3D.GetExpansions();
    LibUtilities::BasisKey bkey = graph3D.GetBasisKey(expansions[0],0);
    int nmodes = (int) expansions[0]->m_NumModesEqn.Evaluate();
    cout << "Solving 3D C0 continuous Projection"  << endl; 
    cout << "    Expansion  : " << SpatialDomains::kExpansionTypeStr[expansions[0]->m_ExpansionType] << endl;
    cout << "    No. modes  : " << nmodes << endl;
    cout << endl;
    //----------------------------------------------
   
    //----------------------------------------------
    // Define Expansion 
    Exp = MemoryManager<MultiRegions::ContExpList3D>::AllocateSharedPtr(graph3D);
    //----------------------------------------------  
    
    //----------------------------------------------
    // Set up coordinates of mesh for Forcing function evaluation
    coordim = Exp->GetCoordim(0);
    nq      = Exp->GetTotPoints();
    
    xc0 = Array<OneD,NekDouble>(nq,0.0);
    xc1 = Array<OneD,NekDouble>(nq,0.0);
    xc2 = Array<OneD,NekDouble>(nq,0.0);
    
    switch(coordim)
    {
    case 1:
        Exp->GetCoords(xc0);
        break;
    case 2:
        Exp->GetCoords(xc0,xc1);
        break;
    case 3:
        Exp->GetCoords(xc0,xc1,xc2);
        break;
    }
    //----------------------------------------------

    //----------------------------------------------
    // Define forcing function
    fce = Array<OneD,NekDouble>(nq);  
 
    for(i = 0; i < nq; ++i)
    {
        fce[i] = 0.0;
        for(j = 0; j < nmodes; ++j)
        {
            fce[i] += pow(xc0[i],j);
            fce[i] += pow(xc1[i],j);
            fce[i] += pow(xc2[i],j);
        }
    }
    
    //---------------------------------------------
    // Set up ExpList1D containing the solution 
    Fce = MemoryManager<MultiRegions::ContExpList3D>::AllocateSharedPtr(*Exp);
    Fce->SetPhys(fce);
    //---------------------------------------------

    //---------------------------------------------
    // Project onto Expansion 
    Exp->FwdTrans(*Fce);
    //---------------------------------------------
    
    //-------------------------------------------
    // Backward Transform Solution to get projected values
    Exp->BwdTrans(*Exp);
    //-------------------------------------------  

    //----------------------------------------------
    // Write solution 
    ofstream outfile2("ProjectContFile3D.dat");
    Exp->WriteToFile(outfile2,eTecplot);
    outfile2.close();
    //----------------------------------------------
    
    //--------------------------------------------
    // Calculate L_inf error 
    cout << "L infinity error: " << Exp->Linf(*Fce) << endl;
    cout << "L 2 error:        " << Exp->L2  (*Fce) << endl;
    //--------------------------------------------

    return 0;
}

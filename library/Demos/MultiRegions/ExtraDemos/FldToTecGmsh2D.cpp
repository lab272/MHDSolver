#include <cstdio>
#include <cstdlib>

#include <MultiRegions/ExpList2D.h>

using namespace Nektar;

int main(int argc, char *argv[])
{
    MultiRegions::ExpList2DSharedPtr Exp;
    int     i, j, nq,  coordim;
    Array<OneD,NekDouble>  fce; 
    Array<OneD,NekDouble>  xc0,xc1,xc2; 

    if(argc != 3)
    {
#ifdef TECPLOT
        fprintf(stderr,"Usage: FldToTec2D  meshfile fieldfile\n");
#else
        fprintf(stderr,"Usage: FldToGmsh2D  meshfile fieldfile\n");
#endif
        exit(1);
    }

    //----------------------------------------------
    // Read in mesh from input file
    string meshfile(argv[argc-2]);
    SpatialDomains::MeshGraph2D graph2D; 
    graph2D.ReadGeometry(meshfile);
    //graph2D.ReadExpansions(meshfile);
    //----------------------------------------------
    
    //----------------------------------------------
    // Import field file. 
    string fieldfile(argv[argc-1]);
    vector<SpatialDomains::FieldDefinitionsSharedPtr> fielddef;
    vector<vector<NekDouble> > fielddata;
    graph2D.Import(fieldfile,fielddef,fielddata);
    //----------------------------------------------

    //----------------------------------------------
    // Set up Expansion information
    vector< vector<LibUtilities::PointsType> > pointstype;
    for(i = 0; i < fielddef.size(); ++i)
    {         vector<LibUtilities::PointsType> ptype;
        for(j = 0; j < 2; ++j)
        {
            ptype.push_back(LibUtilities::ePolyEvenlySpaced);
        }
        pointstype.push_back(ptype);
    }
    graph2D.SetExpansions(fielddef,pointstype);
    //----------------------------------------------

    //----------------------------------------------        
    // Define Expansion 
    Exp = MemoryManager<MultiRegions::ExpList2D>::AllocateSharedPtr(graph2D);
    //----------------------------------------------  

    //----------------------------------------------
    // Copy data to file 
    for(int i = 0; i < fielddata.size(); ++i)
    {
        Exp->ExtractDataToCoeffs(fielddef[i],fielddata[i],fielddef[i]->m_Fields[0]);
    }
    Exp->BwdTrans(Exp->GetCoeffs(),Exp->UpdatePhys());    
    //----------------------------------------------
    
    //----------------------------------------------
    // Write solution  depending on #define
#ifdef TECPLOT
    string   outfile(strtok(argv[argc-1],"."));
    string   endfile(".dat");
    outfile += endfile; 
    ofstream outstrm(outfile.c_str());
    
    Exp->WriteToFile(outstrm,eTecplot);
    outstrm.close();
#else
    string   outfile(strtok(argv[argc-1],"."));
    string   endfile(".pos");
    outfile += endfile; 
    ofstream outstrm(outfile.c_str());

    Exp->WriteToFile(outstrm,eGmsh);
    outstrm.close();
#endif
    //----------------------------------------------
    return 0;
}

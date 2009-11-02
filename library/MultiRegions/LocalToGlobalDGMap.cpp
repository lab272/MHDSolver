///////////////////////////////////////////////////////////////////////////////
//
// File LocToGlobalDGMap.cpp
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
// Description: Local to Global Base Class mapping routines
//
///////////////////////////////////////////////////////////////////////////////

#include <MultiRegions/LocalToGlobalDGMap.h>

#include <boost/config.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>


namespace Nektar
{
    namespace MultiRegions
    {
        LocalToGlobalDGMap::LocalToGlobalDGMap():
            m_numDirichletBndPhys(0)
        {
        }
        
        LocalToGlobalDGMap::~LocalToGlobalDGMap()
        {
        }

        LocalToGlobalDGMap::LocalToGlobalDGMap( const SpatialDomains::MeshGraph1D &graph1D,
                                                const boost::shared_ptr<StdRegions::StdExpansionVector> &exp1D,
                                                const GlobalSysSolnType solnType,
                                                const Array<OneD, const LocalRegions::PointExpSharedPtr> &bndCondExp,
                                                const Array<OneD, const SpatialDomains::BoundaryConditionShPtr> &bndCond)
        {
            int i,j; 
            int cnt, vid, gid;
            int nbnd = bndCondExp.num_elements();
            
            // set up Local to Continuous mapping 
            Array<OneD,unsigned int> vmap;
            LocalRegions::SegExpSharedPtr locSegExp;
            
            m_numGlobalBndCoeffs  = exp1D->size()+1;
            m_numGlobalCoeffs = m_numGlobalBndCoeffs;
            m_numLocalBndCoeffs = 2*exp1D->size();
            m_numLocalCoeffs = m_numLocalBndCoeffs;
            m_localToGlobalBndMap   = Array<OneD, int>(m_numLocalBndCoeffs,-1);
            m_localToGlobalBndSign  = Array<OneD, NekDouble>(m_numLocalBndCoeffs,1.0);
            m_signChange = true;
            m_solnType = solnType;
            m_staticCondLevel = 0;
            m_numPatches = exp1D->size();
            m_numLocalBndCoeffsPerPatch =  Array<OneD, unsigned int>(m_numPatches);
            m_numLocalIntCoeffsPerPatch =  Array<OneD, unsigned int>(m_numPatches);
            for(i = 0; i < m_numPatches; ++i) 
            {
                m_numLocalBndCoeffsPerPatch[i] = (unsigned int) (*exp1D)[i]->NumDGBndryCoeffs();
                m_numLocalIntCoeffsPerPatch[i] = (unsigned int) 0;
            }

            map<int, int> MeshVertToLocalVert;            
            // Order the Dirichlet vertices first.
            gid = 0;
            for(i = 0; i < nbnd; i++)
            {
                if(bndCond[i]->GetBoundaryConditionType() ==SpatialDomains::eDirichlet)
                {
                    m_numDirichletBndPhys++;
                    vid = ((bndCondExp[i])->GetVertex())->GetVid();

                    MeshVertToLocalVert[vid] = gid++;
                }
            }
            
            // set up simple map based on vertex and edge id's
            cnt = 0;
            for(i = 0; i < exp1D->size(); ++i)
            {
                if(locSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>((*exp1D)[i]))
                {
                    locSegExp->GetBoundaryMap(vmap);
                        
                    for(j = 0; j < locSegExp->GetNverts(); ++j)
                    {   
                        vid = (locSegExp->GetGeom1D())->GetVid(j);
                        
                        if(MeshVertToLocalVert.count(vid) == 0)
                        {
                            MeshVertToLocalVert[vid] = gid++;
                        }   
                        
                        m_localToGlobalBndMap[cnt + j] = 
                            MeshVertToLocalVert.find(vid)->second;
                    }    
                    cnt += locSegExp->NumBndryCoeffs();
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a segment expansion failed");
                }
            }  

            // Set up boundary mapping
            m_bndCondCoeffsToGlobalCoeffsMap = Array<OneD, int >(nbnd);
            m_bndCondCoeffsToGlobalCoeffsSign = Array<OneD, NekDouble >(nbnd,1.0);
            m_numLocalDirBndCoeffs = 0;
            m_numDirichletBndPhys = 0;

            for(i = 0; i < nbnd; ++i)
            {
                vid = ((bndCondExp[i])->GetVertex())->GetVid();
                m_bndCondCoeffsToGlobalCoeffsMap[i] = MeshVertToLocalVert.find(vid)->second;
                
                if(bndCond[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                {
                    m_numLocalDirBndCoeffs += 1;
                    m_numDirichletBndPhys  += 1;
                }
            }

            m_numGlobalDirBndCoeffs = m_numLocalDirBndCoeffs;
            CalculateBndSystemBandWidth();
        }

        LocalToGlobalDGMap::LocalToGlobalDGMap(SpatialDomains::MeshGraph2D &graph2D, 
                                               const GenExpList1DSharedPtr &trace, 
                                               const boost::shared_ptr<StdRegions::StdExpansionVector> &exp2D,
                                               const GlobalSysSolnType solnType, 
                                               const Array<OneD, MultiRegions::ExpList1DSharedPtr> &bndCondExp,
                                               const Array<OneD, SpatialDomains::BoundaryConditionShPtr> &bndCond, 
                                               const map<int,int> &periodicEdges)
        {
            int i,j,k,cnt,id, id1, order_e,gid;
            int ntrace_exp = trace->GetExpSize();
            int nel        = exp2D->size();
            int nbnd = bndCondExp.num_elements();
            int offset; 
            LocalRegions::SegExpSharedPtr  locSegExp,locSegExp1;
            LocalRegions::QuadExpSharedPtr locQuadExp;
            LocalRegions::TriExpSharedPtr  locTriExp;
            SpatialDomains::Geometry1DSharedPtr SegGeom;
            
            map<int, int> MeshEdgeId;

            m_signChange = true;
            
            // determine mapping from geometry edges to trace
            for(i = 0; i < ntrace_exp; ++i)
            {
                if(locSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(trace->GetExp(i)))
                {
                    id = (locSegExp->GetGeom1D())->GetEid();
                    
                    if(periodicEdges.count(id) > 0)
                    {
                        if(MeshEdgeId.count(id) == 0)
                        {
                            id1 = periodicEdges.find(id)->second;
                            MeshEdgeId[id] = i;
                            MeshEdgeId[id1] = i;
                        }
                    }
                    else
                    {
                        MeshEdgeId[id] = i;
                    }
                }
                else
                {
                    ASSERTL0(false,"Dynamics cast to segment expansion failed");
                }
            }

            // Count total number of edges 
            cnt = 0;
            for(i = 0; i < nel; ++i)
            {
                cnt += (*exp2D)[i]->GetNedges();
            }
            
            Array<OneD, StdRegions::StdExpansion1DSharedPtr> edgemap(cnt);
            m_elmtToTrace = Array<OneD, Array<OneD,StdRegions::StdExpansion1DSharedPtr> >(nel);

            // set up edge expansions links; 
            cnt = 0;
            for(i = 0; i < nel; ++i)
            {
                m_elmtToTrace[i] = edgemap + cnt; 

                if(locQuadExp = boost::dynamic_pointer_cast<LocalRegions::QuadExp>((*exp2D)[i]))
                {
                    for(j = 0; j < locQuadExp->GetNedges(); ++j)
                    {   
                        SegGeom = (locQuadExp->GetGeom2D())->GetEdge(j);
                        
                        id = SegGeom->GetEid();
                        
                        if(MeshEdgeId.count(id) > 0)
                        {
                            m_elmtToTrace[i][j] = boost::dynamic_pointer_cast< LocalRegions::GenSegExp> ((*trace).GetExp(MeshEdgeId.find(id)->second));

                        }
                        else
                        {
                            ASSERTL0(false,"Failed to find edge map");
                        }
                    }
                }
                else if(locTriExp = boost::dynamic_pointer_cast<LocalRegions::TriExp>((*exp2D)[i]))
                {
                    for(j = 0; j < locTriExp->GetNedges(); ++j)
                    {    
                        SegGeom = (locTriExp->GetGeom2D())->GetEdge(j);

                        id = SegGeom->GetEid();
                        
                        if(MeshEdgeId.count(id) > 0)
                        {
                            m_elmtToTrace[i][j] = boost::dynamic_pointer_cast< LocalRegions::GenSegExp> ((*trace).GetExp((MeshEdgeId.find(id))->second));

                        }
                        else
                        {
                            ASSERTL0(false,"Failed to find edge map");
                        }
                    }
                
                }
                else
                {
                    ASSERTL0(false,"dynamic cast to a local 2D expansion failed");
                }                    
                cnt += (*exp2D)[i]->GetNedges();
            }

            // Set up boundary mapping
            cnt = 0;
            for(i = 0; i < nbnd; ++i)
            {
                cnt += bndCondExp[i]->GetExpSize();
            }

            m_bndCondCoeffsToGlobalCoeffsMap = Array<OneD,int >(cnt);
            m_bndExpAdjacentOrient = Array<OneD, AdjacentTraceOrientation > (cnt);
            m_numLocalDirBndCoeffs = 0;
            m_numDirichletBndPhys   = 0;
            cnt = 0;
            for(i = 0; i < bndCondExp.num_elements(); ++i)
            {
                for(j = 0; j < bndCondExp[i]->GetExpSize(); ++j)
                {
                    if(locSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(bndCondExp[i]->GetExp(j)))
                    {
                        SegGeom = locSegExp->GetGeom1D();
                        
                        
#if OLDMAP
                        id = SegGeom->GetEid();
                        if(MeshEdgeId.count(id) > 0)
                        {
                            m_bndCondCoeffsToGlobalCoeffsMap[cnt+j] = MeshEdgeId.find(id)->second;
                        }
                        else
                        {
                            ASSERTL0(false,"Failed to find edge map");
                        }
#endif
                        // Check to see which way boundary edge is
                        // orientated with respect to connecting
                        // element counter-clockwise convention.

                        SpatialDomains::ElementEdgeVectorSharedPtr con_elmt
                            = graph2D.GetElementsFromEdge(SegGeom);
                        
                        if((boost::dynamic_pointer_cast<SpatialDomains::Geometry2D>((*con_elmt)[0]->m_Element))->GetEorient((*con_elmt)[0]->m_EdgeIndx) == StdRegions::eForwards)
                        {
                            m_bndExpAdjacentOrient[cnt+j] = eAdjacentEdgeIsForwards;
                        }
                        else
                        {
                            m_bndExpAdjacentOrient[cnt+j] = eAdjacentEdgeIsBackwards;
                        }
                        
                    }
                    else
                    {
                        ASSERTL0(false,"dynamic cast to a local Segment expansion failed");
                    }
                    
                    if(bndCond[i]->GetBoundaryConditionType() == SpatialDomains::eDirichlet)
                    {
                        m_numLocalDirBndCoeffs  += locSegExp->GetNcoeffs();
                        m_numDirichletBndPhys   += locSegExp->GetTotPoints();
                    }
                }
                cnt += j;
            }
            
            // Set up integer mapping array and sign change for each
            // degree of freedom + initialise some more data members
            m_solnType = solnType;
            m_staticCondLevel = 0;
            m_numPatches = nel;
            m_numLocalBndCoeffsPerPatch =  Array<OneD, unsigned int>(nel);
            m_numLocalIntCoeffsPerPatch =  Array<OneD, unsigned int>(nel);
            int nbndry = 0;
            for(i = 0; i < nel; ++i) // count number of elements in array
            {
                nbndry += (*exp2D)[i]->NumDGBndryCoeffs();
                m_numLocalBndCoeffsPerPatch[i] = (unsigned int) (*exp2D)[i]->NumDGBndryCoeffs();
                m_numLocalIntCoeffsPerPatch[i] = (unsigned int) 0;
            }

            m_numGlobalDirBndCoeffs = m_numLocalDirBndCoeffs;
            m_numLocalBndCoeffs = nbndry;
            m_numLocalCoeffs = nbndry;
            m_localToGlobalBndMap  = Array<OneD, int > (nbndry);
            m_localToGlobalBndSign = Array<OneD, NekDouble > (nbndry,1);

            // Set up array for potential mesh optimsation
            Array<OneD,int> TraceElmtGid(ntrace_exp,-1);
            int nDir = 0;
            cnt = 0;

            // We are now going to construct a graph of the mesh
            // which can be reordered depending on the type of solver we would
            // like to use.
            typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> BoostGraph;
            typedef boost::graph_traits<BoostGraph>::vertex_descriptor BoostVertex;

            BoostGraph boostGraphObj;
            int trace_id,trace_id1;
            
            // make trace edge renumbering map where first solved
            // edge starts at 0 so we can set up graph.
            for(i = 0; i < ntrace_exp; ++i)
            {
                if(trace->GetCoeff_Offset(i) >= m_numLocalDirBndCoeffs)
                {
                    // Initial put in element ordering (starting
                    // from zero) into TraceElmtGid
                    boost::add_vertex(boostGraphObj);
                    TraceElmtGid[i] = cnt++;             
                }
                else
                {
                    // Use existing offset for Dirichlet edges
                    TraceElmtGid[i] = trace->GetCoeff_Offset(i);
                    nDir++;
                }
            }
            
            // Set up boost Graph
            for(i = 0; i < nel; ++i)
            {
                nbndry += (*exp2D)[i]->NumDGBndryCoeffs();
                
                for(j = 0; j < (*exp2D)[i]->GetNedges(); ++j)
                {   
                    locSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(m_elmtToTrace[i][j]);
                    SegGeom = locSegExp->GetGeom1D();
                    
                    // Add edge to boost graph for non-Dirichlet Boundary 
                    id  = SegGeom->GetEid();
                    trace_id = MeshEdgeId.find(id)->second;
                    if(trace->GetCoeff_Offset(trace_id) >= m_numLocalDirBndCoeffs) 
                    {
                        for(k = j+1; k < (*exp2D)[i]->GetNedges(); ++k)
                        {   
                            locSegExp1 = boost::dynamic_pointer_cast<LocalRegions::SegExp>(m_elmtToTrace[i][k]);
                            SegGeom = locSegExp1->GetGeom1D();
                            
                            id1  = SegGeom->GetEid();
                            trace_id1 = MeshEdgeId.find(id1)->second;
                            if(trace->GetCoeff_Offset(trace_id1)
                               >= m_numLocalDirBndCoeffs)
                            {
                                boost::add_edge( (size_t) TraceElmtGid[trace_id], (size_t) TraceElmtGid[trace_id1], boostGraphObj);                       
                            }
                        }
                    }
                }
            }
            
            
            int nGraphVerts = ntrace_exp-nDir;
            Array<OneD, int> perm(nGraphVerts);
            Array<OneD, int> iperm(nGraphVerts);
            BottomUpSubStructuredGraphSharedPtr bottomUpGraph;
            Array<OneD, int> vwgts(nGraphVerts);
            for(i = 0; i < nGraphVerts; ++i)
            {
                vwgts[i] = trace->GetExp(i+nDir)->GetNcoeffs();   
            }

            if(nGraphVerts)
            {
                switch(solnType)
                {
                case eDirectFullMatrix:
                    {
                        NoReordering(boostGraphObj,perm,iperm);
                    }
                    break;
                case eDirectStaticCond:
                    {
                        CuthillMckeeReordering(boostGraphObj,perm,iperm);
                    }
                    break;
                case eDirectMultiLevelStaticCond: 
                    {
                        MultiLevelBisectionReordering(boostGraphObj,vwgts,perm,iperm,bottomUpGraph);
                    }
                    break;
                default:
                    {
                        ASSERTL0(false,"Unrecognised solution type");
                    }
                }  
            }
 
            // Recast the permutation so that it can be
            // used as a map Form old trace edge ID to new trace
            // edge ID
            cnt = m_numLocalDirBndCoeffs;
            for(i = 0; i < ntrace_exp-nDir; ++i)
            {
                TraceElmtGid[perm[i]+nDir]=cnt;
                cnt += trace->GetExp(perm[i]+nDir)->GetNcoeffs();
            }  

            // Now have trace edges Gid position
            nbndry = cnt = 0;
            for(i = 0; i < nel; ++i)
            {
                nbndry += (*exp2D)[i]->NumDGBndryCoeffs();

                for(j = 0; j < (*exp2D)[i]->GetNedges(); ++j)
                {   
                    locSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(m_elmtToTrace[i][j]);
                    SegGeom = locSegExp->GetGeom1D();
                    
                    id  = SegGeom->GetEid();
                    gid = TraceElmtGid[MeshEdgeId.find(id)->second];
                    
                    
                    order_e = locSegExp->GetNcoeffs();
                    
                    if((*exp2D)[i]->GetEorient(j) == StdRegions::eForwards)
                    {
                        for(k = 0; k < order_e; ++k)
                        {
                            m_localToGlobalBndMap[k+cnt] = gid + k;
                        }
                    }
                    else // backwards orientated
                    {
                        switch(locSegExp->GetBasisType(0))
                        {
                        case LibUtilities::eModified_A:
                            // reverse vertex order
                            m_localToGlobalBndMap[cnt] = gid + 1;
                            m_localToGlobalBndMap[cnt+1] = gid;
                            for(k = 2; k < order_e; ++k)
                            {
                                m_localToGlobalBndMap[k+cnt] = gid + k;
                            }

                            // negate odd modes
                            for(k = 3; k < order_e; k+=2)
                            {
                                m_localToGlobalBndSign[cnt+k] = -1.0;
                            }
                            
                            
                            break;
                        case LibUtilities::eGLL_Lagrange:
                            // reverse  order
                            for(k = 0; k < order_e; ++k)
                            {
                                m_localToGlobalBndMap[cnt+order_e-k-1] = gid + k;
                            }
                            break;
                        default:
                            ASSERTL0(false,"Boundary type not permitted");
                            
                        }
                    }
                    cnt += order_e;
                }
            }

            // set up m_bndCondCoeffsToGlobalCoeffsMap to align with map
            cnt = 0;
            for(i = 0; i < nbnd; ++i)
            {
                cnt += bndCondExp[i]->GetNcoeffs();
            }
            
            m_bndCondCoeffsToGlobalCoeffsMap = Array<OneD,int >(cnt);
            
            for(cnt = i = 0; i < nbnd; ++i)
            {
                for(j = 0; j < bndCondExp[i]->GetExpSize(); ++j)
                {
                    if(locSegExp = boost::dynamic_pointer_cast<LocalRegions::SegExp>(bndCondExp[i]->GetExp(j)))
                    {
                        SegGeom = locSegExp->GetGeom1D();
                        id      = SegGeom->GetEid();
                        gid     = TraceElmtGid[MeshEdgeId.find(id)->second];
                        
                        order_e = locSegExp->GetNcoeffs();

                        // Since boundary information is defined to be
                        // aligned with the geometry just use forward
                        // defintiion for gid's
                        for(k = 0; k < order_e; ++k)
                        {
                            m_bndCondCoeffsToGlobalCoeffsMap[cnt++] = gid + k;
                        }
                    }
                }
            }

            m_numGlobalBndCoeffs = trace->GetNcoeffs();
            m_numGlobalCoeffs = m_numGlobalBndCoeffs;

            CalculateBndSystemBandWidth();

            if( (solnType == eDirectMultiLevelStaticCond) && nGraphVerts )
            {
                if(m_staticCondLevel < (bottomUpGraph->GetNlevels()-1))
                {
                    m_nextLevelLocalToGlobalMap = MemoryManager<LocalToGlobalBaseMap>::
                        AllocateSharedPtr(this,bottomUpGraph);
                }
            }
        }
        
    }


}

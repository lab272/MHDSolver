////////////////////////////////////////////////////////////////////////////////
//
// File: CwipiExchange.cpp
//
// For more information, please see: http://www.nektar.info/
//
// The MIT License
//
// Copyright (c) 2015 Kilian Lackhove
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
// Description: CWIPI Exchange class
//
////////////////////////////////////////////////////////////////////////////////

#include "CwipiExchange.h"

#include <LibUtilities/BasicUtils/PtsField.h>

#include <MultiRegions/ContField1D.h>
#include <MultiRegions/ContField2D.h>
#include <MultiRegions/ContField3D.h>

#include <LibUtilities/Foundations/Interp.h>
#include <LibUtilities/Foundations/PhysGalerkinProject.h>

#include <cwipi.h>

#define RECV_ORDER_DIFF 0

namespace Nektar
{
namespace SolverUtils
{

CwipiCoupling::CwipiCoupling(MultiRegions::ExpListSharedPtr field,
                                   string name, string distAppname, int outputFreq, double geomTol) :
    Coupling(field, name),
    m_distAppname(distAppname),
    m_outputFormat("Ensight Gold"),
    m_outputFormatOption("text"),
    m_outputFreq(outputFreq),
    m_geomTol(geomTol),
    m_coords(NULL),
    m_connecIdx(NULL),
    m_connec(NULL)
{
    cwipi_dump_application_properties();

    //  Init Coupling
    cwipi_solver_type_t solver_type = CWIPI_SOLVER_CELL_VERTEX;

    SpatialDomains::MeshGraphSharedPtr graph = m_evalField->GetGraph();
    int spacedim = graph->GetSpaceDimension();

    SpatialDomains::MeshGraphSharedPtr recvGraph = SpatialDomains::MeshGraph::Read(m_evalField->GetSession());
    recvGraph->SetExpansionsToPolyOrder(RECV_ORDER_DIFF + m_evalField->GetExp(0)->GetBasisNumModes(0));

    // HACK: 6
    // TODO: DeclareCoeffPhysArrays
    m_recvFields = Array<OneD, MultiRegions::ExpListSharedPtr>(6);
    switch (spacedim)
    {
        case 1:
        {
            for (int i = 0; i < 6; ++i)
            {
                m_recvFields[i] = MemoryManager<MultiRegions::ExpList1D>::
                    AllocateSharedPtr(m_evalField->GetSession(), recvGraph);
            }
            break;
        }

        case 2:
        {
            for (int i = 0; i < 6; ++i)
            {
                m_recvFields[i] = MemoryManager<MultiRegions::ExpList2D>::
                    AllocateSharedPtr(m_evalField->GetSession(), recvGraph);
            }
            break;
        }

        case 3:
        {
            for (int i = 0; i < 6; ++i)
            {
                m_recvFields[i] = MemoryManager < MultiRegions::ExpList3D >::
                    AllocateSharedPtr(m_evalField->GetSession(), recvGraph);
            }
            break;
        }

        default:
        {
            ASSERTL0(false, "Expansion dimension not recognised");
            break;
        }
    }

    cwipi_create_coupling(m_name.c_str(),
                          CWIPI_COUPLING_PARALLEL_WITH_PARTITIONING,
                          m_distAppname.c_str(),
                          spacedim,
                          m_geomTol,
                          CWIPI_STATIC_MESH,
                          solver_type,
                          m_outputFreq,
                          m_outputFormat.c_str(),
                          m_outputFormatOption.c_str());

    // get Elements
    SpatialDomains::SegGeomMap   seggeom;
    SpatialDomains::TriGeomMap   trigeom;
    SpatialDomains::QuadGeomMap  quadgeom;
    SpatialDomains::TetGeomMap   tetgeom;
    SpatialDomains::PyrGeomMap   pyrgeom;
    SpatialDomains::PrismGeomMap prismgeom;
    SpatialDomains::HexGeomMap   hexgeom;
    if (spacedim == 1)
    {
        seggeom   = graph->GetAllSegGeoms();
    }
    else if (spacedim == 2)
    {
        trigeom   = graph->GetAllTriGeoms();
        quadgeom  = graph->GetAllQuadGeoms();
    }
    else if (spacedim == 3)
    {
        tetgeom   = graph->GetAllTetGeoms();
        pyrgeom   = graph->GetAllPyrGeoms();
        prismgeom = graph->GetAllPrismGeoms();
        hexgeom   = graph->GetAllHexGeoms();
    };

    int nVerts = graph->GetNvertices();
    int nElts = seggeom.size() + trigeom.size() + quadgeom.size() +
                tetgeom.size() + pyrgeom.size() + prismgeom.size() +
                hexgeom.size();

    // allocate CWIPI arrays
    m_coords = (double *) malloc(sizeof(double) * 3 * nVerts);
    ASSERTL1(m_coords != NULL, "malloc failed for m_coords");
    int tmp = 2 * seggeom.size() + 3 * trigeom.size() + 4 * quadgeom.size() +
              4 * tetgeom.size() + 5 * pyrgeom.size() + 6 * prismgeom.size() +
              8 * hexgeom.size();
    m_connec = (int *) malloc(sizeof(int) * tmp);
    ASSERTL1(m_connec != NULL, "malloc failed for m_connec");
    m_connecIdx = (int *) malloc(sizeof(int) * (nElts + 1));
    ASSERTL1(m_connecIdx != NULL, "malloc failed for m_connecIdx");

    m_connecIdx[0] = 0;
    int coordsPos = 0;
    int connecPos = 0;
    int conidxPos = 0;

    AddElementsToMesh(seggeom,   coordsPos, connecPos, conidxPos);
    AddElementsToMesh(trigeom,   coordsPos, connecPos, conidxPos);
    AddElementsToMesh(quadgeom,  coordsPos, connecPos, conidxPos);
    AddElementsToMesh(tetgeom,   coordsPos, connecPos, conidxPos);
    AddElementsToMesh(pyrgeom,   coordsPos, connecPos, conidxPos);
    AddElementsToMesh(prismgeom, coordsPos, connecPos, conidxPos);
    AddElementsToMesh(hexgeom,   coordsPos, connecPos, conidxPos);


    // output the mesh in tecplot format. If this works, CWIPI will be able
    // to process it, too
    /*
    cout << "VARIABLES = \"X\", \"Y\", \"Z\", \"U\"" <<  endl;
    cout << "ZONE NODES=" <<  nVerts << ",ELEMENTS=" <<  nElts;
    cout << ",DATAPACKING=POINT, ZONETYPE=FETETRAHEDRON" <<  endl;
    for (int i = 0; i < nVerts; ++i)
    {
        cout << m_coords[3*i + 0] << " " << m_coords[3*i + 1] << " ";
        cout << m_coords[3*i + 2] << " " <<  1.0 << endl;
    }
    for (int i = 0; i < nElts; ++i)
    {
        cout << m_connec[i*4 + 0] << " " << m_connec[i*4 + 1] << " ";
        cout << m_connec[i*4 + 2] <<  " " <<  m_connec[i*4 + 3] <<  endl;
    }
    */

    cwipi_define_mesh(m_name.c_str(),
                      nVerts,
                      nElts,
                      m_coords,
                      m_connecIdx,
                      m_connec);

    // define the quadrature points at which we want to receive data
    m_nPoints = m_recvFields[0]->GetTotPoints();
    Array <OneD,  Array<OneD,  NekDouble> > coords(3);
    coords[0] = Array<OneD, NekDouble>(m_nPoints);
    coords[1] = Array<OneD, NekDouble>(m_nPoints);
    coords[2] = Array<OneD, NekDouble>(m_nPoints);
    m_recvFields[0]->GetCoords(coords[0], coords[1], coords[2]);

    m_points = (double *) malloc(sizeof(double) * 3 * m_nPoints);
    ASSERTL1(m_points != NULL, "malloc failed for m_points");

    for (int i = 0; i < m_nPoints; ++i)
    {
        m_points[3 * i + 0] = double(coords[0][i]);

        if (spacedim > 1)
        {
            m_points[3 * i + 1] = double(coords[1][i]);
        }
        else
        {
            m_points[3 * i + 1] = 0.0;
        }

        if (spacedim > 2)
        {
            m_points[3 * i + 2] = double(coords[2][i]);
        }
        else
        {
            m_points[3 * i + 2] = 0.0;
        }
    }

    cwipi_set_points_to_locate(m_name.c_str(), m_nPoints, m_points);
}


CwipiCoupling::~CwipiCoupling()
{
    free(m_coords);
    free(m_points);
    free(m_connec);
    free(m_connecIdx);
}


void CwipiCoupling::v_FinalizeCoupling(void)
{
    cwipi_delete_coupling(m_name.c_str());
}


template <typename T>
void CwipiCoupling::AddElementsToMesh(T geom, int &coordsPos, int &connecPos,
        int &conidxPos)
{
    // helper variables
    Array<OneD, NekDouble>  x(3);
    SpatialDomains::PointGeomSharedPtr vert;
    int vertID;

    int kNverts = T::mapped_type::element_type::kNverts;

    // iterate over all elements
    typename T::iterator it;

    for (it = geom.begin(); it != geom.end(); it++)
    {
        //  iterate over the elements vertices
        for (int j = 0; j < kNverts; ++j)
        {
            vert = it->second->GetVertex(j);
            vertID = vert->GetVid();

            // check if we already stored the vertex
            if (m_vertMap.count(vertID) ==  0)
            {
                //  store the vertex
                vert->GetCoords(x[0], x[1], x[2]);
                m_coords[3 * coordsPos + 0] = double(x[0]);
                m_coords[3 * coordsPos + 1] = double(x[1]);
                m_coords[3 * coordsPos + 2] = double(x[2]);

                // store the vertex position in the m_coords array
                m_vertMap[vertID] = coordsPos;
                coordsPos++;
            }

            m_connec[connecPos] = m_vertMap[vertID] + 1;
            connecPos++;
        }

        m_connecIdx[conidxPos + 1] = m_connecIdx[conidxPos] + kNverts;
        conidxPos++;
    }
}


void CwipiCoupling::PrintProgressbar(const int position, const int goal) const
{
    // print only every 2 percent
    if (int(100 * position / goal) % 2 ==  0)
    {
        cout << "." <<  flush;
    }
}


CwipiExchange::CwipiExchange(SolverUtils::CouplingSharedPointer coupling, string name,
                     int nEVars) :
    Exchange(coupling,  name),
    m_nEVars(nEVars),
    m_rValsInterl(NULL)
{
    int nPoints = m_coupling->GetNPoints();

    m_rValsInterl = (double *) malloc(sizeof(double) * nPoints * m_nEVars);
    ASSERTL1(m_rValsInterl != NULL, "malloc failed for m_rValsInterl");
}


CwipiExchange::~CwipiExchange()
{
    free(m_rValsInterl);
}


void CwipiExchange::v_SendFields(const int step, const NekDouble time,
                                    Array<OneD, Array<OneD, NekDouble> > &field)
{
    ASSERTL0(false, "not implemented yet")
}


void CwipiExchange::v_ReceiveFields(const int step, const NekDouble time,
                                       Array<OneD, Array<OneD, NekDouble> > &field)
{
    static NekDouble lastUdate = -1;
    ASSERTL0(time > lastUdate, "CwipiExchange::v_ReceiveFields called twice in this timestep")
    lastUdate = time;

    int nPoints = m_coupling->GetNPoints();
    ASSERTL1(m_nEVars ==  field.num_elements(), "field size mismatch");

    cout << "receiving fields at i = " << step << ", t = " << time << endl;

    Array<OneD, MultiRegions::ExpListSharedPtr> recvFields = m_coupling->GetRecvFields();
    MultiRegions::ExpListSharedPtr evalField = m_coupling->GetEvalField();

    // TODO: interpolate from evalField to recvField

    int nNotLoc = 0;

    // workaround a bug in cwipi: receiving_field_name should be const char* but is char*
    char recFN[m_recvFieldName.length() + 1];
    strcpy(recFN, m_recvFieldName.c_str());

    cwipi_exchange(m_coupling->GetName().c_str(),
                   m_name.c_str(),
                   m_nEVars,
                   step,
                   time,
                   "",
                   NULL,
                   recFN,
                   m_rValsInterl,
                   &nNotLoc);

    int tmp = -1;
    const int *notLoc = &tmp;
    if (nNotLoc !=  0)
    {
        cout << "WARNING: relocating " << nNotLoc << " of " << nPoints << " points" <<  endl;
        notLoc = cwipi_get_not_located_points(m_coupling->GetName().c_str());
    }

    int locPos = 0;
    int intPos = 0;
    for (int j = 0; j < m_nEVars; ++j)
    {
        locPos = 0;
        intPos = 0;

        for (int i = 0; i < nPoints; ++i)
        {
            // cwipi indices start from 1
            if (notLoc[locPos] - 1 == i)
            {
                // keep the original value of field[j][i]
                locPos++;
            }
            else
            {
//                 recvFields[j]->UpdatePhys()[i] = m_rValsInterl[intPos * m_nEVars + j];
                field[j][i] = m_rValsInterl[intPos * m_nEVars + j];
                intPos++;
            }
        }
    }

    // TODO: copy coeffs from recvFields to field
}


}
}

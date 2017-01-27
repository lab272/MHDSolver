////////////////////////////////////////////////////////////////////////////////
//
//  File: Curvemesh.h
//
//  For more information, please see: http://www.nektar.info/
//
//  The MIT License
//
//  Copyright (c) 2006 Division of Applied Mathematics, Brown University (USA),
//  Department of Aeronautics, Imperial College London (UK), and Scientific
//  Computing and Imaging Institute, University of Utah (USA).
//
//  License for the specific language governing rights and limitations under
//  Permission is hereby granted, free of charge, to any person obtaining a
//  copy of this software and associated documentation files (the "Software"),
//  to deal in the Software without restriction, including without limitation
//  the rights to use, copy, modify, merge, publish, distribute, sublicense,
//  and/or sell copies of the Software, and to permit persons to whom the
//  Software is furnished to do so, subject to the following conditions:
//
//  The above copyright notice and this permission notice shall be included
//  in all copies or substantial portions of the Software.
//
//  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
//  OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
//  THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
//  FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
//  DEALINGS IN THE SOFTWARE.
//
//  Description: object for individual curve meshes.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_MESHUTILS_SURFACEMESHING_CURVEMESH_H
#define NEKTAR_MESHUTILS_SURFACEMESHING_CURVEMESH_H

#include <boost/shared_ptr.hpp>

#include <NekMeshUtils/MeshElements/Mesh.h>
#include <NekMeshUtils/CADSystem/CADVert.h>
#include <NekMeshUtils/CADSystem/CADCurve.h>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

/**
 * @brief class for meshing individual curves (1d meshing)
 */
namespace Nektar
{
namespace NekMeshUtils
{

class CurveMesh
{
public:
    friend class MemoryManager<CurveMesh>;

    /**
     * @brief default constructor
     */
    CurveMesh(int id, MeshSharedPtr m)
        : m_id(id), m_mesh(m)
    {
        m_cadcurve = m_mesh->m_cad->GetCurve(m_id);
    }

    /**
     * @brief alternative constructor with mesh points already created
     */
    CurveMesh(int id, MeshSharedPtr m, std::vector<NodeSharedPtr> n)
        : m_id(id), m_mesh(m), m_meshpoints(n)
    {
        m_cadcurve = m_mesh->m_cad->GetCurve(m_id);
        CreateEdges();
    }

    /**
     * @brief execute meshing
     */
    void Mesh();

    /**
     * @brief get id of first node
     */
    NodeSharedPtr GetFirstPoint()
    {
        return m_meshpoints[0];
    }

    /**
     * @brief get id of last node
     */
    NodeSharedPtr GetLastPoint()
    {
        return m_meshpoints.back();
    }

    /**
     * @brief get list of mesh nodes
     */
    std::vector<NodeSharedPtr> GetMeshPoints()
    {
        return m_meshpoints;
    }

    std::vector<EdgeSharedPtr> GetMeshEdges()
    {
        return m_meshedges;
    }

    /**
     * @brief get the number of points in the curve
     */
    int GetNumPoints()
    {
        return m_meshpoints.size();
    }

    /**
     * @brief get the length of the curve
     */
    NekDouble GetLength()
    {
        return m_curvelength;
    }

    /**
     * @brief get the CAD curve
     */
    CADCurveSharedPtr GetCADCurve()
    {
        return m_cadcurve;
    }

private:
    /**
     * @brief get node spacing sampling function
     */
    void GetSampleFunction();

    /**
     * @brief get node spacing phi function
     */
    void GetPhiFunction();

    /**
     * @brief evaluate paramter ds at curve location s
     */
    NekDouble EvaluateDS(NekDouble s);

    /**
     * @brief evaluate paramter ps at curve location s
     */
    NekDouble EvaluatePS(NekDouble s);

    /**
     * @brief create edges using m_meshpoints and fill m_edgeSet
     */
    void CreateEdges();

    /// CAD curve
    CADCurveSharedPtr m_cadcurve;
    /// length of the curve in real space
    NekDouble m_curvelength;
    /// number of sampling points used in algorithm
    int m_numSamplePoints;
    /// coords of the ends of the parametric curve
    Array<OneD, NekDouble> m_bounds;
    /// array of function ds evaluations
    std::vector<std::vector<NekDouble> > m_dst;
    /// array of function ps evaluations
    std::vector<std::vector<NekDouble> > m_ps;
    /// spacing function evaluation
    NekDouble Ae;
    /// ds
    NekDouble ds;
    /// number of edges to be made in the curve as defined by the spacing
    /// funtion
    int Ne;
    /// paramteric coordiates of the mesh nodes
    std::vector<NekDouble> meshsvalue;
    /// list of mesh edges in the curvemesh
    std::vector<EdgeSharedPtr> m_meshedges;
    /// id of the curvemesh
    int m_id;
    ///
    MeshSharedPtr m_mesh;
    /// ids of the mesh nodes
    std::vector<NodeSharedPtr> m_meshpoints;
};

typedef boost::shared_ptr<CurveMesh> CurveMeshSharedPtr;
}
}

#endif

////////////////////////////////////////////////////////////////////////////////
//
//  File: TetMesh.h
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
//  Description: class for tet meshing
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKTAR_MESHUTILS_TETMESHING_TETMESH_H
#define NEKTAR_MESHUTILS_TETMESHING_TETMESH_H

#include <boost/shared_ptr.hpp>

#include <LibUtilities/BasicUtils/SharedArray.hpp>
#include <LibUtilities/Memory/NekMemoryManager.hpp>

#include <MeshUtils/MeshElements/MeshElements.h>
#include <MeshUtils/CADSystem/CADSystem.h>
#include <MeshUtils/Octree/Octree.h>
#include <MeshUtils/SurfaceMeshing/SurfaceMesh.h>
#include <MeshUtils/MeshElem.hpp>

namespace Nektar{
namespace MeshUtils{

class TetMesh
{
public:
    friend class MemoryManager<TetMesh>;

    /**
     *@brief default constructor
     */
    TetMesh(MeshSharedPtr m, OctreeSharedPtr oct)
                : m_mesh(m), m_octree(oct)
    {
    };

    /**
     *@brief execute tet meshing
     */
    void Mesh();

private:

    /// print stuff to screen?
    MeshSharedPtr m_mesh;
    /// octree object
    OctreeSharedPtr m_octree;
    /// number of tetrahedra
    int numtet;
    /// conncetivity of the tets from the interface
    std::vector<Array<OneD, int> > tetconnect;
    /// list of all nodes
    std::map<int, MeshNodeSharedPtr> Nodes;
    /// list of all surface edges
    std::map<int, MeshEdgeSharedPtr> Edges;
    /// list of all surface triangles
    std::map<int, MeshTriSharedPtr> Tris;
    /// list of all tets
    std::map<int, MeshTetSharedPtr> Tets;
};

typedef boost::shared_ptr<TetMesh> TetMeshSharedPtr;

}
}

#endif

////////////////////////////////////////////////////////////////////////////////
//
//  File: Node.h
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
//  Description: Mesh manipulation objects.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef NEKMESHUTILS_MESHELEMENTS_NODE
#define NEKMESHUTILS_MESHELEMENTS_NODE

#include <NekMeshUtils/NekMeshUtilsDeclspec.h>

#ifdef NEKTAR_USE_MESHGEN
#include <SpatialDomains/PointGeom.h>
#include <NekMeshUtils/CADSystem/CADSystem.h>
#endif

namespace Nektar
{
namespace NekMeshUtils
{
class Node;
typedef boost::shared_ptr<Node> NodeSharedPtr;

/**
 * @brief Represents a point in the domain.
 *
 * Such points may either be element vertices, or simply control
 * points on high-order edges/faces, although this information is not
 * contained within this class.
 */
class Node
{
public:
    /// Create a new node at a specified coordinate.
    NEKMESHUTILS_EXPORT Node(int pId, NekDouble pX, NekDouble pY, NekDouble pZ)
        : m_id(pId), m_x(pX), m_y(pY), m_z(pZ), m_geom()
    {
    }
    /// Copy an existing node.
    // Node(const Node& pSrc)
    //    : m_id(pSrc.m_id), m_x(pSrc.m_x), m_y(pSrc.m_y),
    //      m_z(pSrc.m_z), m_geom() {}
    NEKMESHUTILS_EXPORT Node() : m_id(0), m_x(0.0), m_y(0.0), m_z(0.0), m_geom()
    {
    }
    NEKMESHUTILS_EXPORT ~Node()
    {
    }

    /// Reset the local id;
    NEKMESHUTILS_EXPORT void SetID(int pId)
    {
        m_id = pId;
    }

    /// Get the local id;
    NEKMESHUTILS_EXPORT int GetID(void)
    {
        return m_id;
    }

    /// Define node ordering based on ID.
    NEKMESHUTILS_EXPORT bool operator<(const Node &pSrc)
    {
        return (m_id < pSrc.m_id);
    }
    /// Define node equality based on coordinate.
    NEKMESHUTILS_EXPORT bool operator==(const Node &pSrc)
    {
        return m_x == pSrc.m_x && m_y == pSrc.m_y && m_z == pSrc.m_z;
    }

    NEKMESHUTILS_EXPORT Node operator+(const Node &pSrc) const
    {
        return Node(m_id, m_x + pSrc.m_x, m_y + pSrc.m_y, m_z + pSrc.m_z);
    }

    NEKMESHUTILS_EXPORT Node operator-(const Node &pSrc) const
    {
        return Node(m_id, m_x - pSrc.m_x, m_y - pSrc.m_y, m_z - pSrc.m_z);
    }

    NEKMESHUTILS_EXPORT Node operator*(const Node &pSrc) const
    {
        return Node(m_id, m_x * pSrc.m_x, m_y * pSrc.m_y, m_z * pSrc.m_z);
    }

    NEKMESHUTILS_EXPORT Node operator*(const NekDouble &alpha) const
    {
        return Node(m_id, alpha * m_x, alpha * m_y, alpha * m_z);
    }

    NEKMESHUTILS_EXPORT Node operator/(const NekDouble &alpha) const
    {
        return Node(m_id, m_x / alpha, m_y / alpha, m_z / alpha);
    }

    NEKMESHUTILS_EXPORT void operator+=(const Node &pSrc)
    {
        m_x += pSrc.m_x;
        m_y += pSrc.m_y;
        m_z += pSrc.m_z;
    }

    NEKMESHUTILS_EXPORT void operator*=(const NekDouble &alpha)
    {
        m_x *= alpha;
        m_y *= alpha;
        m_z *= alpha;
    }

    NEKMESHUTILS_EXPORT void operator/=(const NekDouble &alpha)
    {
        m_x /= alpha;
        m_y /= alpha;
        m_z /= alpha;
    }

    NEKMESHUTILS_EXPORT NekDouble abs2() const
    {
        return m_x * m_x + m_y * m_y + m_z * m_z;
    }

    NEKMESHUTILS_EXPORT NekDouble dot(const Node &pSrc) const
    {
        return m_x * pSrc.m_x + m_y * pSrc.m_y + m_z * pSrc.m_z;
    }

    NEKMESHUTILS_EXPORT Node curl(const Node &pSrc) const
    {
        return Node(m_id,
                    m_y * pSrc.m_z - m_z * pSrc.m_y,
                    m_z * pSrc.m_x - m_x * pSrc.m_z,
                    m_x * pSrc.m_y - m_y * pSrc.m_x);
    }

    NEKMESHUTILS_EXPORT Array<OneD, NekDouble> GetLoc()
    {
        Array<OneD, NekDouble> out(3);
        out[0] = m_x;
        out[1] = m_y;
        out[2] = m_z;
        return out;
    }

    /// Generate a %SpatialDomains::PointGeom for this node.
    NEKMESHUTILS_EXPORT SpatialDomains::PointGeomSharedPtr GetGeom(int coordDim)
    {
        SpatialDomains::PointGeomSharedPtr ret =
            MemoryManager<SpatialDomains::PointGeom>::AllocateSharedPtr(
                coordDim, m_id, m_x, m_y, m_z);

        return ret;
    }

    NEKMESHUTILS_EXPORT NekDouble Distance(NodeSharedPtr &p)
    {
        return sqrt((m_x - p->m_x) * (m_x - p->m_x) +
                    (m_y - p->m_y) * (m_y - p->m_y) +
                    (m_z - p->m_z) * (m_z - p->m_z));
    }

    NEKMESHUTILS_EXPORT NekDouble Angle(NodeSharedPtr &a, NodeSharedPtr &b)
    {
        Array<OneD, NekDouble> va(3), vb(3), cn(3);
        va[0] = a->m_x - m_x;
        va[1] = a->m_y - m_y;
        va[2] = a->m_z - m_z;
        vb[0] = b->m_x - m_x;
        vb[1] = b->m_y - m_y;
        vb[2] = b->m_z - m_z;

        NekDouble lva = sqrt(va[0] * va[0] + va[1] * va[1] + va[2] * va[2]);
        NekDouble lvb = sqrt(vb[0] * vb[0] + vb[1] * vb[1] + vb[2] * vb[2]);

        NekDouble aw = 1.0 / (lva * lvb);

        NekDouble cosw = (va[0] * vb[0] + va[1] * vb[1] + va[2] * vb[2]) * aw;

        cn[0] = vb[1] * va[2] - vb[2] * va[1];
        cn[1] = vb[2] * va[0] - vb[0] * va[2];
        cn[2] = vb[0] * va[1] - vb[1] * va[0];

        NekDouble lcn = sqrt(cn[0] * cn[0] + cn[1] * cn[1] + cn[2] * cn[2]);

        NekDouble sinw = aw * lcn;

        NekDouble an = atan2(sinw, cosw);

        if (an < 0)
            an += 6.2831853071796;

        return an;
    }

#ifdef NEKTAR_USE_MESHGEN // fucntions for cad information

    void SetCADCurve(int i, CADCurveSharedPtr c, NekDouble t)
    {
        CADCurveList[i] = std::pair<CADCurveSharedPtr, NekDouble>(c, t);
    }

    void SetCADSurf(int i, CADSurfSharedPtr s, Array<OneD, NekDouble> uv)
    {
        CADSurfList[i] =
            std::pair<CADSurfSharedPtr, Array<OneD, NekDouble> >(s, uv);
    }

    CADCurveSharedPtr GetCADCurve(int i)
    {
        std::map<int, std::pair<CADCurveSharedPtr, NekDouble> >::iterator
            search = CADCurveList.find(i);
        ASSERTL0(search->first == i, "node not on this curve");

        return search->second.first;
    }

    CADSurfSharedPtr GetCADSurf(int i)
    {
        std::map<int, std::pair<CADSurfSharedPtr, Array<OneD, NekDouble> > >::
            iterator search = CADSurfList.find(i);
        ASSERTL0(search->first == i, "surface not found");

        return search->second.first;
    }

    NekDouble GetCADCurveInfo(int i)
    {
        std::map<int, std::pair<CADCurveSharedPtr, NekDouble> >::iterator
            search = CADCurveList.find(i);
        ASSERTL0(search->first == i, "node not on this curve");

        return search->second.second;
    }

    Array<OneD, NekDouble> GetCADSurfInfo(int i)
    {
        std::map<int, std::pair<CADSurfSharedPtr, Array<OneD, NekDouble> > >::
            iterator search = CADSurfList.find(i);
        ASSERTL0(search->first == i, "surface not found");

        return search->second.second;
    }

    std::vector<std::pair<int, CADSurfSharedPtr> > GetCADSurfs()
    {
        std::vector<std::pair<int, CADSurfSharedPtr> > lst;
        std::map<int, std::pair<CADSurfSharedPtr, Array<OneD, NekDouble> > >::
            iterator s;
        for (s = CADSurfList.begin(); s != CADSurfList.end(); s++)
        {
            lst.push_back(
                std::pair<int, CADSurfSharedPtr>(s->first, s->second.first));
        }
        return lst;
    }

    int GetNumCadCurve()
    {
        return CADCurveList.size();
    }

    int GetNumCADSurf()
    {
        return CADSurfList.size();
    }

    void Move(Array<OneD, NekDouble> l, int s, Array<OneD, NekDouble> uv)
    {
        m_x                 = l[0];
        m_y                 = l[1];
        m_z                 = l[2];
        CADSurfSharedPtr su = CADSurfList[s].first;
        CADSurfList[s] =
            std::pair<CADSurfSharedPtr, Array<OneD, NekDouble> >(su, uv);
    }

#endif

    /// ID of node.
    int m_id;
    /// X-coordinate.
    NekDouble m_x;
    /// Y-coordinate.
    NekDouble m_y;
    /// Z-coordinate.
    NekDouble m_z;

#ifdef NEKTAR_USE_MESHGEN // tag to tell the meshelemnets to include cad
                          // information
    /// list of cadcurves the node lies on
    std::map<int, std::pair<CADCurveSharedPtr, NekDouble> > CADCurveList;
    /// list of cadsurfs the node lies on
    std::map<int, std::pair<CADSurfSharedPtr, Array<OneD, NekDouble> > >
        CADSurfList;
#endif

private:
    SpatialDomains::PointGeomSharedPtr m_geom;
};
/// Shared pointer to a Node.

NEKMESHUTILS_EXPORT bool operator==(NodeSharedPtr const &p1,
                                    NodeSharedPtr const &p2);
NEKMESHUTILS_EXPORT bool operator<(NodeSharedPtr const &p1,
                                   NodeSharedPtr const &p2);
NEKMESHUTILS_EXPORT bool operator!=(NodeSharedPtr const &p1,
                                    NodeSharedPtr const &p2);
std::ostream &operator<<(std::ostream &os, const NodeSharedPtr &n);

/**
 * @brief Defines a hash function for nodes.
 *
 * The hash of a node is straight-forward; a combination of the x, y,
 * and z co-ordinates in this order.
 */
struct NodeHash : std::unary_function<NodeSharedPtr, std::size_t>
{
    std::size_t operator()(NodeSharedPtr const &p) const
    {
        std::size_t seed = 0;
        boost::hash_combine(seed, p->m_x);
        boost::hash_combine(seed, p->m_y);
        boost::hash_combine(seed, p->m_z);
        return seed;
    }
};
typedef boost::unordered_set<NodeSharedPtr, NodeHash> NodeSet;
}
}

#endif
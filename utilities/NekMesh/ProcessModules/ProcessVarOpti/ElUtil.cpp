////////////////////////////////////////////////////////////////////////////////
//
//  File: ElUtil.cpp
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
//  Description: Calculate jacobians of elements.
//
////////////////////////////////////////////////////////////////////////////////

#include "ElUtil.h"
#include "ProcessVarOpti.h"

#include <LibUtilities/Foundations/ManagerAccess.h>

using namespace std;

namespace Nektar
{
namespace Utilities
{

boost::mutex mtx2;

ElUtil::ElUtil(ElementSharedPtr e, DerivUtilSharedPtr d,
               ResidualSharedPtr r, int n)
{
    m_el = e;
    derivUtil = d;
    res = r;
    m_mode = n;
    m_dim = m_el->GetDim();
    vector<NodeSharedPtr> ns;
    m_el->GetCurvedNodes(ns);
    nodes.resize(ns.size());
    for (int i = 0; i < ns.size(); ++i)
    {
        nodes[i].resize(m_dim);
        nodes[i][0] = &ns[i]->m_x;

        if (m_dim >= 2)
        {
            nodes[i][1] = &ns[i]->m_y;
        }

        if (m_dim >= 3)
        {
            nodes[i][2] = &ns[i]->m_z;
        }

        m_idmap[ns[i]->m_id] = i;
    }
    maps = MappingIdealToRef();
}

vector<Array<OneD, NekDouble> > ElUtil::MappingIdealToRef()
{
    vector<Array<OneD, NekDouble> > ret;

    LibUtilities::ShapeType st = m_el->GetConf().m_e;

    if(m_el->GetConf().m_e == LibUtilities::eQuadrilateral)
    {
        vector<Array<OneD, NekDouble> > xyz(6);
        vector<NodeSharedPtr> ns = m_el->GetVertexList();
        for(int i = 0; i < 4; i++)
        {
            Array<OneD, NekDouble> x(3);
            x[0] = ns[i]->m_x;
            x[1] = ns[i]->m_y;
            x[2] = ns[i]->m_z;
            xyz[i] = x;
        }

        vector<DNekMat> tmp;

        for (int i = 0; i < derivUtil->ptsHigh; ++i)
        {
            NekDouble a1  = 0.5 * (1 - derivUtil->ptx[i]);
            NekDouble a2  = 0.5 * (1 + derivUtil->ptx[i]);
            NekDouble b1  = 0.5 * (1 - derivUtil->pty[i]);
            NekDouble b2  = 0.5 * (1 + derivUtil->pty[i]);

            DNekMat J(2, 2, 1.0, eFULL);

            J(0,0) = - 0.5 * b1 * xyz[0][0] + 0.5 * b1 * xyz[1][0]
                     + 0.5 * b2 * xyz[2][0] - 0.5 * b2 * xyz[3][0];
            J(1,0) = - 0.5 * b1 * xyz[0][1] + 0.5 * b1 * xyz[1][1]
                     + 0.5 * b2 * xyz[2][1] - 0.5 * b2 * xyz[3][1];

            J(0,1) = - 0.5 * a1 * xyz[0][0] - 0.5 * a2 * xyz[1][0]
                     + 0.5 * a2 * xyz[2][0] + 0.5 * a1 * xyz[3][0];
            J(1,1) = - 0.5 * a1 * xyz[0][1] - 0.5 * a2 * xyz[1][1]
                     + 0.5 * a2 * xyz[2][1] + 0.5 * a1 * xyz[3][1];

            J.Invert();

            tmp.push_back(J);

            Array<OneD, NekDouble> r(10,0.0); //store det in 10th entry

            r[9] = 1.0 / (J(0,0) * J(1,1) - J(0,1) * J(1,0));

            r[0] = J(0,0);
            r[1] = J(1,0);
            r[2] = 0.0;
            r[3] = J(0,1);
            r[4] = J(1,1);
            r[5] = 0.0;
            r[6] = 0.0;
            r[7] = 0.0;
            r[8] = 0.0;
            ret.push_back(r);
        }
    }
    else if(m_el->GetConf().m_e == LibUtilities::eTriangle)
    {
        DNekMat J(2,2,0.0);
        J(0,0) = (*nodes[1][0] - *nodes[0][0]);
        J(1,0) = (*nodes[1][1] - *nodes[0][1]);
        J(0,1) = (*nodes[2][0] - *nodes[0][0]);
        J(1,1) = (*nodes[2][1] - *nodes[0][1]);

        J.Invert();

        DNekMat R(2,2,0.0);
        R(0,0) = 2.0;
        R(1,1) = 2.0;

        J = J * R;

        for(int i = 0 ; i < derivUtil->ptsHigh; i++)
        {
            Array<OneD, NekDouble> r(10,0.0); //store det in 10th entry

            r[9] = 1.0 / (J(0,0) * J(1,1) - J(0,1) * J(1,0));
            r[0] = J(0,0);
            r[1] = J(1,0);
            r[2] = 0.0;
            r[3] = J(0,1);
            r[4] = J(1,1);
            r[5] = 0.0;
            r[6] = 0.0;
            r[7] = 0.0;
            r[8] = 0.0;
            ret.push_back(r);
        }
    }
    else if(m_el->GetConf().m_e == LibUtilities::eTetrahedron)
    {
        DNekMat J(3,3,0.0);
        J(0,0) = (*nodes[1][0] - *nodes[0][0]);
        J(1,0) = (*nodes[1][1] - *nodes[0][1]);
        J(2,0) = (*nodes[1][2] - *nodes[0][2]);
        J(0,1) = (*nodes[2][0] - *nodes[0][0]);
        J(1,1) = (*nodes[2][1] - *nodes[0][1]);
        J(2,1) = (*nodes[2][2] - *nodes[0][2]);
        J(0,2) = (*nodes[3][0] - *nodes[0][0]);
        J(1,2) = (*nodes[3][1] - *nodes[0][1]);
        J(2,2) = (*nodes[3][2] - *nodes[0][2]);

        J.Invert();

        DNekMat R(3,3,0.0);
        R(0,0) = 2.0;
        R(1,1) = 2.0;
        R(2,2) = 2.0;

        J = J * R;

        for(int i = 0 ; i < derivUtil->ptsHigh; i++)
        {
            Array<OneD, NekDouble> r(10,0.0); //store det in 10th entry

            r[9] = 1.0/(J(0,0)*(J(1,1)*J(2,2)-J(2,1)*J(1,2))
                       -J(0,1)*(J(1,0)*J(2,2)-J(2,0)*J(1,2))
                       +J(0,2)*(J(1,0)*J(2,1)-J(2,0)*J(1,1)));

            r[0] = J(0,0);
            r[1] = J(1,0);
            r[2] = J(2,0);
            r[3] = J(0,1);
            r[4] = J(1,1);
            r[5] = J(2,1);
            r[6] = J(0,2);
            r[7] = J(1,2);
            r[8] = J(2,2);
            ret.push_back(r);
        }
    }
    else if(m_el->GetConf().m_e == LibUtilities::ePrism)
    {
        vector<Array<OneD, NekDouble> > xyz(6);
        vector<NodeSharedPtr> ns = m_el->GetVertexList();
        for(int i = 0; i < 6; i++)
        {
            Array<OneD, NekDouble> x(3);
            x[0] = ns[i]->m_x;
            x[1] = ns[i]->m_y;
            x[2] = ns[i]->m_z;
            xyz[i] = x;
        }

        for (int i = 0; i < derivUtil->ptsHigh; ++i)
        {
            NekDouble a1  = 0.5 * (1 - derivUtil->ptx[i]);
            NekDouble a2  = 0.5 * (1 + derivUtil->ptx[i]);
            NekDouble b1  = 0.5 * (1 - derivUtil->pty[i]);
            NekDouble b2  = 0.5 * (1 + derivUtil->pty[i]);
            NekDouble c1  = 0.5 * (1 - derivUtil->ptz[i]);
            NekDouble c2  = 0.5 * (1 + derivUtil->ptz[i]);
            NekDouble d   = 0.5 * (derivUtil->ptx[i] + derivUtil->ptz[i]);

            DNekMat J(3, 3, 1.0, eFULL);

            J(0,0) = - 0.5 * b1 * xyz[0][0] + 0.5 * b1 * xyz[1][0]
                     + 0.5 * b2 * xyz[2][0] - 0.5 * b2 * xyz[3][0];
            J(1,0) = - 0.5 * b1 * xyz[0][1] + 0.5 * b1 * xyz[1][1]
                     + 0.5 * b2 * xyz[2][1] - 0.5 * b2 * xyz[3][1];
            J(2,0) = - 0.5 * b1 * xyz[0][2] + 0.5 * b1 * xyz[1][2]
                     + 0.5 * b2 * xyz[2][2] - 0.5 * b2 * xyz[3][2];

            J(0,1) =   0.5 * d  * xyz[0][0] - 0.5 * a2 * xyz[1][0]
                     + 0.5 * a2 * xyz[2][0] - 0.5 * d  * xyz[3][0]
                     - 0.5 * c2 * xyz[4][0] + 0.5 * c2 * xyz[5][0];
            J(1,1) =   0.5 * d  * xyz[0][1] - 0.5 * a2 * xyz[1][1]
                     + 0.5 * a2 * xyz[2][1] - 0.5 * d  * xyz[3][1]
                     - 0.5 * c2 * xyz[4][1] + 0.5 * c2 * xyz[5][1];
            J(2,1) =   0.5 * d  * xyz[0][2] - 0.5 * a2 * xyz[1][2]
                     + 0.5 * a2 * xyz[2][2] - 0.5 * d  * xyz[3][2]
                     - 0.5 * c2 * xyz[4][2] + 0.5 * c2 * xyz[5][2];

            J(0,2) = - 0.5 * b1 * xyz[0][0] - 0.5 * b2 * xyz[3][0]
                     + 0.5 * b1 * xyz[4][0] + 0.5 * b2 * xyz[5][0];
            J(1,2) = - 0.5 * b1 * xyz[0][1] - 0.5 * b2 * xyz[3][1]
                     + 0.5 * b1 * xyz[4][1] + 0.5 * b2 * xyz[5][1];
            J(2,2) = - 0.5 * b1 * xyz[0][2] - 0.5 * b2 * xyz[3][2]
                     + 0.5 * b1 * xyz[4][2] + 0.5 * b2 * xyz[5][2];

            J.Invert();

            Array<OneD, NekDouble> r(10,0.0); //store det in 10th entry

            r[9] = 1.0/(J(0,0)*(J(1,1)*J(2,2)-J(2,1)*J(1,2))
                       -J(0,1)*(J(1,0)*J(2,2)-J(2,0)*J(1,2))
                       +J(0,2)*(J(1,0)*J(2,1)-J(2,0)*J(1,1)));

            r[0] = J(0,0);
            r[1] = J(1,0);
            r[2] = J(2,0);
            r[3] = J(0,1);
            r[4] = J(1,1);
            r[5] = J(2,1);
            r[6] = J(0,2);
            r[7] = J(1,2);
            r[8] = J(2,2);
            ret.push_back(r);
        }

    }
    else if (m_el->GetConf().m_e == LibUtilities::eHexahedron)
    {
        vector<Array<OneD, NekDouble> > xyz(8);
        vector<NodeSharedPtr> ns = m_el->GetVertexList();
        for(int i = 0; i < 8; i++)
        {
            Array<OneD, NekDouble> x(3);
            x[0] = ns[i]->m_x;
            x[1] = ns[i]->m_y;
            x[2] = ns[i]->m_z;
            xyz[i] = x;
        }

        for (int i = 0; i < derivUtil->ptsHigh; ++i)
        {
            NekDouble a1  = 0.5 * (1 - derivUtil->ptx[i]);
            NekDouble a2  = 0.5 * (1 + derivUtil->ptx[i]);
            NekDouble b1  = 0.5 * (1 - derivUtil->pty[i]);
            NekDouble b2  = 0.5 * (1 + derivUtil->pty[i]);
            NekDouble c1  = 0.5 * (1 - derivUtil->ptz[i]);
            NekDouble c2  = 0.5 * (1 + derivUtil->ptz[i]);

            DNekMat J(3, 3, 1.0, eFULL);

            J(0,0) = - 0.5 * b1 * c1 * xyz[0][0] + 0.5 * b1 * c1 * xyz[1][0]
                     + 0.5 * b2 * c1 * xyz[2][0] - 0.5 * b2 * c1 * xyz[3][0]
                     - 0.5 * b1 * c2 * xyz[5][0] + 0.5 * b1 * c2 * xyz[5][0]
                     + 0.5 * b2 * c2 * xyz[6][0] - 0.5 * b2 * c2 * xyz[7][0];
            J(1,0) = - 0.5 * b1 * c1 * xyz[0][1] + 0.5 * b1 * c1 * xyz[1][1]
                     + 0.5 * b2 * c1 * xyz[2][1] - 0.5 * b2 * c1 * xyz[3][1]
                     - 0.5 * b1 * c2 * xyz[5][1] + 0.5 * b1 * c2 * xyz[5][1]
                     + 0.5 * b2 * c2 * xyz[6][1] - 0.5 * b2 * c2 * xyz[7][1];
            J(2,0) = - 0.5 * b1 * c1 * xyz[0][2] + 0.5 * b1 * c1 * xyz[1][2]
                     + 0.5 * b2 * c1 * xyz[2][2] - 0.5 * b2 * c1 * xyz[3][2]
                     - 0.5 * b1 * c2 * xyz[5][2] + 0.5 * b1 * c2 * xyz[5][2]
                     + 0.5 * b2 * c2 * xyz[6][2] - 0.5 * b2 * c2 * xyz[7][2];

            J(0,1) = - 0.5 * a1 * c1 * xyz[0][0] - 0.5 * a2 * c1 * xyz[1][0]
                     + 0.5 * a2 * c1 * xyz[2][0] + 0.5 * a1 * c1 * xyz[3][0]
                     - 0.5 * a1 * c2 * xyz[5][0] - 0.5 * a2 * c2 * xyz[5][0]
                     + 0.5 * a2 * c2 * xyz[6][0] + 0.5 * a1 * c2 * xyz[7][0];
            J(1,1) = - 0.5 * a1 * c1 * xyz[0][1] - 0.5 * a2 * c1 * xyz[1][1]
                     + 0.5 * a2 * c1 * xyz[2][1] + 0.5 * a1 * c1 * xyz[3][1]
                     - 0.5 * a1 * c2 * xyz[5][1] - 0.5 * a2 * c2 * xyz[5][1]
                     + 0.5 * a2 * c2 * xyz[6][1] + 0.5 * a1 * c2 * xyz[7][1];
            J(2,1) = - 0.5 * a1 * c1 * xyz[0][2] - 0.5 * a2 * c1 * xyz[1][2]
                     + 0.5 * a2 * c1 * xyz[2][2] + 0.5 * a1 * c1 * xyz[3][2]
                     - 0.5 * a1 * c2 * xyz[5][2] - 0.5 * a2 * c2 * xyz[5][2]
                     + 0.5 * a2 * c2 * xyz[6][2] + 0.5 * a1 * c2 * xyz[7][2];

            J(0,0) = - 0.5 * b1 * a1 * xyz[0][0] - 0.5 * b1 * a2 * xyz[1][0]
                     - 0.5 * b2 * a2 * xyz[2][0] - 0.5 * b2 * a1 * xyz[3][0]
                     + 0.5 * b1 * a1 * xyz[5][0] + 0.5 * b1 * a2 * xyz[5][0]
                     + 0.5 * b2 * a2 * xyz[6][0] + 0.5 * b2 * a1 * xyz[7][0];
            J(1,0) = - 0.5 * b1 * a1 * xyz[0][1] - 0.5 * b1 * a2 * xyz[1][1]
                     - 0.5 * b2 * a2 * xyz[2][1] - 0.5 * b2 * a1 * xyz[3][1]
                     + 0.5 * b1 * a1 * xyz[5][1] + 0.5 * b1 * a2 * xyz[5][1]
                     + 0.5 * b2 * a2 * xyz[6][1] + 0.5 * b2 * a1 * xyz[7][1];
            J(2,0) = - 0.5 * b1 * a1 * xyz[0][2] - 0.5 * b1 * a2 * xyz[1][2]
                     - 0.5 * b2 * a2 * xyz[2][2] - 0.5 * b2 * a1 * xyz[3][2]
                     + 0.5 * b1 * a1 * xyz[5][2] + 0.5 * b1 * a2 * xyz[5][2]
                     + 0.5 * b2 * a2 * xyz[6][2] + 0.5 * b2 * a1 * xyz[7][2];

            J.Invert();

            Array<OneD, NekDouble> r(10,0.0); //store det in 10th entry

            r[9] = 1.0/(J(0,0)*(J(1,1)*J(2,2)-J(2,1)*J(1,2))
                       -J(0,1)*(J(1,0)*J(2,2)-J(2,0)*J(1,2))
                       +J(0,2)*(J(1,0)*J(2,1)-J(2,0)*J(1,1)));

            r[0] = J(0,0);
            r[1] = J(1,0);
            r[2] = J(2,0);
            r[3] = J(0,1);
            r[4] = J(1,1);
            r[5] = J(2,1);
            r[6] = J(0,2);
            r[7] = J(1,2);
            r[8] = J(2,2);
            ret.push_back(r);
        }
    }
    else
    {
        ASSERTL0(false,"not coded");
    }

    return ret;
}

void ElUtil::Evaluate()
{
    NekDouble mx = -1.0 * numeric_limits<double>::max();
    NekDouble mn =  numeric_limits<double>::max();

    ASSERTL0(nodes.size() == derivUtil->ptsLow,"node count wrong");

    if(m_dim == 2)
    {
        NekVector<NekDouble> X(nodes.size()),Y(nodes.size());
        for(int j = 0; j < nodes.size(); j++)
        {
            X(j) = *nodes[j][0];
            Y(j) = *nodes[j][1];
        }

        NekVector<NekDouble> x1i(nodes.size()),y1i(nodes.size()),
                             x2i(nodes.size()),y2i(nodes.size());

        x1i = derivUtil->VdmDL[0]*X;
        y1i = derivUtil->VdmDL[0]*Y;
        x2i = derivUtil->VdmDL[1]*X;
        y2i = derivUtil->VdmDL[1]*Y;

        for(int j = 0; j < nodes.size(); j++)
        {
            NekDouble jacDet = x1i(j) * y2i(j) - x2i(j)*y1i(j);
            mx = max(mx,jacDet);
            mn = min(mn,jacDet);
        }
    }
    else if(m_dim == 3)
    {
        NekVector<NekDouble> X(nodes.size()),Y(nodes.size()),Z(nodes.size());
        for(int j = 0; j < nodes.size(); j++)
        {
            X(j) = *nodes[j][0];
            Y(j) = *nodes[j][1];
            Z(j) = *nodes[j][2];
        }

        NekVector<NekDouble> x1i(nodes.size()),y1i(nodes.size()),z1i(nodes.size()),
                             x2i(nodes.size()),y2i(nodes.size()),z2i(nodes.size()),
                             x3i(nodes.size()),y3i(nodes.size()),z3i(nodes.size());

        x1i = derivUtil->VdmDL[0]*X;
        y1i = derivUtil->VdmDL[0]*Y;
        z1i = derivUtil->VdmDL[0]*Z;
        x2i = derivUtil->VdmDL[1]*X;
        y2i = derivUtil->VdmDL[1]*Y;
        z2i = derivUtil->VdmDL[1]*Z;
        x3i = derivUtil->VdmDL[2]*X;
        y3i = derivUtil->VdmDL[2]*Y;
        z3i = derivUtil->VdmDL[2]*Z;

        for(int j = 0; j < nodes.size(); j++)
        {
            DNekMat dxdz(3,3,1.0,eFULL);
            dxdz(0,0) = x1i(j);
            dxdz(0,1) = x2i(j);
            dxdz(0,2) = x3i(j);
            dxdz(1,0) = y1i(j);
            dxdz(1,1) = y2i(j);
            dxdz(1,2) = y3i(j);
            dxdz(2,0) = z1i(j);
            dxdz(2,1) = z2i(j);
            dxdz(2,2) = z3i(j);

            NekDouble jacDet = dxdz(0,0)*(dxdz(1,1)*dxdz(2,2)-dxdz(2,1)*dxdz(1,2))
                              -dxdz(0,1)*(dxdz(1,0)*dxdz(2,2)-dxdz(2,0)*dxdz(1,2))
                              +dxdz(0,2)*(dxdz(1,0)*dxdz(2,1)-dxdz(2,0)*dxdz(1,1));

            mx = max(mx,jacDet);
            mn = min(mn,jacDet);
        }

        /*NekVector<NekDouble> x1i2(derivUtil->ptsHigh),y1i2(derivUtil->ptsHigh),z1i2(derivUtil->ptsHigh),
                             x2i2(derivUtil->ptsHigh),y2i2(derivUtil->ptsHigh),z2i2(derivUtil->ptsHigh),
                             x3i2(derivUtil->ptsHigh),y3i2(derivUtil->ptsHigh),z3i2(derivUtil->ptsHigh);

        x1i2 = derivUtil->VdmD[0]*X;
        y1i2 = derivUtil->VdmD[0]*Y;
        z1i2 = derivUtil->VdmD[0]*Z;
        x2i2 = derivUtil->VdmD[1]*X;
        y2i2 = derivUtil->VdmD[1]*Y;
        z2i2 = derivUtil->VdmD[1]*Z;
        x3i2 = derivUtil->VdmD[2]*X;
        y3i2 = derivUtil->VdmD[2]*Y;
        z3i2 = derivUtil->VdmD[2]*Z;

        NekDouble mx2 = -1.0 * numeric_limits<double>::max();
        NekDouble mn2 =  numeric_limits<double>::max();

        for(int j = 0; j < derivUtil->ptsHigh; j++)
        {
            DNekMat dxdz(3,3,1.0,eFULL);
            dxdz(0,0) = x1i2(j);
            dxdz(0,1) = x2i2(j);
            dxdz(0,2) = x3i2(j);
            dxdz(1,0) = y1i2(j);
            dxdz(1,1) = y2i2(j);
            dxdz(1,2) = y3i2(j);
            dxdz(2,0) = z1i2(j);
            dxdz(2,1) = z2i2(j);
            dxdz(2,2) = z3i2(j);

            NekDouble jacDet = dxdz(0,0)*(dxdz(1,1)*dxdz(2,2)-dxdz(2,1)*dxdz(1,2))
                              -dxdz(0,1)*(dxdz(1,0)*dxdz(2,2)-dxdz(2,0)*dxdz(1,2))
                              +dxdz(0,2)*(dxdz(1,0)*dxdz(2,1)-dxdz(2,0)*dxdz(1,1));

            mx2 = max(mx2,jacDet);
            mn2 = min(mn2,jacDet);
        }

        cout << mn/mx << " " << mn2/mx2 << endl;*/
    }

    mtx2.lock();
    if(mn < 0)
    {
        res->startInv++;
    }
    res->worstJac = min(res->worstJac,mn/mx);
    mtx2.unlock();

    //maps = MappingIdealToRef();

    minJac = mn;
    scaledJac = mn/mx;
}

ElUtilJob* ElUtil::GetJob()
{
    return new ElUtilJob(this);
}

}
}

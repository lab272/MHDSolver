////////////////////////////////////////////////////////////////////////////////
//
//  File: OptimseFunctions.cpp
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
//  Description: functions and gradients for optimisation
//
////////////////////////////////////////////////////////////////////////////////

#include <NekMeshUtils/SurfaceMeshing/SurfaceMesh.h>

using namespace std;
namespace Nektar
{
namespace NekMeshUtils
{

NekDouble Dot(Array<OneD, NekDouble> a, Array<OneD, NekDouble> b)
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}
Array<OneD, NekDouble> Take(Array<OneD, NekDouble> a, Array<OneD, NekDouble> b)
{
    Array<OneD, NekDouble> ret(3);
    ret[0] = a[0] - b[0];
    ret[1] = a[1] - b[1];
    ret[2] = a[2] - b[2];
    return ret;
}
Array<OneD, NekDouble> Times(NekDouble t, Array<OneD, NekDouble> a)
{
    Array<OneD, NekDouble> ret(3);
    ret[0] = a[0] *t;
    ret[1] = a[1] *t;
    ret[2] = a[2] *t;
    return ret;
}
Array<OneD, NekDouble> Add(Array<OneD, NekDouble> a, Array<OneD, NekDouble> b)
{
    Array<OneD, NekDouble> ret(3);
    ret[0] = a[0] + b[0];
    ret[1] = a[1] + b[1];
    ret[2] = a[2] + b[2];
    return ret;
}

void SurfaceMesh::CurveEdgeJac(Array<OneD, NekDouble> t, Array<OneD, NekDouble> z,
                                  CADCurveSharedPtr c, DNekMat &Jac, DNekMat &Hes)
{
    vector<Array<OneD, NekDouble> > r;
    vector<Array<OneD, NekDouble> > dr;
    vector<Array<OneD, NekDouble> > d2r;
    for(int i = 0; i < t.num_elements(); i++)
    {
        Array<OneD, NekDouble> ri(3), dri(3), d2ri(3);
        Array<OneD, NekDouble> d2 = c->D2(t[i]);
        for(int j = 0; j < 3; j++)
        {
            ri[j] = d2[j];
            dri[j] = d2[j+3];
            d2ri[j] = d2[j+6];
        }
        r.push_back(ri);
        dr.push_back(dri);
        d2r.push_back(d2ri);
    }

    DNekMat J(t.num_elements() - 2 , 1, 0.0);
    for(int i = 0; i < t.num_elements() - 2; i++)
    {
        J(i,0) = 2.0/(z[i+1] - z[i]) * Dot(dr[i+1],Take(r[i+1],r[i])) +
                 2.0/(z[i+2] - z[i+1]) * Dot(dr[i+1],Take(r[i+1],r[i+2]));
    }

    DNekMat H(t.num_elements() - 2, t.num_elements() - 2, 0.0);
    for(int i = 0; i < t.num_elements() - 2; i++)
    {
        for(int j = 0; j < t.num_elements() - 2; j++)
        {
            if(i == j)
            {
                H(i,j) = 2.0/(z[i+1] - z[i]) * (Dot(dr[i+1], dr[i+1]) +
                                                Dot(d2r[i+1], Take(r[i+1],r[i]))) +
                         2.0/(z[i+2] - z[i+1]) * (Dot(dr[i+1], dr[i+1]) +
                                                         Dot(d2r[i+1], Take(r[i+1],r[i+2])));
            }
            else if(abs(i-j) == 1)
            {
                int k = max(i,j);
                H(i,j) = -2.0/(z[k+1] - z[k]) * Dot(dr[k], dr[k+1]);
            }
        }
    }
    Jac = J;
    Hes = H;
}

void SurfaceMesh::FaceEdgeJac(Array<OneD, Array<OneD, NekDouble> > uv,
                              Array<OneD, NekDouble> z, CADSurfSharedPtr s,
                              DNekMat &Jac, DNekMat &Hes)
{
    vector<Array<OneD, NekDouble> > r;
    vector<Array<OneD, NekDouble> > dru, drv;
    vector<Array<OneD, NekDouble> > d2ru, d2ruv, d2rv;
    for(int i = 0; i < uv.num_elements(); i++)
    {
        Array<OneD, NekDouble> ri(3), drui(3), drvi(3), d2rui(3), d2rvi(3), d2ruvi(3);
        Array<OneD, NekDouble> d2 = s->D2(uv[i]);
        for(int j = 0; j < 3; j++)
        {
            ri[j] = d2[j];
            drui[j] = d2[j+3];
            drvi[j] = d2[j+6];
            d2rui[j] = d2[j+9];
            d2rvi[j] = d2[j+12];
            d2ruvi[j] = d2[j+15];
        }
        r.push_back(ri);
        dru.push_back(drui);
        drv.push_back(drvi);
        d2ru.push_back(d2rui);
        d2rv.push_back(d2rvi);
        d2ruv.push_back(d2ruvi);
    }

    DNekMat J(2*(uv.num_elements() - 2), 1, 0.0);
    for(int i = 0; i < uv.num_elements() - 2; i++)
    {
        J(2*i+0,0) = 2.0/(z[i+1]-z[i]) * Dot(dru[i+1], Take(r[i+1], r[i])) +
                     2.0/(z[i+2]-z[i+1]) * Dot(dru[i+1], Take(r[i+1], r[i+2]));

        J(2*i+1,0) = 2.0/(z[i+1]-z[i]) * Dot(drv[i+1], Take(r[i+1], r[i])) +
                     2.0/(z[i+2]-z[i+1]) * Dot(drv[i+1], Take(r[i+1], r[i+2]));
    }

    DNekMat H(2*(uv.num_elements() - 2), 2*(uv.num_elements() - 2), 0.0);
    for(int i = 0; i < uv.num_elements() - 2; i++)
    {
        for(int j = 0; j < uv.num_elements() - 2; j++)
        {
            if(i == j)
            {
                H(2*i+0,2*j+0) = 2.0/(z[i+1]-z[i]) * (Dot(d2ru[i+1],Take(r[i+1],r[i])) + Dot(dru[i+1],dru[i+1])) +
                                 2.0/(z[i+2]-z[i+1]) * (Dot(d2ru[i+1],Take(r[i+1],r[i+2])) + Dot(dru[i+1],dru[i+1]));

                H(2*i+1,2*j+1) = 2.0/(z[i+1]-z[i]) * (Dot(d2rv[i+1],Take(r[i+1],r[i])) + Dot(drv[i+1],drv[i+1])) +
                                 2.0/(z[i+2]-z[i+1]) * (Dot(d2rv[i+1],Take(r[i+1],r[i+2])) + Dot(drv[i+1],drv[i+1]));

                H(2*i+0,2*j+1) = 2.0/(z[i+1]-z[i]) * (Dot(d2ruv[i+1],Take(r[i+1],r[i])) + Dot(dru[i+1],drv[i+1])) +
                                 2.0/(z[i+2]-z[i+1]) * (Dot(d2ruv[i+1],Take(r[i+1],r[i+2])) + Dot(dru[i+1],drv[i+1]));

                H(2*i+1,2*j+0) = H(2*i+0,2*j+1);
            }
            else if(i-j == 1)
            {
                //matrix is symmetric about i==j therefore i and j are interchangable and only need to be calculated for the i > j case
                H(2*i+0,2*j+0) = -2.0/(z[i+1]-z[i]) * Dot(dru[i], dru[i+1]);
                H(2*i+1,2*j+1) = -2.0/(z[i+1]-z[i]) * Dot(drv[i], drv[i+1]);
                H(2*j+0,2*i+0) = H(2*i+0,2*j+0);
                H(2*j+1,2*i+1) = H(2*i+1,2*j+1);

                H(2*i+0,2*j+1) = -2.0/(z[i+1]-z[i])*Dot(dru[i+1], drv[i]);
                H(2*j+1,2*i+0) = H(2*i+0,2*j+1);

                H(2*i+1,2*j+0) = -2.0/(z[i+1]-z[i])*Dot(drv[i+1], dru[i]);
                H(2*j+0,2*i+1) = H(2*i+1,2*j+0);
            }
        }
    }

    Jac = J;
    Hes = H;

}

void SurfaceMesh::FaceFaceJac(Array<OneD, Array<OneD, NekDouble> > uv,
                              map<int, vector<NekDouble> > z, map<int, vector<int> > n, CADSurfSharedPtr s,
                              DNekMat &Jac, DNekMat &Hes)
{
    int np = uv.num_elements();
    int nq = 0.5*(sqrt(1+8*np)-1);
    int numint = (nq-3)*(nq-2)/2;

    DNekMat J(numint*2,1,0.0);
    DNekMat H(numint*2,numint*2,0.0);

    vector<Array<OneD, NekDouble> > r;
    vector<Array<OneD, NekDouble> > dru, drv;
    vector<Array<OneD, NekDouble> > d2ru, d2ruv, d2rv;
    for(int i = 0; i < uv.num_elements(); i++)
    {
        Array<OneD, NekDouble> ri(3), drui(3), drvi(3), d2rui(3), d2rvi(3), d2ruvi(3);
        Array<OneD, NekDouble> d2 = s->D2(uv[i]);
        for(int j = 0; j < 3; j++)
        {
            ri[j] = d2[j];
            drui[j] = d2[j+3];
            drvi[j] = d2[j+6];
            d2rui[j] = d2[j+9];
            d2rvi[j] = d2[j+12];
            d2ruvi[j] = d2[j+15];
        }
        r.push_back(ri);
        dru.push_back(drui);
        drv.push_back(drvi);
        d2ru.push_back(d2rui);
        d2rv.push_back(d2rvi);
        d2ruv.push_back(d2ruvi);
    }

    // index the derivatives for each point such that 0 is du1 1 is dv1 2 is du2 and so on.
    map<int, vector<Array<OneD, NekDouble> > > deriv;
    for(int i = 0; i < np - numint; i++)
    {
        vector<Array<OneD, NekDouble> > der;
        for(int j = 0; j < numint; j++)
        {
            Array<OneD, NekDouble> u(3,0.0), v(3,0.0);
            der.push_back(u);
            der.push_back(v);
        }
        deriv[i] = der;
    }
    for(int i = 0; i < numint; i++)
    {
        vector<Array<OneD, NekDouble> > der;
        for(int j = 0; j < numint; j++)
        {
            Array<OneD, NekDouble> u(3,0.0), v(3,0.0);
            if(i == j)
            {
                u = dru[i+numint];
                v = drv[i+numint];
            }
            der.push_back(u);
            der.push_back(v);
        }
        deriv[i+np-numint] = der;
    }

    for(int i = 0; i < numint; i++)
    {
        for(int j = 0; j < numint; j++)
        {
            for(int k = 0; k < 3; k++)
            {
                Array<OneD, NekDouble> drui = Times(0.5, Add(deriv[n[i+np-numint][k*2+0]][2*j+0], deriv[n[i+np-numint][k*2+1]][2*j+0]) );
                Array<OneD, NekDouble> drvi = Times(0.5, Add(deriv[n[i+np-numint][k*2+0]][2*j+1], deriv[n[i+np-numint][k*2+1]][2*j+1]) );
                Array<OneD, NekDouble> ri  = Times(0.5, Add(r[n[i+np-numint][k*2+0]], r[n[i+np-numint][k*2+1]]) );
                J(2*i+0,0) = z[i+np-numint][k*2+0]*Dot(Take(deriv[i+np-numint][2*j+0],drui),Take(r[i+np-numint], ri));
                J(2*i+1,0) = z[i+np-numint][k*2+0]*Dot(Take(deriv[i+np-numint][2*j+1],drvi),Take(r[i+np-numint], ri));
            }
        }
    }



    Jac = J;
    Hes = H;
}

}
}

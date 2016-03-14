////////////////////////////////////////////////////////////////////////////////
//
//  File: SurfaceMeshing.cpp
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
//  Description: surfacemeshing object methods.
//
////////////////////////////////////////////////////////////////////////////////

#include <NekMeshUtils/Optimisation/BGFS-B.h>
#include <LibUtilities/LinearAlgebra/NekMatrix.hpp>

#include <limits>
#include <set>

using namespace std;
namespace Nektar
{
namespace NekMeshUtils
{

// this function will perform an update on the solution vector contained within
// opti
bool BGFSUpdate(OptiObjSharedPtr opti, DNekMat &J, DNekMat &B, DNekMat &H)
{

    Array<OneD, NekDouble> xi = opti->Getxi();
    Array<OneD, NekDouble> ui = opti->Getui();
    Array<OneD, NekDouble> li = opti->Getli();

    set<int> Fset;
    Array<OneD, NekDouble> ti(xi.num_elements());
    for (int i = 0; i < ti.num_elements(); i++)
    {
        if (J(i, 0) < 0)
        {
            ti[i] = (xi[i] - ui[i]) / J(i, 0);
        }
        else if (J(i, 0) > 0)
        {
            ti[i] = (xi[i] - li[i]) / J(i, 0);
        }
        else
        {
            ti[i] = numeric_limits<double>::max();
        }
        if (ti[i] > 0)
        {
            Fset.insert(i);
        }
    }

    // intitialise d
    DNekMat d(xi.num_elements(), 1);
    for (int i = 0; i < xi.num_elements(); i++)
    {
        if (fabs(ti[i]) < 1E-10)
        {
            d(i, 0) = 0;
        }
        else
        {
            d(i, 0) = -1.0 * J(i, 0);
        }
    }

    Array<OneD, NekDouble> xci(xi.num_elements());

    bool hitbounded = false;

    for (int i = 0; i < xci.num_elements(); i++)
    {
        if (xi[i] + d(i, 0) < li[i])
        {
            xci[i] = li[i];
            Fset.erase(i);
            cout << "hit bounded" << endl;
            hitbounded = true;
            continue;
        }
        else
        {
            xci[i] = xi[i] + d(i, 0);
        }

        if (xi[i] + d(i, 0) > ui[i])
        {
            xci[i] = ui[i];
            Fset.erase(i);
            cout << "hit bounded" << endl;
            hitbounded = true;
            continue;
        }
        else
        {
            xci[i] = xi[i] + d(i, 0);
        }
    }

    DNekMat Z(xci.num_elements(), xci.num_elements(), 0.0);

    set<int>::iterator it;
    for (int i = 0; i < xci.num_elements(); i++)
    {
        it = Fset.find(i);
        if (it != Fset.end())
        {
            Z(i, i) = 1.0;
        }
    }

    DNekMat dx(xci.num_elements(), 1, 0.0);
    for (int i = 0; i < xci.num_elements(); i++)
    {
        dx(i, 0) = xci[i] - xi[i];
    }

    DNekMat rg = Z * (J + B * dx);

    DNekMat du = -1.0 * H * rg;

    NekDouble alpha = 1.0;
    for (it = Fset.begin(); it != Fset.end(); it++)
    {
        int i = (*it);
        if (li[i] - xci[i] > alpha * du(i, 0))
        {
            alpha = min(alpha, (li[i] - xci[i]) / du(i, 0));
        }
        else if (ui[i] - xci[i] < alpha * du(i, 0))
        {
            alpha = min(alpha, (ui[i] - xci[i]) / du(i, 0));
        }
    }

    DNekMat grad = alpha * du;

    Array<OneD, NekDouble> dk(xci.num_elements()), xibar(xci.num_elements());
    for (int i = 0; i < xci.num_elements(); i++)
    {
        set<int>::iterator f = Fset.find(i);
        if (f != Fset.end())
        {
            xibar[i] = xci[i] + grad(i, 0);
        }
        else
        {
            xibar[i] = xci[i];
        }
    }

    /*if(hitbounded)
    {
        cout << endl << endl << Z << endl << endl;
        cout << rg << endl << endl;
        cout << du << endl << endl;
        cout << alpha << endl << endl;
        cout << grad << endl << endl;

        for(int i = 0; i < xi.num_elements(); i++)
        {
            cout << xci[i] << " " << xibar[i] << endl;
        }

        exit(-1);
    }*/

    Vmath::Vsub(xci.num_elements(), &xibar[0], 1, &xi[0], 1, &dk[0], 1);

    NekDouble c = 0.0;
    NekDouble r = 0.0;
    NekDouble l = 0.0;
    for (int i = 0; i < dk.num_elements(); i++)
    {
        c += 1E-4 * J(i, 0) * dk[i];
        r += J(i, 0) * dk[i];
    }

    /*cout << endl << J << endl << endl;

    for(int i = 0; i < dk.num_elements(); i++)
    {
        cout << dk[i] << endl;
    }
    cout << endl;*/

    // this section needs a case evaluation for edges on faces
    NekDouble lam = 2.0;
    int iterct    = 0;
    NekDouble fo  = opti->F(xi);
    NekDouble fn;
    Array<OneD, NekDouble> tst(xi.num_elements());
    do
    {
        if (iterct > 20)
        {
            // cout << "failed line search" << endl;
            return false;
        }
        iterct++;

        lam *= 0.5;

        for (int i = 0; i < xi.num_elements(); i++)
        {
            tst[i] = xi[i] + lam * dk[i];
        }

        fn = opti->F(tst);

        DNekMat jn = opti->dF(tst);

        l = 0.0;
        for (int i = 0; i < dk.num_elements(); i++)
        {
            l += jn(i, 0) * dk[i];
        }

    } while (fn > fo + c || fabs(l) > 0.9 * fabs(r));
    // wolfe conditions

    // tst at this point is the new all vector
    // now need to update hessians
    DNekMat Jn = opti->dF(tst);
    DNekMat y  = Jn - J;
    DNekMat yT = Jn - J;
    yT.Transpose();
    DNekMat s(dk.num_elements(), 1, 0.0);
    for (int i = 0; i < dk.num_elements(); i++)
    {
        s(i, 0) = lam * dk[i];
    }
    DNekMat sT = s;
    sT.Transpose();

    DNekMat d1 = yT * s;
    DNekMat d2 = sT * B * s;
    DNekMat d3 = sT * y;
    DNekMat n1 = yT * H * y;

    NekDouble ynorm = 0.0;
    for (int i = 0; i < dk.num_elements(); i++)
    {
        ynorm += y(i, 0) * y(i, 0);
    }

    if (d3(0, 0) > 2.2E-16 * ynorm)
    {
        B = B + y * yT * (1.0 / d1(0, 0)) - B * s * sT * B * (1.0 / d2(0, 0));
        H = H + (d3(0, 0) + n1(0, 0)) / d3(0, 0) / d3(0, 0) * s * sT -
            1.0 / d3(0, 0) * (H * y * sT + s * yT * H);
    }

    J = Jn;
    opti->Update(tst);

    return true;
}
}
}
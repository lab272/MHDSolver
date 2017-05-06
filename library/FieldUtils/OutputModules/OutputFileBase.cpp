////////////////////////////////////////////////////////////////////////////////
//
//  File: OutputFileBase.cpp
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
//  Description: Base class for outputting to a file
//
////////////////////////////////////////////////////////////////////////////////

#include <set>
#include <string>
using namespace std;

#include "OutputFileBase.h"
#include <LibUtilities/BasicUtils/FileSystem.h>
#include <boost/format.hpp>
#include <iomanip>

namespace Nektar
{
namespace FieldUtils
{

OutputFileBase::OutputFileBase(FieldSharedPtr f) : OutputModule(f)
{
    m_requireEquiSpaced = false;
}

OutputFileBase::~OutputFileBase()
{
}

void OutputFileBase::Process(po::variables_map &vm)
{
    string filename = m_config["outfile"].as<string>();

    if(m_f->m_fieldPts != LibUtilities::NullPtsField)
    {
        ASSERTL0(!m_f->m_writeBndFld,
                "Boundary can't be obtained from pts.");
        if( WriteFile(filename, vm))
        {
            OutputFromPts(vm);
        }
    }
    else if(m_f->m_exp.size())
    {
        // reset expansion definition to use equispaced points if required.
        if (m_requireEquiSpaced && (vm.count("noequispaced") == 0 ) &&
            m_f->m_exp[0]->GetNumElmts() != 0)
        {
            // Information to create new expansion
            int numFields   = m_f->m_exp.size();
            m_f->m_fielddef = m_f->m_exp[0]->GetFieldDefinitions();

            // Set points to equispaced
            int nPointsNew  = 0;
            if (vm.count("output-points"))
            {
                nPointsNew = vm["output-points"].as<int>();
            }
            m_f->m_graph->SetExpansionsToEvenlySpacedPoints(nPointsNew);

            // Save original expansion
            vector<MultiRegions::ExpListSharedPtr> expOld = m_f->m_exp;
            // Create new expansion
            m_f->m_exp[0] = m_f->SetUpFirstExpList(m_f->m_numHomogeneousDir,
                                            true);
            for(int i = 1; i < numFields; ++i)
            {
                m_f->m_exp[i] =
                        m_f->AppendExpList(m_f->m_numHomogeneousDir);
            }
            // Extract result to new expansion
            for(int i = 0; i < numFields; ++i)
            {
                m_f->m_exp[i]->ExtractCoeffsToCoeffs(
                        expOld[i],
                        expOld[i]->GetCoeffs(),
                        m_f->m_exp[i]->UpdateCoeffs());
                m_f->m_exp[i]->BwdTrans(
                        m_f->m_exp[i]->GetCoeffs(),
                        m_f->m_exp[i]->UpdatePhys());
            }
            // Extract boundary expansion if needed
            if (m_f->m_writeBndFld)
            {
                Array<OneD, const MultiRegions::ExpListSharedPtr> BndExpOld;
                MultiRegions::ExpListSharedPtr                    BndExp;
                for(int i = 0; i < numFields; ++i)
                {
                    BndExpOld = expOld[i]->GetBndCondExpansions();
                    for(int j = 0; j < BndExpOld.num_elements(); ++j)
                    {
                        BndExp = m_f->m_exp[i]->UpdateBndCondExpansion(j);

                        BndExp->ExtractCoeffsToCoeffs(
                            BndExpOld[j],
                            BndExpOld[j]->GetCoeffs(),
                            BndExp->UpdateCoeffs());
                    }
                }
            }
        }

        if (m_f->m_writeBndFld)
        {
            if (m_f->m_verbose && m_f->m_comm->TreatAsRankZero())
            {
                cout << "\t" << GetModuleName()
                     << ": Writing boundary file(s): ";
                for (int i = 0; i < m_f->m_bndRegionsToWrite.size(); ++i)
                {
                    cout << m_f->m_bndRegionsToWrite[i];
                    if (i < m_f->m_bndRegionsToWrite.size() - 1)
                    {
                        cout << ", ";
                    }
                }
                cout << endl;
            }
        }
        if (m_f->m_writeBndFld && m_f->m_exp[0]->GetNumElmts() != 0)
        {
            int nfields = m_f->m_exp.size();
            int normdim = m_f->m_graph->GetMeshDimension();

            // Prepare for normals output
            if (m_f->m_addNormals)
            {
                // Prepare for creating expansions for normals
                m_f->m_exp.resize(nfields + normdim);;

                // Include normal name in m_variables
                string normstr[3] = {"Norm_x", "Norm_y", "Norm_z"};
                for (int j = 0; j < normdim; ++j)
                {
                    m_f->m_exp[nfields+j] =
                        m_f->AppendExpList(m_f->m_numHomogeneousDir);
                    m_f->m_variables.push_back(normstr[j]);
                }
            }

            // Move m_exp to a new expansion
            vector<MultiRegions::ExpListSharedPtr> exp(m_f->m_exp.size());
            exp.swap(m_f->m_exp);

            // Extract data to boundaryconditions
            if (m_f->m_fldToBnd)
            {
                for (int i = 0; i < exp.size(); ++i)
                {
                    exp[i]->FillBndCondFromField();
                }
            }

            Array<OneD, Array<OneD, const MultiRegions::ExpListSharedPtr> >
                BndExp(exp.size());
            for (int i = 0; i < exp.size(); ++i)
            {
                BndExp[i] = exp[i]->GetBndCondExpansions();
            }

            // get hold of partition boundary regions so we can match it to 
            // desired region extraction
            SpatialDomains::BoundaryConditions bcs(m_f->m_session,
                                                   exp[0]->GetGraph());
            const SpatialDomains::BoundaryRegionCollection bregions =
                bcs.GetBoundaryRegions();
            SpatialDomains::BoundaryRegionCollection::const_iterator breg_it;
            map<int, int> BndRegionMap;
            int cnt = 0;
            for (breg_it = bregions.begin(); breg_it != bregions.end();
                 ++breg_it, ++cnt)
            {
                BndRegionMap[breg_it->first] = cnt;
            }

            // find ending of output file and insert _b1, _b2
            int dot     = filename.find_last_of('.') + 1;
            string ext  = filename.substr(dot, filename.length() - dot);
            string name = filename.substr(0, dot - 1);

            for (int i = 0; i < m_f->m_bndRegionsToWrite.size(); ++i)
            {
                string outname =
                    name + "_b" +
                    boost::lexical_cast<string>(m_f->m_bndRegionsToWrite[i]) + "." +
                    ext;

                if (BndRegionMap.count(m_f->m_bndRegionsToWrite[i]) == 1)
                {
                    if(!WriteFile(outname, vm))
                    {
                        continue;
                    }
                    RegisterConfig("outfile", outname);

                    int Border = BndRegionMap[m_f->m_bndRegionsToWrite[i]];

                    for (int j = 0; j < exp.size(); ++j)
                    {
                        m_f->m_exp[j] = BndExp[j][Border];
                        m_f->m_exp[j]->BwdTrans(
                                    m_f->m_exp[j]->GetCoeffs(),
                                    m_f->m_exp[j]->UpdatePhys());
                    }

                    if (m_f->m_addNormals)
                    {
                        // Get normals
                        Array<OneD, Array<OneD, NekDouble> > NormPhys;
                        exp[0]->GetBoundaryNormals(Border, NormPhys);

                        // add normal coefficients to expansions
                        for (int j = 0; j < normdim; ++j)
                        {
                            m_f->m_exp[nfields+j] = BndExp[nfields+j][Border];
                            Vmath::Vcopy(m_f->m_exp[nfields+j]->GetTotPoints(),
                                NormPhys[j], 1,
                                m_f->m_exp[nfields+j]->UpdatePhys(), 1);
                            m_f->m_exp[nfields+j]->FwdTrans_IterPerExp(
                                    m_f->m_exp[nfields+j]->GetPhys(),
                                    m_f->m_exp[nfields+j]->UpdateCoeffs());
                        }
                    }
                    OutputFromExp(vm);
                }
            }
        }
        else
        {
            if( WriteFile(filename, vm))
            {
                OutputFromExp(vm);
            }
        }
    }
    else if(m_f->m_data.size())
    {
        ASSERTL0(!m_f->m_writeBndFld,
                "Boundary extraction requires xml file.");
        if( WriteFile(filename, vm))
        {
            OutputFromData(vm);
        }
    }
}

bool OutputFileBase::WriteFile(std::string &filename, po::variables_map &vm)
{
    fs::path outFile(filename);
    int writeFile = 1;
    if (fs::exists(outFile) && (vm.count("forceoutput") == 0))
    {
        LibUtilities::CommSharedPtr comm;
        if (m_f->m_comm)
        {
            comm = m_f->m_comm;
        }
        else
        {
            comm = LibUtilities::GetCommFactory().CreateInstance(
                "Serial", 0, 0);
        }

        writeFile = 0; // set to zero for reduce all to be correct.

        if (comm->TreatAsRankZero())
        {
            string answer;
            cout << "Did you wish to overwrite " << filename << " (y/n)? ";
            getline(cin, answer);
            if (answer.compare("y") == 0)
            {
                writeFile = 1;
            }
            else
            {
                cout << "Not writing file " << filename
                     << " because it already exists" << endl;
            }
        }

        comm->AllReduce(writeFile, LibUtilities::ReduceSum);
    }
    return (writeFile == 0) ? false : true;
}

}
}

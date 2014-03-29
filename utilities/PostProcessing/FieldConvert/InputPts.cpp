////////////////////////////////////////////////////////////////////////////////
//
//  File: InputPts.cpp
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
//  Description: Read xml file of a series of points and hold
//
////////////////////////////////////////////////////////////////////////////////

#include <string>
#include <iostream>
using namespace std;

#include <tinyxml/tinyxml.h>

#include "InputPts.h"

namespace Nektar
{
    namespace Utilities
    {
        ModuleKey InputPts::m_className[5] = {
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eInputModule, "pts"), InputPts::create,
                "Reads Pts file."),
            GetModuleFactory().RegisterCreatorFunction(
                ModuleKey(eInputModule, "pts.gz"), InputPts::create,
                "Reads Pts file."),
        };

        /**
         * @brief Set up InputPts object.
         *
         */
        InputPts::InputPts(FieldSharedPtr f) : InputModule(f)
        {
            m_allowedFiles.insert("pts");
            f->m_fieldPts = MemoryManager<FieldPts>::AllocateSharedPtr();
        }
        
        InputPts::~InputPts()
        {
        }

        /**
         *
         */
        void InputPts::Process(po::variables_map &vm)
        {
            if(m_f->m_verbose)
            {
                cout << "Processing input pts file" << endl;
            }

            string pts_ending = "pts";

            TiXmlDocument docInput;
            if (!docInput.LoadFile((m_f->m_inputfiles["pts"][0]).c_str()))
            {
                std::cerr << "Unable to open file '" << m_f->m_inputfiles["pts"][0] << "'." << std::endl;
                exit(1);
            }

            TiXmlElement *nektar = docInput.FirstChildElement("NEKTAR");

            TiXmlElement *points = nektar->FirstChildElement("POINTS");
            
            TiXmlAttribute *dimAttr = points->FirstAttribute();
            std::string attrName(dimAttr->Name());
            
            int dim; 
            int err = dimAttr->QueryIntValue(&dim);
            ASSERTL0(err == TIXML_SUCCESS, "Unable to read attribute DIM.");
            
            m_f->m_fieldPts->m_ptsDim = dim;
            m_f->m_fieldPts->m_pts = Array<OneD, Array<OneD, NekDouble> >(dim);            
      
            TiXmlNode *pointsBody = points->FirstChild();
            
            std::istringstream pointsDataStrm(pointsBody->ToText()->Value());

            vector<NekDouble> pts; 
            NekDouble       in_pts; 
            try
            {
                while(!pointsDataStrm.fail())
                {
                    pointsDataStrm >> in_pts; 

                    pts.push_back(in_pts);
                }
            }
            catch(...)
            {
                ASSERTL0(false, "Unable to read Points data.");
            }

            int npts = pts.size()/dim;
            if(m_f->m_verbose)
            {
                cout << " Read " << npts << " points of dimension " << dim << endl;
            }
            
            for(int i = 0; i < dim; ++i)
            {
                m_f->m_fieldPts->m_pts[i] = Array<OneD, NekDouble>(npts);
            }
            
            for(int i = 0; i < npts; ++i)
            {
                for(int j = 0; j < dim; ++j)
                {
                    m_f->m_fieldPts->m_pts[j][i] = pts[i*dim +j];
                }
            }
        }
    }
}

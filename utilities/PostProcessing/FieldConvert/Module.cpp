////////////////////////////////////////////////////////////////////////////////
//
//  File: Module.cpp
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
//  Description: Abstract input/output modules.
//
////////////////////////////////////////////////////////////////////////////////

#include <iomanip>

#include "Module.h"

using namespace std;

namespace Nektar
{
    namespace Utilities
    {
        /**
         * Returns an instance of the module factory, held as a singleton.
         */
        ModuleFactory& GetModuleFactory()
        {
            typedef Loki::SingletonHolder<ModuleFactory,
                Loki::CreateUsingNew,
                Loki::NoDestroy > Type;
            return Type::Instance();
        }

        /**
         * Prints a given module key to a stream.
         */
        std::ostream& operator<<(std::ostream& os, const ModuleKey& rhs)
        {
            return os << ModuleTypeMap[rhs.first] << ": " << rhs.second;
        }

        InputModule::InputModule(FieldSharedPtr m) : Module(m)
        {
            config["infile"] = ConfigOption(false, "", "Input filename.");
        }
        
        OutputModule::OutputModule(FieldSharedPtr m) : Module(m)
        {
            config["outfile"] = ConfigOption(false, "", "Output filename.");
        }

        /**
         * @brief Open a file for input.
         */
        void InputModule::OpenStream()
        {
            string fname = config["infile"].as<string>();
            mshFile.open(fname.c_str());
            if (!mshFile.good())
            {
                cerr << "Error opening file: " << fname << endl;
                abort();
            }
        }

        /**
         * @brief Open a file for output.
         */
        void OutputModule::OpenStream()
        {
            string fname = config["outfile"].as<string>();
            mshFile.open(fname.c_str());
            if (!mshFile.good())
            {
                cerr << "Error opening file: " << fname << endl;
                abort();
            }
        }
        
        /**
         * @brief Register a configuration option with a module.
         */
        void Module::RegisterConfig(string key, string val)
        {
            map<string, ConfigOption>::iterator it = config.find(key);
            if (it == config.end())
            {
                cerr << "WARNING: Unrecognised config option " << key
                     << ", proceeding anyway." << endl;
            }

            it->second.beenSet = true;
            
            if (it->second.isBool)
            {
                it->second.value = "1";
            }
            else
            {
                it->second.value = val;
            }
        }
        
        /**
         * @brief Print out all configuration options for a module.
         */
        void Module::PrintConfig()
        {
            map<string, ConfigOption>::iterator it;
            
            if (config.size() == 0)
            {
                cerr << "No configuration options for this module." << endl;
                return;
            }
            
            for (it = config.begin(); it != config.end(); ++it)
            {
                cerr << setw(10) << it->first << ": " << it->second.desc 
                     << endl;
            }
        }
        
        /**
         * @brief Sets default configuration options for those which have not
         * been set.
         */
        void Module::SetDefaults()
        {
            map<string, ConfigOption>::iterator it;
            
            for (it = config.begin(); it != config.end(); ++it)
            {
                if (!it->second.beenSet)
                {
                    it->second.value = it->second.defValue;
                }
            }
        }

        /**
         * @brief Print a brief summary of information.
         */
        void InputModule::PrintSummary()
        {
            cout << "Field size = " << f->data[0].num_elements() * sizeof(NekDouble) << endl;
        }
    }
}

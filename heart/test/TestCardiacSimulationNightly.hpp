/*

Copyright (c) 2005-2013, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef TESTCARDIACSIMULATIONNIGHTLY_HPP_
#define TESTCARDIACSIMULATIONNIGHTLY_HPP_

#include <cxxtest/TestSuite.h>

#include "CardiacSimulation.hpp"
#include "PetscSetupAndFinalize.hpp"

#include "CompareHdf5ResultsFiles.hpp"

class TestCardiacSimulationNightly : public CxxTest::TestSuite
{
public:
    void TestCardiacSimulationBasicBidomain() throw(Exception)
    {
        // run a bidomain simulation
        CardiacSimulation simulation("heart/test/data/xml/base_bidomain.xml");
        std::string foldername = "BaseBidomainNightly/";

        // compare the files, using the CompareFilesViaHdf5DataReader() method
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "base_bidomain_results", false,
                                                 foldername, "SimulationResults", true, 1e-9));
    }

    void TestCardiacSimulationBasicMonodomain() throw(Exception)
    {
        // run a bidomain simulation
        CardiacSimulation simulation("heart/test/data/xml/base_monodomain.xml");
        std::string foldername = "BaseMonodomainNightly/";

        // compare the files, using the CompareFilesViaHdf5DataReader() method
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "base_monodomain_results", false,
                                                 foldername, "SimulationResults", true, 1e-9));
    }

    void TestCardiacSimulationPostprocessMonodomain() throw(Exception)
    {
        // run a bidomain simulation
        CardiacSimulation simulation("heart/test/data/xml/postprocess_monodomain.xml");
        std::string foldername = "PostprocessMonodomainNightly/";

        // compare the files, using the CompareFilesViaHdf5DataReader() method
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "postprocess_monodomain_results", false,
                                                 foldername, "SimulationResults", true, 1e-9));
    }

    void TestCardiacSimulationSaveBidomain() throw(Exception)
    {
        // run a bidomain simulation
        CardiacSimulation simulation("heart/test/data/xml/save_bidomain.xml");
        std::string foldername = "SaveBidomainNightly";

        // compare the files, using the CompareFilesViaHdf5DataReader() method
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "save_bidomain_results", false,
                                                 foldername, "SimulationResults", true, 1e-9));

        FileFinder file(foldername + "_checkpoints/10ms/" + foldername + "_10ms/archive.arch.0",
                        RelativeTo::ChasteTestOutput);
        TS_ASSERT(file.Exists());
    }

    void TestCardiacSimulationSaveMonodomain() throw(Exception)
    {
        // run a bidomain simulation
        CardiacSimulation simulation("heart/test/data/xml/save_monodomain.xml");
        std::string foldername = "SaveMonodomainNightly";

        // compare the files, using the CompareFilesViaHdf5DataReader() method
        TS_ASSERT( CompareFilesViaHdf5DataReader("heart/test/data/cardiac_simulations", "save_monodomain_results", false,
                                                 foldername, "SimulationResults", true,
                                                 1e-9 /*tol*/));

        FileFinder file(foldername + "_checkpoints/10ms/" + foldername + "_10ms/archive.arch.0",
                        RelativeTo::ChasteTestOutput);
        TS_ASSERT(file.Exists());
    }
};

#endif /*TESTCARDIACSIMULATIONNIGHTLY_HPP_*/

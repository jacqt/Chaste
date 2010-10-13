/*

Copyright (C) University of Oxford, 2005-2010

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/


#include "AbstractHdf5Converter.hpp"
#include "HeartConfig.hpp"
#include "Version.hpp"


template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractHdf5Converter<ELEMENT_DIM, SPACE_DIM>::AbstractHdf5Converter(std::string inputDirectory,
                          std::string fileBaseName,
                          AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM> *pMesh,
                          std::string subdirectoryName) :
                    mFileBaseName(fileBaseName),
                    mpMesh(pMesh)
{
    // store directory, mesh and filenames and create the reader
    this->mpReader = new Hdf5DataReader(inputDirectory, this->mFileBaseName);
    //Create new directory in which to store everything
    mpOutputFileHandler = new OutputFileHandler(HeartConfig::Instance()->GetOutputDirectory() + "/" + subdirectoryName, false);
    // check the data file read has one, two or three variables
    std::vector<std::string> variable_names = this->mpReader->GetVariableNames();
    mNumVariables = variable_names.size();
    if(mNumVariables==0 || mNumVariables>3)
    {
        delete mpReader;
        delete mpOutputFileHandler;
        EXCEPTION("Data has zero or more than three variables - doesn't appear to be mono, bidomain or extended bidomain");//see #1499 and #1502
    }

    if (mpReader->GetNumberOfRows() != mpMesh->GetNumNodes())
    {
        delete mpReader;
        delete mpOutputFileHandler;
        EXCEPTION("Mesh and HDF5 file have a different number of nodes");
    }

    //Write an info file
    if (PetscTools::AmMaster())
    {
        //Note that we don't want the child processes to write info files

        out_stream p_file = this->mpOutputFileHandler->OpenOutputFile(this->mFileBaseName + "_times.info");
        unsigned num_timesteps = this->mpReader->GetUnlimitedDimensionValues().size();
        *p_file << "Number of timesteps " << num_timesteps << std::endl;
        *p_file << "timestep " << HeartConfig::Instance()->GetPrintingTimeStep() << std::endl;
        double first_timestep=this->mpReader->GetUnlimitedDimensionValues().front();
        *p_file << "First timestep " << first_timestep << std::endl;
        double last_timestep=this->mpReader->GetUnlimitedDimensionValues().back();
        *p_file << "Last timestep " << last_timestep << std::endl;
		*p_file << ChasteBuildInfo::GetProvenanceString();
        p_file->close();

    }
    //Write the parameters out
    HeartConfig::Instance()->Write(true, false, subdirectoryName);

}
template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractHdf5Converter<ELEMENT_DIM,SPACE_DIM>::~AbstractHdf5Converter()
{
    delete mpReader;
    delete mpOutputFileHandler;
}

/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class AbstractHdf5Converter<1,1>;
template class AbstractHdf5Converter<1,2>;
template class AbstractHdf5Converter<2,2>;
template class AbstractHdf5Converter<1,3>;
template class AbstractHdf5Converter<2,3>;
template class AbstractHdf5Converter<3,3>;


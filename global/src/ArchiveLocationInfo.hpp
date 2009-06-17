/*

Copyright (C) University of Oxford, 2005-2009

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
#ifndef ARCHIVELOCATIONINFO_HPP_
#define ARCHIVELOCATIONINFO_HPP_

#include <string>

/**
 * Mini-class to help with 'archiving' various classes that don't write their
 * data directly to the archive file.  They thus need to know information
 * about where the archive is being written to, in order to write their own
 * files into the same folder.  The main methods are GetArchiveDirectory and
 * SetArchiveDirectory.
 * 
 * This functionality is used by the meshes, LinearSystem and HeartConfig.
 * 
 * For the benefit of the meshes (and the cancer code), there are also
 * shortcut methods SetMeshPathname and GetMeshPathname, allowing you to
 * specify the base file name for the mesh.  This is needed because the cancer
 * code adds timestamp information to the file name.
 */
class ArchiveLocationInfo
{
private:

    /** Directory that archives are being written to. */
    static std::string mDirPath;
    
    /** Mesh filename (relative to #mDirPath). */
    static std::string mMeshFilename;
    
public:
    /**
     * Get the base path for the mesh files (minus file extension).
     */
    static std::string GetMeshPathname()
    {
        return mDirPath + mMeshFilename;
    }
    
    /**
     * Set the location to write mesh files.
     * @param rDirectory  the directory to write to.
     * @param rFilename  the base name (minus extension) for the mesh files.
     */
    static void SetMeshPathname(const std::string &rDirectory, const std::string &rFilename)
    {
        SetArchiveDirectory(rDirectory);
        mMeshFilename = rFilename;
    }
    
    /**
     * Get the directory that archives are being written to.
     * Will always end in a '/'.
     */
    static std::string GetArchiveDirectory()
    {
        return mDirPath;
    }
    
    /**
     * Set the directory that archives are being written to.
     * @param the directory in question.
     */
    static void SetArchiveDirectory(const std::string &rDirectory)
    {
        mDirPath = rDirectory;
        if (! ( *(mDirPath.end()-1) == '/'))
        {
            mDirPath = mDirPath + "/";
        }
    }
};

#endif /*ARCHIVELOCATIONINFO_HPP_*/

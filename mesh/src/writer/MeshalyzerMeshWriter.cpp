#ifndef _MESHALYZERMESHWRITER_CPP_
#define _MESHALYZERMESHWRITER_CPP_

#include "MeshalyzerMeshWriter.hpp"

template<int ELEMENT_DIM, int SPACE_DIM>
MeshalyzerMeshWriter<ELEMENT_DIM, SPACE_DIM>::MeshalyzerMeshWriter(const std::string &rDirectory, 
                                           const std::string &rBaseName, 
                                           const bool &rSetCoolGraphics)
    : AbstractMeshWriter<ELEMENT_DIM, SPACE_DIM>(rDirectory, rBaseName)
{
	if ( ELEMENT_DIM != 3 || SPACE_DIM !=3 )
    {
        throw Exception("MeshalyzerMeshWriter(): Can only write 3D data");
    }
    
    if (rSetCoolGraphics) 
    {
		this->mIndexFromZero=false;
		this->mWriteMetaFile=true;
	} 
    else 
    {
		this->mIndexFromZero=true;
		this->mWriteMetaFile=false;
	}
}

template<int ELEMENT_DIM, int SPACE_DIM>
void MeshalyzerMeshWriter<ELEMENT_DIM, SPACE_DIM>::WriteFiles()
{
	//std::string comment="#Generated by Chaste mesh file writer";
	
    //Write node file
	std::string node_file_name = this->mBaseName+".pts";
    out_stream p_node_file = this->mpOutputFileHandler->OpenOutputFile(node_file_name);
	
	//Write the node header
	int num_nodes = this->GetNumNodes();
	*p_node_file<< num_nodes << "\n";
	
	//Write each node's data
	for (int item_num=0; item_num<num_nodes; item_num++)
	{
		std::vector<double> current_item = this->mNodeData[item_num];
		for (unsigned int i=0;i<SPACE_DIM;i++)
		{
			*p_node_file<<current_item[i]<<"\t";
		}
		*p_node_file<<"\n";
		
	}
	p_node_file->close();
	
    //Write Element file
	std::string element_file_name = this->mBaseName+".tetras";
	out_stream p_element_file = this->mpOutputFileHandler->OpenOutputFile(element_file_name);

	//Write the element header
	int num_elements = this->GetNumElements();
	
	*p_element_file<< num_elements << "\n";
		
	//Write each element's data
	int nodes_per_element = 4;
	for (int item_num=0; item_num<num_elements; item_num++)
	{
		std::vector<int> current_item = this->mElementData[item_num];
		for (int i=0;i<nodes_per_element;i++)
		{
			if (this->mIndexFromZero)
			{
				*p_element_file<<current_item[i]<<"\t";
			} 
			else
			{
				*p_element_file<<current_item[i]+1<<"\t";
			}
		}
		*p_element_file<<"\n";
		
	}
	p_element_file->close();

	//Write boundary face file
	std::string face_file_name = this->mBaseName+".tris";
    out_stream p_face_file = this->mpOutputFileHandler->OpenOutputFile(face_file_name);

	//Write the boundary face header
	int num_faces = this->GetNumBoundaryFaces();
	
	*p_face_file<< num_faces << "\n";
		
	//Write each face's data
	double material_property= 0.0;
	for (int item_num=0; item_num<num_faces; item_num++)
	{
		std::vector<int> current_item = this->mBoundaryFaceData[item_num];
		for (unsigned int i=0;i<ELEMENT_DIM;i++)
		{
			if (this->mIndexFromZero)
			{
				*p_face_file<<current_item[i]<<"\t";
			} 
			else
			{
				*p_face_file<<current_item[i]+1<<"\t";
			}
		}
		*p_face_file<<material_property<<"\n";
		
	}
	p_face_file->close();

	if (this->mWriteMetaFile)
    {
		std::string meta_file_name = this->mBaseName+".cg_in";
        out_stream p_meta_file = this->mpOutputFileHandler->OpenOutputFile(meta_file_name);
        
		*p_meta_file<< "1\n"<< "0\n";
		*p_meta_file<< face_file_name<<"\n";
		p_meta_file->close();
	}
}

template<int ELEMENT_DIM, int SPACE_DIM>
MeshalyzerMeshWriter<ELEMENT_DIM, SPACE_DIM>::~MeshalyzerMeshWriter()
{
}

#endif //_MESHALYZERMESHWRITER_CPP_

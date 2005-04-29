/** 
 * Concrete version of the AbstractVisual class.
 * A MatlabVisualizer takes the base name and modify a little of 
 * the .node file and .out file so that all the comments are in 
 * Matlab style.
 */

#ifndef _MATLABVISUALIZER_HPP_
#define _MATLABVISUALIZER_HPP_

#include "AbstractVisualizer.hpp"

class MatlabVisualizer: public AbstractVisualizer
{		
public:
	MatlabVisualizer(std::string pathBaseName, int dimension);
	virtual ~MatlabVisualizer();
	
	void CreateNodesFileForVisualization();	     
	void CreateOutputFileForVisualization();
	std::vector<std::string> GetRawDataFromFile(std::string fileName);
};

#endif //_MATLABVISUALIZER_HPP_

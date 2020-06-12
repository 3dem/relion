#ifndef PDB_HELPER_H
#define PDB_HELPER_H

#include <src/jaz/gravis/t3Vector.h>
#include <vector>
#include <map>
#include <string>
#include <src/assembly.h>


class PdbHelper
{
	public:
		
		static std::map<std::string,std::vector<gravis::d3Vector>>
			splitAtoms(
				const std::string& pdbFile,
				bool byElement);
		
		static std::string getElement(const std::string& atomName);
		static gravis::d3Vector getPosition(const Atom& atom);
		
};

#endif

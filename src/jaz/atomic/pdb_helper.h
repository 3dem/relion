#ifndef PDB_HELPER_H
#define PDB_HELPER_H

#include <src/jaz/gravis/t3Vector.h>
#include <src/assembly.h>
#include <vector>
#include <map>
#include <string>
#include <set>


class PdbHelper
{
	public:

		static std::map<std::string,std::vector<gravis::d3Vector>>
			groupAtomsByName(
				const Assembly& assembly);

		static std::map<std::string,std::vector<gravis::d3Vector>>
			groupAtomsByName(
				const Assembly& assembly,
				const std::set<std::string>& includeExternal);
		
		static std::map<std::string,std::vector<gravis::d3Vector>>
			groupAtomsByElement(
				const Assembly& assembly);
		
		static std::map<std::string,std::map<std::string,std::vector<gravis::d3Vector>>>
			groupAtomsByNameByResidue(
				const Assembly& assembly);
		
		static std::map<std::string,std::vector<gravis::d3Vector>>
			groupAtomsByResidue(
				const Assembly& assembly, std::string element = "");
		
		static std::string getElement(const std::string& atomName);
		static gravis::d3Vector getPosition(const Atom& atom);
		
};

#endif

#include "pdb_helper.h"

using namespace gravis;


std::map<std::string,std::vector<d3Vector>> PdbHelper::splitAtomsByElement(
		const std::string& pdbFile)
{
	std::map<std::string,std::vector<d3Vector>> out;
	
	Assembly assembly;
	assembly.readPDB(pdbFile);
	
	Molecule& molecule = assembly.molecules[0];
	
	const int rc = molecule.residues.size();
	
	for (int r = 0; r < rc; r++)
	{
		Residue& residue = molecule.residues[r];
		
		const int ac = residue.atoms.size();
		
		for (int a = 0; a < ac; a++)
		{
			Atom& atom = residue.atoms[a];
			out[getElement(atom.name)].push_back(getPosition(atom));
		}
	}
	
	const int lac = assembly.looseAtoms.size();
	
	for (int la = 0; la < lac; la++)
	{
		Atom& atom = assembly.looseAtoms[la];
		out[getElement(atom.name)].push_back(getPosition(atom));
	}
	
	return out;
}

std::string PdbHelper::getElement(const std::string& atomName)
{
	if (atomName[0] == ' ') return atomName.substr(1, 1);
	else return atomName.substr(0, 2);
}

d3Vector PdbHelper::getPosition(const Atom &atom)
{
	const Matrix1D<RFLOAT> c = atom.coords;
	
	return d3Vector(c(0), c(1), c(2));
}

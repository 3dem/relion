/*
 * assembly.cpp
 *
 *  Created on: Apr 16, 2013
 *      Author: "Sjors H.W. Scheres"
 */

#include "src/assembly.h"

void Atom::clear()
{
    name = "";
    occupancy = bfactor = 0.;
    coords.clear();
}

Matrix1D<RFLOAT> Atom::getCoordinates()
{
	return coords;
}

void Residue::clear()
{
	number = -1;
	name = "";
	atoms.clear();
}

long int Residue::addAtom(std::string atomname, RFLOAT x, RFLOAT y, RFLOAT z, RFLOAT occ, RFLOAT bfac)
{
	Atom atom(atomname);
	atom.coords = vectorR3(x,y,z);
	atom.occupancy = occ;
	atom.bfactor = bfac;
	long int result = atoms.size();
	atoms.push_back(atom);
	return result;
}

void Molecule::clear()
{
	name = "";
	alt_name = "";
	residues.clear();
}

long int Molecule::insertResidue(Residue &res, int pos)
{
	long int result = residues.size();
	residues.insert(residues.begin() + pos, res);
	return result;
}

long int Molecule::addResidue(Residue &res)
{
	long int result = residues.size();
	residues.push_back(res);
	return result;
}

long int Molecule::addResidue(std::string name, int resnum)
{
	Residue residue(name, resnum);
	long int result = residues.size();
	residues.push_back(residue);
	return result;
}

void Molecule::insertResidues(Molecule add, int residue_start, int residue_end)
{

	int add_nResidues, this_nResidues;
	int ires_start = -1, ires_end = -1;

	add_nResidues = add.numberOfResidues();

	if (residue_start < 0 && residue_end < 0)
	{
		// Add whole chain
		ires_start = 0;
		ires_end = add_nResidues-1;
	}
	else
	{
		// Find beginning and ending ires
		for (int ires = 0; ires < add_nResidues; ires++)
		{
			if (residue_start == add.residues[ires].number)
				ires_start = ires;
			if (residue_end == add.residues[ires].number)
				ires_end = ires;
		}
	}

	if (ires_start < 0 || ires_end < 0)
	{
		std::cerr << " ires_start= " << ires_start << " ires_end= " << ires_end << std::endl;
		std::cerr << " residue_start= " << residue_start << " residue_end= " << residue_end << std::endl;
		REPORT_ERROR("OrigamiBuilder::insertBases ERROR: negative ires_start or ires_end");
	}

	for (int ires = ires_start; ires <= ires_end; ires++)
	{
		int my_res = add.residues[ires].number;
		bool have_inserted=false;
		for (int ii = 0; ii < numberOfResidues(); ii++)
		{
			if (residues[ii].number > my_res)
			{
				insertResidue(add.residues[ires], ii);
				have_inserted = true;
				break;
			}
		}
		if (!have_inserted)
			addResidue(add.residues[ires]);
	}

}

void Assembly::clear()
{
     name = "";
     molecules.clear();
}

long int Assembly::addMolecule(std::string _name, std::string alt_name)
{
	Molecule molecule(_name, alt_name);
	long int result = molecules.size();
	molecules.push_back(molecule);
	return result;
}

long int Assembly::addMolecule(Molecule &toadd)
{

	std::string ori_toadd_name = simplify(toadd.name);
	// Check whether the name of this molecule is unique
	// If not add a suffix to it
	bool is_uniq = false;
	int suffix = 0;
	while (!is_uniq)
	{
		suffix++;
		is_uniq = true;
		for (int imol = 0; imol < molecules.size(); imol++)
		{
			if (molecules[imol].name == toadd.name)
			{
				is_uniq = false;
				break;
			}
		}
		if (!is_uniq)
		{
			toadd.name = ori_toadd_name + integerToString(suffix);
		}

	}
	long int result = molecules.size();
	molecules.push_back(toadd);
	return result;

}

long int Assembly::numberOfMolecules() const
{
	return molecules.size();
}

long int Assembly::numberOfResidues() const
{
	long int sum = 0;
	for (int imol = 0; imol < molecules.size(); imol++)
		sum += molecules[imol].residues.size();
	return sum;
}

long int Assembly::numberOfAtoms() const
{
	long int sum = 0;
	for (int imol = 0; imol < molecules.size(); imol++)
		for (int ires = 0; ires < molecules[imol].residues.size(); ires++)
			sum += molecules[imol].residues[ires].atoms.size();
	return sum;
}

void Assembly::printInformation(std::ostream& out) const
{
	out << " Assembly: " << name << std::endl;
	out << " - Number of molecules : " << numberOfMolecules() << std::endl;
	out << " - Number of residues  : " << numberOfResidues() << std::endl;
	out << " - Number of atoms     : " << numberOfAtoms() << std::endl;

}


void Assembly::readPDB(std::string filename, bool use_segid_instead_of_chainid, bool do_sort)
{

	// Clear existing object
	clear();

	std::ifstream fh(filename.c_str(), std::ios_base::in);

	if (fh.fail())
		REPORT_ERROR( (std::string) "Assembly::read: File " + filename + " does not exists" );

	char line[100];
	bool is_sorted = true;
	int old_resnum = -1;
	long int mol_id = -1;
	long int res_id = -1;
	std::string molname, alt_molname, old_molname="";

	fh.seekg(0);

	// Loop over all lines
	while (fh.getline (line, 600))
	{
		// Only look at lines with an ATOM label
		std::string record(line,0,6);
		if (record == "ATOM  " || record == "HETATM")
		{
			const bool isConnected = record == "ATOM  ";
			
			// ======================== OLD VERSION ========================
			/*
			  char snum[5]={'\0'};
		      char atomname[5]={'\0'};
		      char resname[4]={'\0'};
		      char chainID[1]={'\0'};
		      int resnum=-1;
		      char insertion_residue_code;
		      float x,y,z;
		      float occupancy, bfactor;
		      char segID[5]={'\0'}, element[3]={'\0'}, charge[3]={'\0'};

		      int nr= sscanf(line, "ATOM  %5s %4s %3s %1s%4d%1c   %8f%8f%8f%6f%6f      %4s%2s%2s",
		                          snum, atomname, resname, chainID, &resnum, &insertion_residue_code,
		                          &x, &y, &z, &occupancy, &bfactor, segID, element, charge);
		    */
			// ======================== OLD VERSION ========================

			// ============= May 7, 2015 - Shaoda - Modified according to wwPDB Format v3.3, v3.2 and v2.3 ================
			// sscanf: spaces are ignored! May not get variables at correct subscripts in the string!
			// last 3 objects are incorrect in the old version? No segID in PDB documentation...
			if(strlen(line) < 20)
			{
		    	std::string str(line);
		    	REPORT_ERROR("Assembly::readPDB ERROR: too few entries on ATOM/HETATM line:" + str);
			}
			char snum[6] = "", atomname[5] = "", altLoc[2] = "", resname[4] = "", chainID[2] = "";
			int resnum = -1;
			char insertion_residue_code[2] = "";
			float x, y, z, occupancy, bfactor;
			char segID[5] = "", element[3] = "", charge[3] = "";
			/*
			int nr= sscanf(line, "ATOM  %5s %4s%1s%3s %1s%4d%1s   %8f%8f%8f%6f%6f      %4s%2s%2s",
					snum, atomname, altLoc, resname, chainID, &resnum, insertion_residue_code,
					&x, &y, &z, &occupancy, &bfactor, segID, element, charge);
			*/
			int nr = 0;
			nr += sscanf(line + 6, "%5[^\n]s", snum);
			nr += sscanf(line + 12, "%4[^\n]s", atomname);
			nr += sscanf(line + 16, "%1[^\n]s", altLoc);
			nr += sscanf(line + 17, "%3[^\n]s", resname);
			nr += sscanf(line + 21, "%1[^\n]s", chainID);
			nr += sscanf(line + 22, "%4d", &resnum);
			nr += sscanf(line + 26, "%1[^\n]s", insertion_residue_code);
			nr += sscanf(line + 30, "%8f%8f%8f%6f%6f", &x, &y, &z, &occupancy, &bfactor);
			nr += sscanf(line + 72, "%4[^\n]s", segID);
			nr += sscanf(line + 76, "%2[^\n]s", element);
			nr += sscanf(line + 78, "%2[^\n]s", charge);

			snum[5] = '\0'; atomname[4] = '\0'; altLoc[1] = '\0'; resname[3] = '\0'; chainID[1] = '\0';
			insertion_residue_code[1] = '\0';
			segID[4] = '\0'; element[2] = '\0'; charge[2] = '\0';

			/*
	    	std::cout << "snum = " << snum << ", "
	    			<< "atomname = " << atomname << ", "
					<< "altLoc = " << altLoc << ", "
					<< "resname = " << resname << ", "
					<< "chainID = " << chainID << ", "
					<< "resnum = " << resnum << ", "
					<< "insertion_residue_code = " << insertion_residue_code << ", "
					<< "x = " << x << ", "
					<< "y = " << y << ", "
					<< "z = " << z << ", "
					<< "occupancy = " << occupancy << ", "
					<< "bfactor = " << bfactor << ", "
					<< "segID = " << segID << ", "
					<< "element = " << element << ", "
					<< "charge = " << charge << ", " << std::endl;
			*/
			// ============= May 7, 2015 - Shaoda - Modified according to wwPDB Format v3.3, v3.2 and v2.3 ================

//#define DEBUG
		      if (nr < 5)
		      {
		    	  std::string str(line);
		    	  REPORT_ERROR("Assembly::readPDB ERROR: too few entries on ATOM line:" + str);
		      }

		      if (resnum < 0)
		      {
		    	  REPORT_ERROR("Assembly::readPDB ERROR: negative residue number encountered");
		      }

		      std::string str_chainID(chainID);
		      std::string str_segID(segID);
		      std::string str_atomname(atomname);
		      std::string str_resname(resname);

		      // 1. Get mol_id: to which molecule does this atom belong?
		      // Allow for non-ordered atoms belonging to the same molecule...
		      // To speed up things: first check whether the chainID/segID is the same as the previous line
		      molname = (use_segid_instead_of_chainid) ? str_segID : str_chainID;
		      alt_molname = (use_segid_instead_of_chainid) ? str_chainID : str_segID;
#ifdef DEBUG
		      std::cerr << " molname= " << molname << " alt_molname= " << alt_molname << " str_chainID= " << str_chainID << " chainID= "<< chainID<<std::endl;
#endif
			  if (isConnected)
			  {
				  if (molname != old_molname)
				  {
					  // Check whether a molecule with the same name already exists
					  mol_id = -1;
					  for (long int imol = 0; imol < molecules.size(); imol++)
					  {
						  if (molname == molecules[imol].name)
						  {
							  mol_id = imol;
							  break;
						  }
					  }
	
					  // If not found, add a new molecule
					  if (mol_id < 0)
					  {
	#ifdef DEBUG
						  std::cerr << " snum= " << snum << " atomname= " << atomname << " resname= " << resname ;
						  std::cerr << " chainID= " << chainID << " resnum= " << resnum << " x= " << x ;
						  std::cerr << " y= " << y << " z= " << z << " bfactor= " << bfactor << std::endl;
	#endif
						  mol_id = addMolecule(molname, alt_molname);
						  old_resnum = -1;
					  }
					  else
					  {
						  int lastres = molecules[mol_id].residues.size()-1;
						  old_resnum = (molecules[mol_id]).residues[lastres].number;
					  }
				  }
	
				  // 2. Check whether this is a new residue
				  // All Atoms inside one Residue are assumed to be together!
				  if (resnum != old_resnum)
				  {
					  res_id = molecules[mol_id].addResidue(str_resname, resnum);
					  if (resnum < old_resnum && molname == old_molname && do_sort)
					  {
						  std::cout <<" Warning unsorted residues: " << resnum << " molname= " << molname << " old_resnum= " << old_resnum << " old_molname= " <<old_molname << ""<<std::endl;
						  is_sorted = false;
					  }
				  }
	
				  // 3. Add the actual atom
				  molecules[mol_id].residues[res_id].addAtom(atomname, x, y, z, occupancy, bfactor);
				  
				  // For use in next line
				  old_molname = molname;
				  old_resnum = resnum;
			  }
			  else
			  {
				  Atom atom(atomname);
				  atom.coords = vectorR3(x,y,z);
				  atom.occupancy = occupancy;
				  atom.bfactor = bfactor;
				  looseAtoms.push_back(atom);
			  }
		}
		else
		{
			// At TER statements: reset old_molname to prevent sorting
			old_molname = "";
			old_resnum = 0;
		}

	}

	fh.close();

	// If the Residues were not sorted, do this now:
	if (!is_sorted && do_sort)
	{
		std::cout << " Resorting residues in input PDB file " << filename << " ... " << std::endl;
		sortResidues();
	}

}

// Write PDB format
void Assembly::writePDB(std::string filename)
{
	FILE *file;
	file = fopen(filename.c_str(), "w");
	if (file==NULL)
		REPORT_ERROR("Error writing file: " + filename);

	fprintf(file, "%s\n", "REMARK Created by DsigNA");
	long int atomnum = 0;
	for (int imol = 0; imol < molecules.size(); imol++)
	{
		for (int ires = 0; ires < molecules[imol].residues.size(); ires++)
		{
			for (int iatom = 0; iatom < molecules[imol].residues[ires].atoms.size(); iatom++, atomnum++)
			{
				// If more than 100,000 atoms: just start counting at one again.
				if (atomnum > 99999)
					atomnum -= 99999;
				char chainID = molecules[imol].name[0];

				// ======================== OLD VERSION ========================
 				/*
				fprintf(file, "ATOM  %5d %-4s %3s %1c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f      %4s\n",
						                          atomnum, molecules[imol].residues[ires].atoms[iatom].name.c_str(), molecules[imol].residues[ires].name.c_str(),
						                          chainID, molecules[imol].residues[ires].number,
						                          XX(molecules[imol].residues[ires].atoms[iatom].coords),
						                          YY(molecules[imol].residues[ires].atoms[iatom].coords),
						                          ZZ(molecules[imol].residues[ires].atoms[iatom].coords),
						                          molecules[imol].residues[ires].atoms[iatom].occupancy,
						                          molecules[imol].residues[ires].atoms[iatom].bfactor,
						                          molecules[imol].name.c_str());
				*/
				// ======================== OLD VERSION ========================

				// ============= May 7, 2015 - Shaoda - Modified according to wwPDB Format v3.3, v3.2 and v2.3 ================
				// last 3 objects are incorrect in the old version? No segID in PDB documentation...
				char atomname[5] = "", element[2] = "";
				strcpy(atomname, molecules[imol].residues[ires].atoms[iatom].name.c_str());
				element[0] = ' '; element[1] = '\0';
				for(int ii = 0; ii < 4; ii++)
				{
					if(atomname[ii] != ' ')
					{
						element[0] = atomname[ii];
						break;
					}
				}
				fprintf(file, "ATOM  %5ld %-4s %3s %1c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f          %2s  \n",
						                          atomnum, molecules[imol].residues[ires].atoms[iatom].name.c_str(), molecules[imol].residues[ires].name.c_str(),
						                          chainID, molecules[imol].residues[ires].number,
						                          XX(molecules[imol].residues[ires].atoms[iatom].coords),
						                          YY(molecules[imol].residues[ires].atoms[iatom].coords),
						                          ZZ(molecules[imol].residues[ires].atoms[iatom].coords),
						                          molecules[imol].residues[ires].atoms[iatom].occupancy,
						                          molecules[imol].residues[ires].atoms[iatom].bfactor,
												  element);
				// ============= May 7, 2015 - Shaoda - Modified according to wwPDB Format v3.3, v3.2 and v2.3 ================

			}
		}
		if (imol + 1 < molecules.size())
			fprintf(file, "%s\n", "TER");
	}

	fprintf(file, "%s\n", "END");
	fclose(file);
}

void Assembly::join(Assembly &tojoin)
{

	for (int imol = 0; imol < tojoin.molecules.size(); imol++)
		addMolecule(tojoin.molecules[imol]);

}

void Assembly::sortResidues()
{

	// Loop over all molecules
	for (int imol = 0; imol < molecules.size(); imol++)
	{

		// A. Sort all Residues
		std::vector<std::pair<int, int> > vp;
		for (int ires = 0; ires < molecules[imol].residues.size(); ires++)
		{
			vp.push_back(std::make_pair(molecules[imol].residues[ires].number, ires));
		}
		std::sort(vp.begin(), vp.end());
		std::vector<Residue> new_residues;
		for (int ires = 0; ires < molecules[imol].residues.size(); ires++)
		{
			new_residues.push_back(molecules[imol].residues[vp[ires].second]);
		}
		molecules[imol].residues = new_residues;
	}

}

void Assembly::applyTransformation(Matrix2D<RFLOAT> &mat, Matrix1D<RFLOAT> &shift)
{
	for (int imol = 0; imol < molecules.size(); imol++)
	{
		for (int ires = 0; ires < molecules[imol].residues.size(); ires++)
		{
			for (int iatom = 0; iatom < molecules[imol].residues[ires].atoms.size(); iatom++)
			{
				(molecules[imol].residues[ires].atoms[iatom]).coords = mat * (molecules[imol].residues[ires].atoms[iatom]).coords;
				(molecules[imol].residues[ires].atoms[iatom]).coords += shift;
			}
		}
	}

}

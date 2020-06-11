/*
 * assembly.h
 *
 *  Created on: Apr 16, 2013
 *      Author: "Sjors H.W. Scheres"
 */

#ifndef ASSEMBLY_H_
#define ASSEMBLY_H_
#include <cstdlib>
#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include "src/args.h"
#include "src/matrix2d.h"

/*
 * Hierarchical model for a macromolecular assembly, e.g. a DNA origami object
 *
 *  Assembly
 *   ->	Molecule
 *   	->	Residue
 *   		->	Atom  (either a true one or a coarse-grain pseudo-atom)
 *
 *
 */

class Atom
{
public:
        // Name of this Atom
        std::string name;

        // Coordinates
        Matrix1D<RFLOAT> coords;

        // Occupancy
        RFLOAT occupancy;

        // B-factor
        RFLOAT bfactor;

        // Empty constructor
        Atom()
        {
        	clear();
        }

        // Named constructor
        Atom(std::string in_name)
        {
        	clear();
        	name = in_name;
        }

        // Destructor needed for work with vectors
        ~Atom()
        {
        	clear();
        }

        // Initialize
        void clear();

        // Get the 3D corrdinates as a POint3D
        Matrix1D<RFLOAT> getCoordinates();
};


class Residue
{
public:
        // Name of this Residue
        std::string name;

        // Number of this Residue
        int number;

        // All the Atoms in this Residue
        std::vector<Atom> atoms;

        // Empty Constructor
        Residue()
        {
        	clear();
        }

        // Constructor
        Residue(std::string in_name, int in_number)
        {
        	clear();
        	name = in_name;
        	number = in_number;
        }

        // Destructor needed for work with vectors
        ~Residue()
        {
        	clear();
        }

        // Initialize
        void clear();

        // Add an Atom to this Residue;
        long int addAtom(std::string atomname, RFLOAT x, RFLOAT y, RFLOAT z, RFLOAT occ = 1.0, RFLOAT bfac = 0.0);

        int numberOfAtoms()
        {
        	return atoms.size();
        }

};

class Molecule
{
public:
        // Name of this Molecule
        std::string name;

        // Alternative name of this Molecule (either chainID or segID)
        std::string alt_name;

        // All the Residues in this Molecule
        std::vector<Residue> residues;

        // Empty Constructor
        Molecule()
        {
        	clear();
        }

        // Constructor
        Molecule(std::string in_name, std::string in_alt_name="")
        {
        	clear();
        	name = in_name;
        	alt_name = in_alt_name;
        }

        // Destructor needed for work with vectors
        ~Molecule()
        {
        	clear();
        }

        // Initialize
        void clear();

        // Number of residues in the molecule
        long int numberOfResidues()
        {
        	return residues.size();
        }

        // Insert a Residue at the specified position in this Molecule
        long int insertResidue(Residue &res, int pos);

        // Add a Residue to this Molecule
        long int addResidue(Residue &res);

        // Add a Residue to this Molecule
        long int addResidue(std::string name, int resnum);

        // Insert a stretch of residues from another Molecule based on consecutive residue numbering
        // If start and end residues are negative: just add the entire molecule
        void insertResidues(Molecule add, int residue_start = -1, int residue_end = -1);

};

class Assembly
{
public:
        // Name of this Assembly
        std::string name;

        // All the Molecules in this Assembly
        std::vector<Molecule> molecules;
		
		// Additional atoms (denoted by HETATM in the PDB file)
		std::vector<Atom> looseAtoms;

        // Empty Constructor
        Assembly()
        {
        	clear();
        }

        // Named Constructor
        Assembly(std::string in_name)
        {
        	clear();
        	name = in_name;
        }

    	// Copy constructor
        Assembly(const Assembly& op)
    	{
    		clear();
    		*this = op;
    	}

    	// Destructor needed for work with vectors
        ~Assembly()
        {
        	clear();
        }

        // Initialize
        void clear();

        // Add a Molecule to this Assembly
        long int addMolecule(std::string name, std::string alt_name);

        // Add a Molecule to this Assembly
        long int addMolecule(Molecule &toadd);

        // return number of Molecules in the Assembly
        long int numberOfMolecules() const;

        // Total number of Atoms
        long int numberOfAtoms() const;

        // Total number of Residues
        long int numberOfResidues() const;

        // Print some information about the assembly
        void printInformation(std::ostream& out = std::cout) const;

        // Read PDB format
        void readPDB(std::string filename, bool use_segid_instead_of_chainid = false, bool do_sort = true);

        // Write the Assembly to a PDB file
        void writePDB(std::string filename);

        // Combine this Assembly with another one
        // If there are identical Molecule.name instances, add a number-suffix to the new Assembly's Molecule.name (in the segID)
        void join(Assembly &tojoin);

        // Make sure that all Residues within each Molecule are in order w.r.t. their residue number
        void sortResidues();

        // Break Molecules into separate ones if a break larger than maximum_residue_break occurs in the residue numbering
        // TODO
        void checkBreaksInResidueNumbering(int maximum_residue_break = 500);

        // Apply a transformation (first rotation, then shift)
        void applyTransformation(Matrix2D<RFLOAT> &mat, Matrix1D<RFLOAT> &shift);

};


#endif /* ASSEMBLY_H_ */

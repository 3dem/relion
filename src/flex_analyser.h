/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#ifndef SRC_FLEX_ANALYSER_H_
#define SRC_FLEX_ANALYSER_H_

#include "src/exp_model.h"
#include "src/ml_model.h"
#include "src/ctf.h"
#include "src/time.h"
#include "src/parallel.h"
#include "src/mpi.h"

class FlexAnalyser
{
public:
	// I/O Parser
	IOParser parser;

	//verbosity
	int verb;

	// Output rootname
	FileName fn_out;

	// The model and the data from the refinement to be analysed
	FileName fn_model, fn_data;
	MlModel model;
	Experiment data;

	// The body STAR file
	FileName fn_bodies;

	// Write out 3D models
	bool do_3dmodels;

	// Box size of the output 3D models
	int size_3dmodels;

	// Rescale factor for output 3D models
	RFLOAT rescale_3dmodels;

	// Perform a PCA on the multibody orientations
	bool do_PCA_orient;

	// Normalisation of the rotations and translations for PCA normalisation
	std::vector<float> norm_pca;

	// Generate maps for movies along principal components
	bool do_generate_maps;

	// How many components to make movies from?
	int nr_components;

	// How much variance to explain with the movies?
	double explain_variance;

	// How many maps to use for the movie of each principal component?
	int nr_maps_per_component;

	// How many bins in a histogram
	int nr_bins;

	// Select particles based on this eigenvalue
	int select_eigenvalue;

	// Select particles based on this eigenvalue minimim
	float select_eigenvalue_min;

	// Select particles based on this eigenvalue minimim
	float select_eigenvalue_max;

	// Write out text file with eigenvalues for all particles
	bool do_write_all_pca_projections;

	// center of mass of the above
	Matrix1D<RFLOAT> com_mask;

	// Pre-calculated rotation matrix for (0,90,0) rotation, and its transpose
	Matrix2D<RFLOAT> A_rot90, A_rot90T;

	MetaDataTable DFo;

	void read(int argc, char **argv);

	void initialise();

	void run(int rank = 0, int size = 1);

	void setupSubtractionMasksAndProjectors();

	void setup3DModels();

	void loopThroughParticles(int rank = 0, int size = 1);

	void subtractOneParticle(long int part_id, long int imgno, int rank = 0, int size = 1);
	void make3DModelOneParticle(long int part_id, long int imgno, std::vector<double> &datarow, int rank = 0, int size = 1);

	// Output logfile.pdf with histograms of all eigenvalues
	void makePCAhistograms(std::vector< std::vector<double> > &projected_input,
	                       std::vector<double> &eigenvalues, std::vector<double> &means);

	// Generate maps to make movies of the variance along the most significant eigenvectors
	void make3DModelsAlongPrincipalComponents(std::vector< std::vector<double> > &projected_input,
	                                          std::vector< std::vector<double> > &eigenvectors, std::vector<double> &means);

	// Dump all projections to a text file
	void writeAllPCAProjections(std::vector< std::vector<double> > &projected_input);

	// Output a particle.star file with a selection based on eigenvalues
	void outputSelectedParticles(std::vector< std::vector<double> > &projected_input);

};

void principalComponentsAnalysis(const std::vector< std::vector<double> > &input,
                                 std::vector< std::vector<double> > &eigenvectors,
                                 std::vector<double> &eigenvalues, std::vector<double> &means,
                                 std::vector< std::vector<double> > &projected_input);

#endif /* SRC_FLEX_ANALYSER_H_ */

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

#include "flex_analyser.h"

void FlexAnalyser::read(int argc, char **argv)
{
	parser.setCommandLine(argc, argv);

	int gen_section = parser.addSection("General options");
	fn_data = parser.getOption("--data", "The _data.star file with the orientations to be analysed", "");
	fn_model = parser.getOption("--model", " The corresponding _model.star file with the refined model", "");
	fn_bodies = parser.getOption("--bodies", "The corresponding star file with the definition of the bodies", "");
	fn_out = parser.getOption("--o", "Output rootname", "analyse");

	int model_section = parser.addSection("3D model options");
	do_3dmodels = parser.checkOption("--3dmodels", "Generate a 3D model for each experimental particles");
	size_3dmodels = textToInteger(parser.getOption("--size_3dmodels", "Output size of the 3D models (default is same as input particles)", "-1"));

	int pca_section = parser.addSection("PCA options");
	do_PCA_orient = parser.checkOption("--PCA_orient", "Perform a principal components analysis on the multibody orientations");
	do_generate_maps = parser.checkOption("--do_maps", "Generate maps along the principal components");
	nr_components = textToInteger(parser.getOption("--k", "Number of principal components to generate maps for", "-1"));
	explain_variance  = textToFloat(parser.getOption("--v", "Or use as many principal components to explain this fraction of variance (<0,1])", "0.75"));
	nr_maps_per_component = textToInteger(parser.getOption("--maps_per_movie", "Number of maps to use for the movie of each principal component", "10"));
	nr_bins = textToInteger(parser.getOption("--bins", "Number of bins in histograms of the eigenvalues for each principal component", "100"));
	select_eigenvalue = textToInteger(parser.getOption("--select_eigenvalue", "Output a selection particle.star file based on eigenvalues along this eigenvector", "-1"));
	select_eigenvalue_min = textToFloat(parser.getOption("--select_eigenvalue_min", "Minimum for eigenvalue to include particles in selection output star file", "-99999."));
	select_eigenvalue_max = textToFloat(parser.getOption("--select_eigenvalue_max", "Maximum for eigenvalue to include particles in selection output star file", "99999."));
	do_write_all_pca_projections = parser.checkOption("--write_pca_projections", "Write out a text file with all PCA projections for all particles");

	// Initialise verb for non-parallel execution
	verb = textToInteger(parser.getOption("--verb", "Verbosity", "1"));

	// Check for errors in the command-line option
	if (parser.checkForErrors())
		REPORT_ERROR("Errors encountered on the command line (see above), exiting...");
}

void FlexAnalyser::initialise()
{
	rescale_3dmodels = 1.0;

	if (verb > 0)
		std::cout << " Reading in data.star file ..." << std::endl;

	if (fn_data == "")
		REPORT_ERROR("ERROR: please provide the --data argument!");
	else
		data.read(fn_data);

	if (verb > 0)
		std::cout << " Reading in model.star file ..." << std::endl;

	if (fn_model == "")
		REPORT_ERROR("ERROR: please provide the --model argument!");
	else
		model.read(fn_model, data.obsModel.numberOfOpticsGroups());

	if (fn_bodies != "")
	{
		if (verb > 0)
			std::cout << " Initialising bodies ..." << std::endl;
		model.initialiseBodies(fn_bodies, fn_out);
	}
	else {
		REPORT_ERROR("ERRPR: please specify the --bodies argument!");
	}

	if (model.nr_bodies != data.nr_bodies)
		REPORT_ERROR("ERROR: Unequal number of bodies in bodies.star and data.star files!");

	if (do_3dmodels && model.nr_bodies == 1)
		REPORT_ERROR("ERROR: --3dmodels option is only valid for multibody refinements.");

	// This creates a rotation matrix for (rot,tilt,psi) = (0,90,0)
	// It will be used to make all Abody orientation matrices relative to (0,90,0) instead of the more logical (0,0,0)
	// This is useful, as psi-priors are ill-defined around tilt=0, as rot becomes the same as -psi!!
	rotation3DMatrix(-90., 'Y', A_rot90, false);
	A_rot90T = A_rot90.transpose();

	if (do_PCA_orient)
	{
		if (model.nr_bodies * 6 > data.numberOfParticles())
			REPORT_ERROR("ERROR: there are not enough particles to perform PCA!");

		if (do_generate_maps)
		{
			if (explain_variance > 1.)
				REPORT_ERROR("ERROR: --v should be expressed as a fraction, i.e. between 0 and 1.");
			if (explain_variance < 0. && nr_components < 0)
				REPORT_ERROR("ERROR: --v or --k should be larger than zero.");
		}

		// Calculate effect of 1 degree rotations and 1 pixel translations on the bodies, in order to normalise vectors for PCA
		norm_pca.clear();

		FileName fn_weights = fn_out + "_pca_weights.dat";
		std::ofstream f_weights(fn_weights);
		std::cout << " Normalisation weights for PCA columns are written to " << fn_weights << std::endl;
		f_weights << "body rot tilt psi offset" << std::endl;
		f_weights << std::scientific;
		for (int ibody = 0; ibody < model.nr_bodies; ibody++)
		{
			MultidimArray<RFLOAT> Mbody, Irefp;
			Irefp = model.Iref[ibody] * model.masks_bodies[ibody];
			// Place each body with its center-of-mass in the center of the box
			selfTranslate(Irefp, -model.com_bodies[ibody], DONT_WRAP);

			f_weights << ibody + 1;
			Matrix2D<RFLOAT> Aresi,  Abody;
			// rot
			Euler_angles2matrix(1., 90., 0., Aresi);
			Abody = (model.orient_bodies[ibody]).transpose() * A_rot90 * Aresi * model.orient_bodies[ibody];
			Abody.resize(4,4);
			MAT_ELEM(Abody, 3, 3) = 1.;
			applyGeometry(Irefp, Mbody, Abody, IS_NOT_INV, DONT_WRAP);
			Mbody -= Irefp;
			norm_pca.push_back(sqrt(Mbody.sum2()));
			f_weights << " " << sqrt(Mbody.sum2());
			// tilt
			Euler_angles2matrix(0., 91., 0., Aresi);
			Abody = (model.orient_bodies[ibody]).transpose() * A_rot90 * Aresi * model.orient_bodies[ibody];
			Abody.resize(4,4);
			MAT_ELEM(Abody, 3, 3) = 1.;
			applyGeometry(Irefp, Mbody, Abody, IS_NOT_INV, DONT_WRAP);
			Mbody -= Irefp;
			norm_pca.push_back(sqrt(Mbody.sum2()));
			f_weights << " " << sqrt(Mbody.sum2());
			// psi
			Euler_angles2matrix(0., 90., 1., Aresi);
			Abody = (model.orient_bodies[ibody]).transpose() * A_rot90 * Aresi * model.orient_bodies[ibody];
			Abody.resize(4,4);
			MAT_ELEM(Abody, 3, 3) = 1.;
			applyGeometry(Irefp, Mbody, Abody, IS_NOT_INV, DONT_WRAP);
			Mbody -= Irefp;
			norm_pca.push_back(sqrt(Mbody.sum2()));
			f_weights << " " << sqrt(Mbody.sum2());
			// translation x & y (considered the same)
			Euler_angles2matrix(0., 90., 0., Aresi);
			Abody = (model.orient_bodies[ibody]).transpose() * A_rot90 * Aresi * model.orient_bodies[ibody];
			Abody.resize(4,4);
			MAT_ELEM(Abody, 0, 3) = 1.;
			MAT_ELEM(Abody, 3, 3) = 1.;
			applyGeometry(Irefp, Mbody, Abody, IS_NOT_INV, DONT_WRAP);
			Mbody -= Irefp;
			norm_pca.push_back(sqrt(Mbody.sum2()));
			f_weights << " " << sqrt(Mbody.sum2()) << std::endl;
		}
		f_weights.close();
	}
}

void FlexAnalyser::run(int rank, int size)
{
	if (size > 1 && do_PCA_orient)
		REPORT_ERROR("PCA analysis (--PCA_orient) must be performed in the non-MPI version.");

	if (do_3dmodels) setup3DModels();

	// Loop through all particles
	loopThroughParticles(rank, size);

	if (size > 1) {
		MPI_Barrier(MPI_COMM_WORLD);
	}
}

void FlexAnalyser::setup3DModels()
{
	for (int ibody = 0; ibody < model.nr_bodies; ibody++)
	{
		// Premultiply the map with the mask (otherwise need to do this again for every particle
		model.Iref[ibody] *= model.masks_bodies[ibody];
		// Place each body with its center-of-mass in the center of the box, as that's where the rotations are around
		selfTranslate(model.Iref[ibody], -model.com_bodies[ibody], DONT_WRAP);
		// And do the same for the masks
		selfTranslate(model.masks_bodies[ibody], -model.com_bodies[ibody], DONT_WRAP);

		if (size_3dmodels < XSIZE(model.Iref[ibody]))
		{
			rescale_3dmodels = (RFLOAT)(size_3dmodels)/(RFLOAT)(XSIZE(model.Iref[ibody]));
			std::cerr << " rescale_3dmodels= " << rescale_3dmodels << std::endl;
			selfScaleToSize(model.Iref[ibody], size_3dmodels, size_3dmodels, size_3dmodels);
			selfScaleToSize(model.masks_bodies[ibody], size_3dmodels, size_3dmodels, size_3dmodels);
			model.Iref[ibody].setXmippOrigin();
			model.masks_bodies[ibody].setXmippOrigin();
		}
	}
}

void FlexAnalyser::loopThroughParticles(int rank, int size)
{
	long int total_nr_particles = data.numberOfParticles();
	// Allow parallelisation
	long int my_first_particle = 0, my_last_particle = total_nr_particles-1;
	if (size > 1)
		divide_equally(total_nr_particles, size, rank, my_first_particle, my_last_particle);
	long int todo_particles = my_last_particle-my_first_particle+1;
	long int update_interval = XMIPP_MAX(1, todo_particles / 60);
	if (verb > 0)
	{
		std::cout << " Processing all particles ... " << std::endl;
		init_progress_bar(todo_particles);
	}

	DFo.clear();
	DFo.setIsList(false);

	std::vector< std::vector<double> > inputdata;
	long int imgno = 0;
	for (long int part_id = my_first_particle; part_id <= my_last_particle; part_id++)
	{

		std::vector<double> datarow;
		if (do_3dmodels || do_PCA_orient)
		{
			make3DModelOneParticle(part_id, imgno, datarow, rank, size);
			if (do_PCA_orient)
				inputdata.push_back(datarow);
		}

		if (imgno%update_interval==0 && verb > 0)
			progress_bar(imgno);
		imgno++;
	}
	if (verb > 0)
		progress_bar(todo_particles);

	if (do_3dmodels)
	{
		FileName fn_star;
		if (size > 1) {
			fn_star.compose(fn_out + "_", rank + 1, "");
			fn_star = fn_star + "_3dmodels.star";
		} else {
			fn_star = fn_out + "_3dmodels.star";
		}
		DFo.write(fn_star);
	}

	if (do_PCA_orient)
	{
		std::vector< std::vector<double> > eigenvectors, projected_data;
		std::vector<double> eigenvalues, means;
		// Do the PCA and make histograms
		principalComponentsAnalysis(inputdata, eigenvectors, eigenvalues, means, projected_data);

		FileName fn_evec = fn_out + "_eigenvectors.dat";
		std::ofstream f_evec(fn_evec);
		std::cout << " Eigenvectors (rotations only):" << std::endl;
		for (int j = 0; j < eigenvectors[0].size(); j++)
		{
			std::string stro = "";
			if (j % 6 == 0)
				stro = "rot";
			else if (j % 6 == 1)
				stro = "tilt";
			else if (j % 6 == 2)
				stro = "psi";
			else if (j % 6 == 3)
				stro = "x";
			else if (j % 6 ==  4)
				stro = "y";
			else if (j % 6 == 5)
				stro = "z";
			if (stro != "")
			{
				stro +=  "-body-" + integerToString(1 + (j / 6));
				f_evec << stro << " ";
				if (j % 6 < 3) {
					std::cout << std::setw(12) << std::right << std::fixed;
					std::cout << stro;
				}
			}
		}
		std::cout << std::endl;
		f_evec << std::endl;
		std::cout << " Full eigenvectors including translations are written to " << fn_evec << std::endl;

		f_evec << std::scientific;
		for (int k = 0; k < eigenvectors.size(); k++)
		{
			for (int j =0; j < eigenvectors[0].size(); j++)
			{
				if (j > 0) f_evec << " ";
				f_evec << eigenvectors[k][j];
			}
			f_evec << std::endl;

			if (k % 6 < 3)
			{
				for (int j =0; j < eigenvectors[0].size(); j++)
				{
					if (j % 6 < 3)
					{
						std::cout << std::setw(12) << std::fixed;
						std::cout << eigenvectors[k][j];
					}
				}
				std::cout << std::endl;
			}
		}

		f_evec.close();

		makePCAhistograms(projected_data, eigenvalues, means);

		// Make movies for the most significant eigenvectors
		if (do_generate_maps)
			make3DModelsAlongPrincipalComponents(projected_data, eigenvectors, means);

		if (do_write_all_pca_projections)
		{
			writeAllPCAProjections(projected_data);
		}

		// Output a particle selection, if requested
		if (select_eigenvalue > 0)
		{
			outputSelectedParticles(projected_data);
		}
	}
}

void FlexAnalyser::make3DModelOneParticle(long int part_id, long int imgno, std::vector<double> &datarow, int rank, int size)
{
	// Get the consensus class, orientational parameters and norm (if present)
	Matrix2D<RFLOAT> Aori;
	RFLOAT rot, tilt, psi, xoff, yoff, zoff;
	data.MDimg.getValue(EMDL_ORIENT_ROT, rot, part_id);
	data.MDimg.getValue(EMDL_ORIENT_TILT, tilt, part_id);
	data.MDimg.getValue(EMDL_ORIENT_PSI, psi, part_id);
	Euler_angles2matrix(rot, tilt, psi, Aori, false);

	RFLOAT my_pixel_size = data.getImagePixelSize(part_id, 0);

	Image<RFLOAT> img;
	MultidimArray<RFLOAT> sumw;
	if (do_3dmodels)
	{
		img().initZeros(model.Iref[0]);
		sumw.initZeros(model.Iref[0]);
	}

	datarow.clear();
 	for (int ibody = 0; ibody < model.nr_bodies; ibody++)
	{
		MultidimArray<RFLOAT> Mbody, Mmask;
		Matrix1D<RFLOAT> body_offset(3), body_offset_3d(3);
		RFLOAT body_rot, body_tilt, body_psi;
		data.MDbodies[ibody].getValue(EMDL_ORIENT_ROT, body_rot, part_id);
		data.MDbodies[ibody].getValue(EMDL_ORIENT_TILT, body_tilt, part_id);
		data.MDbodies[ibody].getValue(EMDL_ORIENT_PSI,  body_psi, part_id);
		data.MDbodies[ibody].getValue(EMDL_ORIENT_ORIGIN_X_ANGSTROM, XX(body_offset), part_id);
		data.MDbodies[ibody].getValue(EMDL_ORIENT_ORIGIN_Y_ANGSTROM, YY(body_offset), part_id);
		if (model.data_dim == 3)
			data.MDbodies[ibody].getValue(EMDL_ORIENT_ORIGIN_Z_ANGSTROM, ZZ(body_offset), part_id);

		// As of v3.1, offsets are in Angstrom: convert back to pixels!
		body_offset /= my_pixel_size;

		// Keep rescaling into account!
		body_offset *= rescale_3dmodels;

		Matrix2D<RFLOAT> Aresi,  Abody, Anew;
		// Aresi is the residual orientation for this ibody
		Euler_angles2matrix(body_rot, body_tilt, body_psi, Aresi);
		// Only apply the residual orientation now!!!
		Abody = (model.orient_bodies[ibody]).transpose() * A_rot90 * Aresi * model.orient_bodies[ibody];

		// Now we have to get back from the 2D refined body_offset to some 3D translation of the body (with one direction non-defined)
		// We will need the original projection direction, Aori for that!!
		// Because one direction is ill-defined, this may not be such a good idea?
		// But anyway, this should bring it closer to truth than not doing anything at all...
		Anew = Aori * Abody;
		body_offset_3d = Anew.inv() * (-body_offset);

		if (do_PCA_orient)
		{
			datarow.push_back(norm_pca[ibody*4+0] * body_rot);
			datarow.push_back(norm_pca[ibody*4+1] * body_tilt);
			datarow.push_back(norm_pca[ibody*4+2] * body_psi);
			datarow.push_back(norm_pca[ibody*4+3] * XX(body_offset_3d));
			datarow.push_back(norm_pca[ibody*4+3] * YY(body_offset_3d));
			datarow.push_back(norm_pca[ibody*4+3] * ZZ(body_offset_3d));
		}

		if (do_3dmodels)
		{
			// Also put back at the centre-of-mass of this body
			body_offset_3d += rescale_3dmodels * model.com_bodies[ibody];
			Abody.resize(4,4);

			MAT_ELEM(Abody, 0, 3) = XX(body_offset_3d);
			MAT_ELEM(Abody, 1, 3) = YY(body_offset_3d);
			MAT_ELEM(Abody, 2, 3) = ZZ(body_offset_3d);
			MAT_ELEM(Abody, 3, 3) = 1.;

			Mbody.resize(model.Iref[ibody]);
			Mmask.resize(model.masks_bodies[ibody]);
			applyGeometry(model.Iref[ibody], Mbody, Abody, IS_NOT_INV, DONT_WRAP);
			applyGeometry(model.masks_bodies[ibody], Mmask, Abody, IS_NOT_INV, DONT_WRAP);

			img() += Mbody;
			sumw += Mmask;
		}

	} // end for ibody

	if (do_3dmodels)
	{
		// Divide the img by sumw to deal with overlapping bodies: just take average
		FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(img())
		{
			if (DIRECT_MULTIDIM_ELEM(sumw, n) > 1.)
				DIRECT_MULTIDIM_ELEM(img(), n) /= DIRECT_MULTIDIM_ELEM(sumw, n);
		}
		// Write the image to disk
		FileName fn_img;
		fn_img.compose(fn_out+"_part", imgno+1,"mrc");
		img.setSamplingRateInHeader(model.pixel_size);
		img.write(fn_img);

		DFo.addObject();
		DFo.setValue(EMDL_MLMODEL_REF_IMAGE, fn_img);
		data.MDimg.getValue(EMDL_IMAGE_NAME, fn_img, part_id);
		DFo.setValue(EMDL_IMAGE_NAME, fn_img);
	}
}

void FlexAnalyser::makePCAhistograms(std::vector< std::vector<double> > &projected_input,
                                     std::vector<double> &eigenvalues, std::vector<double> &means)
{
	std::vector<FileName> all_fn_eps;
	FileName fn_eps = fn_out + "_eigenvalues.eps";
	all_fn_eps.push_back(fn_eps);
	CPlot2D *plot2D=new CPlot2D(fn_eps);
	CDataSet dataSet;
	dataSet.SetDrawMarker(false);
	dataSet.SetDatasetColor(1.0,0.0,0.0);

	// Percentage of variance
	double sum=0.;
	for (int i = 0; i < eigenvalues.size(); i++)
		sum += eigenvalues[i];
	for (int i = 0; i < eigenvalues.size(); i++)
	{
		std::cout << "  + Component " << i+1 << " explains " << eigenvalues[i]*100./sum << "% of variance." << std::endl;
		CDataPoint point1((double)i+0.5, (double)0.);
		CDataPoint point2((double)i+0.5, (double)eigenvalues[i]*100./sum);
		CDataPoint point3((double)i+1.5, (double)eigenvalues[i]*100./sum);
		CDataPoint point4((double)i+1.5, (double)0.);
		dataSet.AddDataPoint(point1);
		dataSet.AddDataPoint(point2);
		dataSet.AddDataPoint(point3);
		dataSet.AddDataPoint(point4);
	}
	plot2D->AddDataSet(dataSet);
	plot2D->SetXAxisTitle("Eigenvalue");
	plot2D->SetYAxisTitle("Variance explained [%]");
	plot2D->OutputPostScriptPlot(fn_eps);
	delete plot2D;

	// Determine how much variance the requested number of components explains
	if (nr_components < 0)
	{
		double cum = 0.;
		for (int i = 0; i < eigenvalues.size(); i++)
		{
			cum += eigenvalues[i]/sum;
			if (cum >= explain_variance)
			{
				nr_components = i + 1;
				break;
			}
		}
	}

	explain_variance = 0.;
	for (int i = 0; i < nr_components; i++)
		explain_variance += eigenvalues[i]*100./sum;

	std::cout << " The first " << nr_components << " eigenvectors explain " << explain_variance << " % of the variance in the data." << std::endl;

	// Output histograms of all eigenvalues
	for (int k = 0; k < eigenvalues.size(); k++)
	{
		// Sort vector of all projected values for this component: divide in nr_maps_per_component bins and take average value
		std::vector<double> project;
		for (long int ipart = 0; ipart < projected_input.size(); ipart++)
			project.push_back(projected_input[ipart][k]);

		// Sort the vector to calculate average of nr_maps_per_component equi-populated bins
		std::sort (project.begin(), project.end());

		// Write the movement plot as well
		FileName fn_eps = fn_out + "_component" + integerToString(k+1, 3) + "_histogram.eps";
		all_fn_eps.push_back(fn_eps);
		CPlot2D *plot2D=new CPlot2D(fn_eps);
		CDataSet dataSet;
		dataSet.SetDrawMarker(false);
		dataSet.SetDatasetColor(1.0,0.0,0.0);

		double minhis = project[0];
		double maxhis = project[project.size()-1];
		double widthhis = (maxhis - minhis) / nr_bins;
		double stophis = minhis + widthhis;
		long int n = 0;
		for (long int ipart = 0; ipart < project.size(); ipart++)
		{
			if (project[ipart] >= stophis)
			{
				CDataPoint point1(stophis-widthhis, (double)0.);
				CDataPoint point2(stophis-widthhis, (double)n);
				CDataPoint point3(stophis, (double)n);
				CDataPoint point4(stophis, (double)0.);
				dataSet.AddDataPoint(point1);
				dataSet.AddDataPoint(point2);
				dataSet.AddDataPoint(point3);
				dataSet.AddDataPoint(point4);
				n = 0;
				stophis += widthhis;
			}
			n++;
		}
		plot2D->AddDataSet(dataSet);
		plot2D->SetXAxisTitle("Eigenvalue");
		plot2D->SetYAxisTitle("Nr particles");
		plot2D->OutputPostScriptPlot(fn_eps);
		delete plot2D;
	}

	joinMultipleEPSIntoSinglePDF(fn_out + "_logfile.pdf", all_fn_eps);
}

void FlexAnalyser::make3DModelsAlongPrincipalComponents(std::vector< std::vector<double> > &projected_input,
                                                        std::vector< std::vector<double> > &eigenvectors, std::vector<double> &means)
{
	// Loop over the principal components
	for (int k = 0; k < nr_components; k++)
	{

		// Sort vector of all projected values for this component: divide in nr_maps_per_component bins and take average value
		std::vector<double> project;
		for (long int ipart = 0; ipart < projected_input.size(); ipart++)
			project.push_back(projected_input[ipart][k]);

		// Sort the vector to calculate average of "nr_maps_per_component" equi-populated bins
		std::sort (project.begin(), project.end());

		long int binwidth = ROUND((double)project.size() / (double)nr_maps_per_component);

		std::cout << " Calculating 3D models for principal component " << k+1 << " ... " << std::endl;

		for (int ibin = 0; ibin < nr_maps_per_component; ibin++)
		{
			long int istart = ibin * binwidth;
			long int istop = (ibin+1) * binwidth - 1;
			if (ibin == nr_maps_per_component - 1)
				istop = project.size() - 1;

			double avg = 0., nn = 0.;
			for (long int ipart = istart; ipart <= istop; ipart++)
			{
				avg += project[ipart];
				nn += 1.;
			}
			if (nn > 0.)
				avg /= nn;

			// Now we have the average value for the PCA values for this bin: make the 3D model...
			std::vector<double> orients;
			for (int j = 0; j < means.size(); j++)
			{
				orients.push_back(avg * eigenvectors[k][j] + means[j]);
				//std::cerr << "j= "<<j<< " orients[j]= " << orients[j]<< " === "<<avg<< "  * " <<eigenvectors[k][j] << "  + " << means[j] << std::endl;
			}

			Image<RFLOAT> img;
			MultidimArray<RFLOAT> sumw;
			img().initZeros(model.Iref[0]);
			sumw.initZeros(model.Iref[0]);
			for (int ibody = 0; ibody < model.nr_bodies; ibody++)
			{

				MultidimArray<RFLOAT> Mbody, Mmask;
				Matrix1D<RFLOAT> body_offset_3d(3);
				RFLOAT body_rot, body_tilt, body_psi;
				body_rot            = orients[ibody * 6 + 0] / norm_pca[ibody*4+0];
				body_tilt           = orients[ibody * 6 + 1] / norm_pca[ibody*4+1];
				body_psi            = orients[ibody * 6 + 2] / norm_pca[ibody*4+2];
				XX(body_offset_3d)  = orients[ibody * 6 + 3] / norm_pca[ibody*4+3];
				YY(body_offset_3d)  = orients[ibody * 6 + 4] / norm_pca[ibody*4+3];
				ZZ(body_offset_3d)  = orients[ibody * 6 + 5] / norm_pca[ibody*4+3];
				//std::cerr << " norm_pca[ibody*4+0]= " << norm_pca[ibody*4+0] << " norm_pca[ibody*4+1]= " << norm_pca[ibody*4+1] << " norm_pca[ibody*4+2]= " << norm_pca[ibody*4+2] << " norm_pca[ibody*4+3]= " << norm_pca[ibody*4+3] << std::endl;
				//std::cerr << " body_rot= " << body_rot << " body_tilt= " << body_tilt << " body_psi= " << body_psi << std::endl;
				//std::cerr << " XX(body_offset_3d)= " << XX(body_offset_3d) << " YY(body_offset_3d)= " << YY(body_offset_3d) << " ZZ(body_offset_3d)= " << ZZ(body_offset_3d) << std::endl;

				Matrix2D<RFLOAT> Aresi,  Abody;
				// Aresi is the residual orientation for this ibody
				Euler_angles2matrix(body_rot, body_tilt, body_psi, Aresi);
				// Only apply the residual orientation now!!!
				Abody = (model.orient_bodies[ibody]).transpose() * A_rot90 * Aresi * model.orient_bodies[ibody];

				// Also put back at the centre-of-mass of this body
				Abody.resize(4,4);
				MAT_ELEM(Abody, 0, 3) = XX(body_offset_3d);
				MAT_ELEM(Abody, 1, 3) = YY(body_offset_3d);
				MAT_ELEM(Abody, 2, 3) = ZZ(body_offset_3d);
				MAT_ELEM(Abody, 3, 3) = 1.;

				Mbody.resize(model.Iref[ibody]);
				Mmask.resize(model.masks_bodies[ibody]);
				applyGeometry(model.Iref[ibody], Mbody, Abody, IS_NOT_INV, DONT_WRAP);
				applyGeometry(model.masks_bodies[ibody], Mmask, Abody, IS_NOT_INV, DONT_WRAP);

				img() += Mbody * Mmask;
				sumw += Mmask;
			}

			// Divide the img by sumw to deal with overlapping bodies: just take average
			FOR_ALL_DIRECT_ELEMENTS_IN_MULTIDIMARRAY(img())
			{
				if (DIRECT_MULTIDIM_ELEM(sumw, n) > 1.)
					DIRECT_MULTIDIM_ELEM(img(), n) /= DIRECT_MULTIDIM_ELEM(sumw, n);
			}

			// Write the image to disk
			FileName fn_img = fn_out + "_component" + integerToString(k+1, 3) + "_bin" + integerToString(ibin+1, 3) + ".mrc";
			img.setSamplingRateInHeader(model.pixel_size);
			img.write(fn_img);

		} // end loop ibin

	} // end loop components
}

void FlexAnalyser::writeAllPCAProjections(std::vector< std::vector<double> > &projected_input)
{
	FileName fnt = fn_out+"_projections_along_eigenvectors_all_particles.txt";
	std::ofstream  fh;
	fh.open((fnt).c_str(), std::ios::out);
	if (!fh)
		REPORT_ERROR( (std::string)" FlexAnalyser::writeAllPCAProjections: cannot write to file: " + fnt);

	for (long int ipart = 0; ipart < projected_input.size(); ipart++)
	{
		data.MDimg.getValue(EMDL_IMAGE_NAME, fnt, ipart);
		fh << fnt << " ";
		for (int ival = 0; ival < projected_input[ipart].size(); ival++)
		{
			fh.width(15);
			fh << projected_input[ipart][ival];

		}
		fh << " \n";
	}

	fh.close();
}

void FlexAnalyser::outputSelectedParticles(std::vector< std::vector<double> > &projected_input)
{
	if (select_eigenvalue <= 0)
		return;

	MetaDataTable MDo;
	for (long int ipart = 0; ipart < projected_input.size(); ipart++)
	{
		if (projected_input[ipart][select_eigenvalue-1] > select_eigenvalue_min &&
		    projected_input[ipart][select_eigenvalue-1] < select_eigenvalue_max)
			MDo.addObject(data.MDimg.getObject(ipart));
	}

	int min = ROUND(select_eigenvalue_min);
	int max = ROUND(select_eigenvalue_max);
	FileName fnt = fn_out+"_eval"+integerToString(select_eigenvalue,3)+"_select";
	if (min > -99998)
		fnt += "_min"+integerToString(min);
	if (max < 99998)
		fnt += "_max"+integerToString(max);
	fnt += ".star";
	data.obsModel.save(MDo, fnt, "particles");
	std::cout << " Written out " << MDo.numberOfObjects() << " selected particles in " << fnt << std::endl;
}

void principalComponentsAnalysis(const std::vector< std::vector<double> > &input,
                                 std::vector< std::vector<double> > &eigenvec,
                                 std::vector<double> &eigenval, std::vector<double> &means,
                                 std::vector< std::vector<double> > &projected_input)
{
	std:: cout << "Calculating PCA ..." << std::endl;

	std::vector<std::vector<double> > a;
	long int datasize = input.size();
	if (datasize == 0)
		REPORT_ERROR("ERROR: empty input vector for PCA!");

	// The dimension (n)
	long int n = input[0].size();
	a.resize(n);

	//Get the mean and variance of the given cluster of vectors
	for (int k = 0; k < n; k++)
	{
		a[k].resize(n);
		double sum = 0.0;
		double nn = 0.;
		for (long int i = 0; i < datasize; i++)
		{
			sum += input[i][k];
			nn += 1.0;
		}
		means.push_back(sum / nn);
	}

	for (int i = 0; i < n;i++)
	{
		for (int j = 0;j <= i; j++)
		{
			double sum = 0.0;
			double nn = 0.;
			for (long int k = 0; k < datasize; k++)
			{
				double d1 = input[k][i] - means[i];
				double d2 = input[k][j] - means[j];
				sum += d1 * d2;
				nn += 1.0;
			}
			if (nn > 0.)
				a[i][j] = a[j][i] = sum / nn;
			else
				a[i][j] = a[j][i] = 0;
		}
	}

	eigenval.resize(n);
	eigenvec.resize(n);

	std::vector<double> b;
	b.resize(n);
	std::vector<double> z;
	z.resize(n);
	std::vector<double> &d = eigenval;
	std::vector< std::vector<double> > &v = eigenvec;

	for (int i = 0; i < n; i++)
	{
		v[i].resize(n);
		v[i][i] = 1.0;
		b[i] = d[i] = a[i][i];
	}

	int nrot = 0;

	// Jacobi method (it=iteration number)
	for (int it = 1; it <= 50; it++)
	{
		double threshold;
		double sm = 0.0;
		for (int ip = 0; ip < n - 1; ip++)
		{
			for (int iq = ip + 1; iq < n; iq++)
				sm += fabs(a[iq][ip]);
		}
		if (sm == 0.0)
		{//Done. Sort vectors
			for (int i = 0; i < n - 1; i++)
			{
				int k = i;
				double p = d[i];

				for (int j = i + 1; j < n; j++)
					if (d[j] >= p)
						p = d[k = j];

				if (k != i)
				{//Swap i<->k
					d[k] = d[i];
					d[i] = p;
					std::vector<double> t = v[i];
					v[i] = v[k];
					v[k] = t;
				}
			}

			// Done with PCA now!
			// Just project all data onto the PCA now and exit
			projected_input = input;
			for (int i = 0; i < n; i++)
			{
				for (long int z = 0; z < datasize; z++)
				{
					double cum = 0;
					for (int j = 0; j < n; j++)
						cum += v[i][j] * (input[z][j] - means[j]);
					projected_input[z][i] = cum;
				}  // z
			} // i

			return;
		}

		if (it < 4)
			threshold = 0.2 * sm / (n * n);
		else
			threshold = 0;
		for (int ip = 0; ip < n - 1; ip++)
		{
			for (int iq = ip + 1; iq < n; iq++)
			{
				double g = 100.0 * fabs(a[iq][ip]);
				if (it > 4
				    && fabs(d[ip]) + g == fabs(d[ip])
				    && fabs(d[iq]) + g == fabs(d[iq]))
					a[iq][ip] = 0.0;
				else if (fabs(a[iq][ip]) > threshold)
				{
					double tau, t, s, c;
					double h = d[iq] - d[ip];
					if (fabs(h) + g == fabs(h))
						t = a[iq][ip] / h;
					else
					{
						double theta = 0.5 * h / a[iq][ip];
						t = 1.0 / (fabs(theta) + sqrt(1.0 + theta * theta));
						if (theta < 0.0)
							t = -t;
					}
					c = 1.0 / sqrt(1 + t * t);
					s = t * c;
					tau = s / (1.0 + c);
					h = t * a[iq][ip];
					z[ip] -= h;
					z[iq] += h;
					d[ip] -= h;
					d[iq] += h;
					a[iq][ip] = 0.0;

#define rotate(a,i,j,k,l) \
    g = a[i][j]; \
    h = a[k][l]; \
    a[i][j] = g - s *(h + g*tau); \
    a[k][l] = h + s*(g - h*tau);

					for (int j = 0; j < ip; j++)
					{
						rotate(a, ip, j, iq, j)
					}
					for (int j = ip + 1; j < iq; j++)
					{
						rotate(a, j, ip, iq, j)
					}
					for (int j = iq + 1; j < n; j++)
					{
						rotate(a, j, ip, j, iq)
					}
					for (int j = 0; j < n; j++)
					{
						rotate(v, ip, j, iq, j)
					}

					nrot += 1;
				}//if
			}//for iq
		}//for ip

		for (int ip = 0; ip < n; ip++)
		{
			b[ip] += z[ip];
			d[ip] = b[ip];
			z[ip] = 0.0;
		}

	}//for it

	REPORT_ERROR("ERROR: too many Jacobi iterations in PCA calculation...");
}

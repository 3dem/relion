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

#include <unistd.h>
#include <string.h>

#include <src/image.h>
#include <src/backprojector.h>
#include <src/funcs.h>
#include <src/memory.h>
#include <src/euler.h>
#include <src/time.h>
#include <src/jaz/single_particle/complex_io.h>



class ext_recons_parameters
{
	public:
   	FileName fn_star, fn_recons, fn_data_real, fn_data_imag, fn_weight, fn_out;
   	MultidimArray<RFLOAT> tau2;
   	RFLOAT tau2_fudge;
   	float padding_factor;
   	int ori_size, current_size, ref_dim;
   	int verb;
   	bool skip_gridding, do_map;


	void read(int argc, char **argv)
	{
		if (argc < 2)
		{
			REPORT_ERROR("  Usage: relion_external_reconstruct input.star");
		}
		FileName fn_star = argv[1];
		if (fn_star.getExtension() != "star")
		{
			REPORT_ERROR(" ERROR: input argument does not have a .star extension.");
		}
		skip_gridding  = checkParameter(argc, argv, "--skip_gridding");
		do_map = !checkParameter(argc, argv, "--no_map");
		fn_out =  (checkParameter(argc, argv, "--o")) ? getParameter(argc, argv, "--o") : "";

		MetaDataTable MDlist, MDtau;
		MDlist.read(fn_star, "external_reconstruct_general");
		MDlist.getValue(EMDL_OPTIMISER_EXTERNAL_RECONS_DATA_REAL, fn_data_real);
		MDlist.getValue(EMDL_OPTIMISER_EXTERNAL_RECONS_DATA_IMAG, fn_data_imag);
		MDlist.getValue(EMDL_OPTIMISER_EXTERNAL_RECONS_WEIGHT, fn_weight);
		MDlist.getValue(EMDL_OPTIMISER_EXTERNAL_RECONS_RESULT, fn_recons);
		MDlist.getValue(EMDL_MLMODEL_TAU2_FUDGE_FACTOR, tau2_fudge);
		MDlist.getValue(EMDL_MLMODEL_PADDING_FACTOR, padding_factor);
		MDlist.getValue(EMDL_MLMODEL_DIMENSIONALITY, ref_dim);
		MDlist.getValue(EMDL_MLMODEL_ORIGINAL_SIZE, ori_size);
		MDlist.getValue(EMDL_MLMODEL_CURRENT_SIZE, current_size);

		if (fn_out != "") fn_recons = fn_out;

		MDtau.read(fn_star, "external_reconstruct_tau2");
		tau2.resize(MDtau.numberOfObjects());
		int idx = 0;
		FOR_ALL_OBJECTS_IN_METADATA_TABLE(MDtau)
		{
			RFLOAT mytau2;
			MDtau.getValue(EMDL_MLMODEL_TAU2_REF, mytau2);
			DIRECT_MULTIDIM_ELEM(tau2, idx) = mytau2;
			idx++;
		}

	}

	void reconstruct()
	{
		BackProjector BP(ori_size, ref_dim, "C1", TRILINEAR, padding_factor);
		BP.initZeros(current_size);

		if (skip_gridding) BP.skip_gridding = skip_gridding;

		Image<Complex> Idata;
		Image<RFLOAT> Iweight;
		std::string fn_ext = "."+fn_data_real.getExtension();
		std::string fn_root = fn_data_real.beforeFirstOf("_real");
		ComplexIO::read(Idata, fn_root, fn_ext);
		Iweight.read(fn_weight);

		// Could there be a 1-pixel different in size? use FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM to be safe
		const int r_max = current_size / 2;
		const int r_max2 = ROUND(r_max * padding_factor) * ROUND(r_max * padding_factor);
		FOR_ALL_ELEMENTS_IN_FFTW_TRANSFORM(Idata())
		{
			if (kp*kp + ip*ip + jp*jp < r_max2)
			{
				A3D_ELEM(BP.data, kp, ip, jp) = DIRECT_A3D_ELEM(Idata(), k, i, j);
				A3D_ELEM(BP.weight, kp, ip, jp) = DIRECT_A3D_ELEM(Iweight(), k, i, j);
			}
		}

		BP.reconstruct(Iweight(), 10, do_map, tau2, tau2_fudge);
		Iweight.write(fn_recons);
	}

};


int main(int argc, char *argv[])
{
	ext_recons_parameters prm;

	try
	{
		prm.read(argc, argv);

		prm.reconstruct();
	}
	catch (RelionError XE)
	{
		std::cerr << XE;
		return RELION_EXIT_FAILURE;
	}

	return RELION_EXIT_SUCCESS;
}

#include <src/args.h>
#include <src/jaz/tomography/imod_import.h>
#include <src/jaz/tomography/tomo_ctf_helper.h>
#include <src/jaz/tomography/tomogram_set.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/util/log.h>
#include <src/jaz/util/zio.h>
#include <src/ctf.h>

enum CtfSource
{
	None = 0,
	CtfFind = 1,
	CtfPlotter = 2
};

void checkMissingMandatoryLabel(
		std::vector<EMDLabel> labelAlternatives,
		MetaDataTable& table,
		std::vector<std::vector<EMDLabel>> missingLabels)
{
	for (int i = 0; i < labelAlternatives.size(); i++)
	{
		if (table.containsLabel(labelAlternatives[i]))
		{
			return;
		}
	}

	missingLabels.push_back(labelAlternatives);
}

void handleMissingOptionalLabel(
		EMDLabel label,
		MetaDataTable& table,
		const std::string& argumentName,
		const std::string& argumentValue,
		std::vector<std::pair<EMDLabel,std::string>> missingLabels)
{
	if (!table.containsLabel(label))
	{
		if (argumentValue == "")
		{
			missingLabels.push_back(std::make_pair(label, argumentName));
		}
		else
		{
			for (int i = 0; i < table.numberOfObjects(); i++)
			{
				table.setValueFromString(label, argumentValue, i);
			}
		}
	}
}

int main(int argc, char *argv[])
{
	try
	{
		IOParser parser;

		double voltage, Cs, Q0, hand, pixelSize;
		std::string inStarFn, inTomoSet, outFn, orderFn, fractionalDoseStr, rootDir;
		ImodImport global_IMOD_import;

		try
		{
			parser.setCommandLine(argc, argv);


			int gen_section = parser.addSection("General options");

			inStarFn = parser.getOption("--i", "Input STAR file containing per-tomogram arguments");
			rootDir = parser.getOption("--root", "Root directory containing all IMOD directories", "");
			inTomoSet = parser.getOption("--t", "Input tomogram set", "");
			outFn = parser.getOption("--o", "Output tomogram set");
			hand = textToDouble(parser.getOption("--hand", "Handedness of the tilt geometry", "-1"));


			int optics_section = parser.addSection("Optics options");

			pixelSize = textToDouble(parser.getOption("--angpix", "Pixel size in Å"));
			voltage = textToDouble(parser.getOption("--voltage", "Voltage in kV", "300"));
			Cs = textToDouble(parser.getOption("--Cs", "Spherical aberration in mm", "2.7"));
			Q0 = textToDouble(parser.getOption("--Q0", "Amplitude contrast", "0.1"));


			int dose_section = parser.addSection("Electron dose options");

			orderFn = parser.getOption("--ol", "Frame-order list", "");
			fractionalDoseStr = parser.getOption("--fd", "Fractional dose in electrons per Å²", "");


			int imod_section = parser.addSection("IMOD import options");

			global_IMOD_import.newstComFn = parser.getOption("--nc", "Input command to IMOD's newstack", "newst.com");
			global_IMOD_import.tltComFn = parser.getOption("--tc", "Input command to IMOD's tilt", "tilt.com");

			global_IMOD_import.thicknessCmd = textToDouble(parser.getOption(
				"--thick",
				"Thickness of original IMOD tomogram (overrides the value in tilt.com)",
				"-1.0"));

			global_IMOD_import.offset3Dx = textToDouble(parser.getOption("--offx", "3D offset, X", "0.0"));
			global_IMOD_import.offset3Dy = textToDouble(parser.getOption("--offy", "3D offset, Y", "0.0"));
			global_IMOD_import.offset3Dz = textToDouble(parser.getOption("--offz", "3D offset, Z", "0.0"));

			global_IMOD_import.flipYZ = parser.checkOption("--flipYZ", "Interchange the Y and Z coordinates");
			global_IMOD_import.flipZ = parser.checkOption("--flipZ", "Change the sign of the Z coordinate");
			global_IMOD_import.flipAngles = parser.checkOption("--flipAng", "Change the sign of all tilt angles");

			global_IMOD_import.ali = parser.checkOption("--ali", "Map to aligned stack (.ali)");
			global_IMOD_import.aliSize = parser.checkOption("--ali_size", "Use the size indicated in newst.com");

			Log::readParams(parser);

			if (parser.checkForErrors())
			{
				exit(RELION_EXIT_FAILURE);
			}
		}
		catch (RelionError XE)
		{
			parser.writeUsage(std::cout);
			std::cerr << XE;
			exit(RELION_EXIT_FAILURE);
		}


		rootDir = ZIO::ensureEndingSlash(rootDir);

		MetaDataTable perTomoArguments;
		perTomoArguments.read(inStarFn);


		bool inputsOkay = true;

		std::vector<std::vector<EMDLabel>> missingMandatoryLabels;

		checkMissingMandatoryLabel({EMDL_TOMO_IMPORT_IMOD_DIR}, perTomoArguments, missingMandatoryLabels);
		checkMissingMandatoryLabel({EMDL_TOMO_TILT_SERIES_NAME}, perTomoArguments, missingMandatoryLabels);
		checkMissingMandatoryLabel({EMDL_TOMO_IMPORT_CTFFIND_FILE, EMDL_TOMO_IMPORT_CTFPLOTTER_FILE}, perTomoArguments, missingMandatoryLabels);
		//checkMissingMandatoryLabel({EMDL_TOMO_IMPORT_CULLED_FILE}, perTomoArguments, missingMandatoryLabels);

		if (missingMandatoryLabels.size() == 1)
		{
			std::vector<EMDLabel> missing = missingMandatoryLabels[0];

			if (missing.size() == 1)
			{
				std::cerr
					<< "The following label needs to be present in the input STAR file: "
					<< EMDL::label2Str(missing[0]) << std::endl;
			}
			else
			{
				std::cerr
					<< "One of the following labels needs to be present in the input STAR file: "
					<< EMDL::label2Str(missing[0]) << " or " << EMDL::label2Str(missing[1]) << std::endl;
			}
		}
		else if (missingMandatoryLabels.size() > 1)
		{
			std::cerr << "The following labels need to be present in the input star file:" << std::endl;

			for (int i = 0; i < missingMandatoryLabels.size(); i++)
			{
				std::vector<EMDLabel> missing = missingMandatoryLabels[i];

				if (missing.size() == 1)
				{
					std::cerr << " - " << EMDL::label2Str(missing[0]) << std::endl;
				}
				else
				{
					std::cerr << " - either "
						<< EMDL::label2Str(missing[0]) << " or "
						<< EMDL::label2Str(missing[1]) << std::endl;
				}
			}
		}

		inputsOkay = missingMandatoryLabels.size() == 0;

		std::vector<std::pair<EMDLabel,std::string>> missingOptionalLabels;

		handleMissingOptionalLabel(
			EMDL_TOMO_IMPORT_ORDER_LIST, perTomoArguments,
			"--ol", orderFn, missingOptionalLabels);

		handleMissingOptionalLabel(
			EMDL_TOMO_IMPORT_FRACT_DOSE, perTomoArguments,
			"--fd", fractionalDoseStr, missingOptionalLabels);

		if (missingOptionalLabels.size() == 1)
		{
			std::pair<EMDLabel,std::string> missing = missingOptionalLabels[0];

			std::cerr
				<< "The following property needs to be either present in the input STAR file or supplied as a command-line argument: "
				<< EMDL::label2Str(missing.first) << " or " << missing.second << std::endl;
		}
		else if (missingOptionalLabels.size() > 1)
		{
			std::cerr << "The following properties need to be either present in the input STAR file or supplied as command-line arguments: " << std::endl;

			for (int i = 0; i < missingOptionalLabels.size(); i++)
			{
				std::pair<EMDLabel,std::string> missing = missingOptionalLabels[i];

				std::cerr
					<< " - "<< EMDL::label2Str(missing.first) << " or " << missing.second << std::endl;
			}
		}

		inputsOkay = inputsOkay && missingOptionalLabels.size() == 0;

		if (!inputsOkay)
		{
			return RELION_EXIT_FAILURE;
		}


		if (perTomoArguments.containsLabel(EMDL_TOMO_IMPORT_CTFFIND_FILE) &&
			perTomoArguments.containsLabel(EMDL_TOMO_IMPORT_CTFPLOTTER_FILE))
		{
			std::cerr
				<< "The CTF can only be imported from CTFFind ("
				<< EMDL::label2Str(EMDL_TOMO_IMPORT_CTFFIND_FILE) << ")"
				<< " or ctfplotter ("
				<< EMDL::label2Str(EMDL_TOMO_IMPORT_CTFPLOTTER_FILE) << ")"
				<< "), not both.";

			return RELION_EXIT_FAILURE;
		}


		TomogramSet tomograms;

		if (inTomoSet != "")
		{
			Log::print("Appending new tomograms to "+inTomoSet);
			tomograms = TomogramSet(inTomoSet);
		}


		Log::print(ZIO::itoa(perTomoArguments.numberOfObjects()) + " tomograms to be imported");

		for (int tomo_index = 0; tomo_index < perTomoArguments.numberOfObjects(); tomo_index++)
		{
			CtfSource ctfSource;

			Log::beginSection("Tomogram " + ZIO::itoa(tomo_index + 1) + " / " + ZIO::itoa(perTomoArguments.numberOfObjects()));

			std::string ctfFindFn, ctfPlotterFn;
			perTomoArguments.getValue(EMDL_TOMO_IMPORT_CTFFIND_FILE, ctfFindFn, tomo_index);
			perTomoArguments.getValue(EMDL_TOMO_IMPORT_CTFPLOTTER_FILE, ctfPlotterFn, tomo_index);


			if (ctfPlotterFn == "")
			{
				ctfSource = CtfFind;
				ctfFindFn = rootDir + ctfFindFn;
				Log::print("The CTF will be read from the CTFFind output file "+ctfFindFn);
			}
			else if (ctfFindFn == "")
			{
				ctfSource = CtfPlotter;
				ctfPlotterFn = rootDir + ctfPlotterFn;
				Log::print("The CTF will be read from the ctfplotter output file "+ctfPlotterFn);
			}

			ImodImport ii = global_IMOD_import;

			ii.inDir = rootDir + perTomoArguments.getString(EMDL_TOMO_IMPORT_IMOD_DIR, tomo_index);

			perTomoArguments.getValue(EMDL_TOMO_IMPORT_OFFSET_X, ii.offset3Dx, tomo_index);
			perTomoArguments.getValue(EMDL_TOMO_IMPORT_OFFSET_Y, ii.offset3Dy, tomo_index);
			perTomoArguments.getValue(EMDL_TOMO_IMPORT_OFFSET_Z, ii.offset3Dz, tomo_index);

			const std::string tsFn0 = rootDir + perTomoArguments.getString(EMDL_TOMO_TILT_SERIES_NAME, tomo_index);
			const std::string orderFn = perTomoArguments.getString(EMDL_TOMO_IMPORT_ORDER_LIST, tomo_index);
			const double fractionalDose = perTomoArguments.getDouble(EMDL_TOMO_IMPORT_FRACT_DOSE, tomo_index);

			ii.tsFn = tsFn0;

			std::string name;

			if (perTomoArguments.containsLabel(EMDL_TOMO_NAME))
			{
				name = perTomoArguments.getString(EMDL_TOMO_NAME, tomo_index);
			}
			else
			{
				name = tsFn0;

				if (name.find_last_of(".") != std::string::npos)
				{
					name = name.substr(0, name.find_last_of("."));
				}
			}

			std::string outFnCrop;
			perTomoArguments.getValue(EMDL_TOMO_IMPORT_CULLED_FILE, outFnCrop, tomo_index);

			std::string tsFn;

			Log::print("Importing frame alignment from "+ii.inDir);


			ImodImport::Mapping mapping = ii.import();

			if (mapping.framesMissing)
			{
				if (outFnCrop == "")
				{
					std::cerr << "According to " << ii.inDir + ii.tltComFn << ", frames have been excluded, "
						<< "but no output filename has been specified for the culled stack ("
						<< EMDL::label2Str(EMDL_TOMO_IMPORT_CULLED_FILE) << ").";

					return RELION_EXIT_FAILURE;
				}

				Log::print("Writing culled frame stack to "+outFnCrop);
				ii.writeCulledStack(mapping.oldFrameIndex, outFnCrop);

				tsFn = outFnCrop;
			}
			else
			{
				tsFn = tsFn0;
			}

			const std::string tltFn = ii.lastTltFn;

			Log::print("Deducing frame sequence from "+orderFn+" and "+tltFn);

			std::vector<std::vector<double>> order = ZIO::readFixedDoublesTable(orderFn, 2, ',');
			std::vector<double> tilts = ZIO::readDoubles(tltFn);

			int fc = mapping.oldFrameIndex.size();

			std::vector<double> cumulativeDose(fc);

			for (int i = 0; i < fc; i++)
			{
				double minDist = 720.0;
				int bestJ = -1;

				// find the tilt angle in *.tlt closest to the angle in the order list
				for (int j = 0; j < order.size(); j++)
				{
					double d = order[j][1] - tilts[mapping.oldFrameIndex[i]];
					double dd = d*d;

					if (dd < minDist)
					{
						bestJ = j;
						minDist = dd;
					}
				}

				cumulativeDose[i] = bestJ * fractionalDose;
			}

			std::vector<CTF> ctfs;
			int fcTilt = tilts.size();

			if (ctfSource == CtfFind)
			{
				ctfs = TomoCtfHelper::loadCtffind4(ctfFindFn, fcTilt, voltage, Cs, Q0);
			}
			else
			{
				ctfs = TomoCtfHelper::loadCtfplotter(ctfPlotterFn, fcTilt, voltage, Cs, Q0);
			}

			if (mapping.framesMissing)
			{
				for (int i = 0; i < fc; i++)
					ctfs[i] = ctfs[mapping.oldFrameIndex[i]];
			}

			std::string opticsGroupName = "opticsGroup1";

			perTomoArguments.getValue(EMDL_IMAGE_OPTICS_GROUP_NAME, opticsGroupName, tomo_index);


			tomograms.addTomogram(
				name, tsFn,
				mapping.projections,
				mapping.w, mapping.h, mapping.d,
				cumulativeDose, fractionalDose,
				ctfs, hand, pixelSize, opticsGroupName);

			Log::endSection();
		}

		tomograms.write(outFn);

		Log::print("Tomogram set written to " + outFn);
	}
	catch (RelionError e)
	{
		return RELION_EXIT_FAILURE;
	}

	return RELION_EXIT_SUCCESS;
}

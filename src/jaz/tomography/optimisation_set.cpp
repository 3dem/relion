#include "optimisation_set.h"


void OptimisationSet::read(
		IOParser &parser,        bool showOptimisationSet,
		bool showParticles,      bool particlesMandatory,
		bool showTomograms,      bool tomogramsMandatory,
		bool showTrajectories,   bool trajectoriesMandatory,
		bool showManifolds,      bool manifoldsMandatory,
		bool showReferenceMap,   bool referenceMandatory)
{
	int glob_section = parser.addSection("Complete optimisation set");

	std::string optimisationSet = "";

	if (showOptimisationSet)
	{
		optimisationSet = parser.getOption("--i", "Optimisation set", "");
	}

	if (showParticles || showTomograms || showTrajectories || showManifolds || showReferenceMap)
	{
		int indiv_section = parser.addSection("Individual optimisation components (overriding the optimisation set)");
	}

	if (showParticles)
	{
		particles = parser.getOption("--p", "Particle set", "");
	}

	if (showTomograms)
	{
		tomograms = parser.getOption("--t", "Tomogram set", "");
	}

	if (showTrajectories)
	{
		trajectories = parser.getOption("--mot", "Particle trajectories set", "");
	}

	if (showManifolds)
	{
		manifolds = parser.getOption("--man", "Manifolds set", "");
	}

	if (showReferenceMap)
	{
		refMap1 = parser.getOption("--ref1", "Reference map, half 1", "");
		refMap2 = parser.getOption("--ref2", "Reference map, half 2", "");
		refMask = parser.getOption("--mask", "Reference mask", "");
		refFSC = parser.getOption("--fsc", "Star file containing the FSC of the reference", "");

		freqCutoff_A =  textToDouble(parser.getOption("--freq_cutoff", "Explicit cutoff frequency (in Å; negative to turn off)", "-1"));
		fscThresholdWidth = textToDouble(parser.getOption("--fsc_thresh_width", "Width of the frq. weight flank", "5"));
		flatWeight = !parser.checkOption("--use_SNR_weight", "Weight each shell proportionally to its reference-map confidence");
	}

	if (optimisationSet != "")
	{
		MetaDataTable table;
		table.read(optimisationSet);

		if (particles == "")
		{
			particles = readFromFile(
				table, EMDL_TOMO_PARTICLES_FILE_NAME, "particles",
				optimisationSet, particlesMandatory, parser);
		}

		if (tomograms == "")
		{
			tomograms = readFromFile(
				table, EMDL_TOMO_TOMOGRAMS_FILE_NAME, "tomograms",
				optimisationSet, tomogramsMandatory, parser);
		}

		if (trajectories == "")
		{
			trajectories = readFromFile(
				table, EMDL_TOMO_TRAJECTORIES_FILE_NAME, "trajectories",
				optimisationSet, trajectoriesMandatory, parser);
		}

		if (manifolds == "")
		{
			manifolds = readFromFile(
				table, EMDL_TOMO_MANIFOLDS_FILE_NAME, "manifolds",
				optimisationSet, manifoldsMandatory, parser);
		}

		if (refMap1 == "")
		{
			refMap1 = readFromFile(
					table, EMDL_TOMO_REFERENCE_MAP_1_FILE_NAME, "", "", false, parser);
		}

		if (refMap2 == "")
		{
			refMap2 = readFromFile(
					table, EMDL_TOMO_REFERENCE_MAP_2_FILE_NAME, "", "", false, parser);
		}

		if (refMask == "")
		{
			refMask = readFromFile(
					table, EMDL_TOMO_REFERENCE_MASK_FILE_NAME, "", "", false, parser);
		}

		if (refFSC == "")
		{
			refFSC = readFromFile(
					table, EMDL_TOMO_REFERENCE_FSC_FILE_NAME, "", "", false, parser);
		}
		if (referenceMandatory && (refMap1 == "" || refMap2 == ""))
		{
			reportError(
				"ERROR: The two reference maps were neither found in the optimisation-set file "
				+optimisationSet+" nor specified on the command line.", parser);
		}
	}
	else
	{
		if (particles == "" && particlesMandatory)
		{
			reportMissing("particles", parser);
		}

		if (tomograms == "" && tomogramsMandatory)
		{
			reportMissing("tomograms", parser);
		}

		if (trajectories == "" && trajectoriesMandatory)
		{
			reportMissing("trajectories", parser);
		}

		if (manifolds == "" && manifoldsMandatory)
		{
			reportMissing("manifolds", parser);
		}

		if (referenceMandatory && (refMap1 == "" || refMap2 == ""))
		{
			reportError(
				"ERROR: The two reference maps were not specified on the command line, and no optimisation-set file has been provided.",
				parser);
		}
	}
}

void OptimisationSet::write(std::string filename)
{
	MetaDataTable table;

	table.setIsList(true);
	table.addObject();

	if (particles != "")
	{
		table.setValue(EMDL_TOMO_PARTICLES_FILE_NAME, particles);
	}

	if (tomograms != "")
	{
		table.setValue(EMDL_TOMO_TOMOGRAMS_FILE_NAME, tomograms);
	}

	if (trajectories != "")
	{
		table.setValue(EMDL_TOMO_TRAJECTORIES_FILE_NAME, trajectories);
	}

	if (manifolds != "")
	{
		table.setValue(EMDL_TOMO_MANIFOLDS_FILE_NAME, manifolds);
	}


	if (refMap1 != "")
	{
		table.setValue(EMDL_TOMO_REFERENCE_MAP_1_FILE_NAME, refMap1);
	}

	if (refMap2 != "")
	{
		table.setValue(EMDL_TOMO_REFERENCE_MAP_2_FILE_NAME, refMap2);
	}

	if (refMask != "")
	{
		table.setValue(EMDL_TOMO_REFERENCE_MASK_FILE_NAME, refMask);
	}

	if (refFSC != "")
	{
		table.setValue(EMDL_TOMO_REFERENCE_FSC_FILE_NAME, refFSC);
	}

	table.write(filename);
}

std::string OptimisationSet::readFromFile(
		MetaDataTable& table, EMDLabel label,
		const std::string& argName,
		const std::string& filename,
		bool mandatory,
		IOParser& parser)
{
	if (table.containsLabel(label))
	{
		return table.getString(label);
	}
	else if (mandatory)
	{
		reportError(
			"ERROR: No "+argName+" file was found in the optimisation-set file "
			+filename+" nor specified on the command line.",
			parser);

		return "";
	}
	else
	{
		return "";
	}
}

void OptimisationSet::reportMissing(
		const std::string& argName,
		IOParser& parser)
{
	reportError(
		"ERROR: Neither a "+argName+ // note: all argNames begin with a consonant
		" file nor an optimisation-set file have been specified on the command line.",
		parser);
}

void OptimisationSet::reportError(const std::string& message, IOParser &parser)
{
	parser.reportError(message);
}

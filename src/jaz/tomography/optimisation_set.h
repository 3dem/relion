#ifndef TOMO_OPTIMISATION_SET_H
#define TOMO_OPTIMISATION_SET_H

#include <src/jaz/tomography/reference_map.h>


class OptimisationSet
{
	public:

		void read(
			IOParser& parser,        bool showOptimisationSet,
			bool showParticles,      bool particlesMandatory,
			bool showTomograms,      bool tomogramsMandatory,
			bool showTrajectories,   bool trajectoriesMandatory,
			bool showManifolds,      bool manifoldsMandatory,
			bool showReferenceMap,   bool referenceMandatory);


			std::string
				particles, tomograms, trajectories, manifolds,
				refMap1, refMap2, refMask, refFSC;

			double fscThresholdWidth, freqCutoff_A;

			bool flatWeight;


		void write(std::string filename);


	private:

		std::string readFromFile(
			MetaDataTable& table, EMDLabel label,
			const std::string& argName,
			const std::string& filename,
			bool mandatory,
			IOParser& parser);

		void reportMissing(const std::string& argName, IOParser& parser);

		void reportError(const std::string& message, IOParser& parser);
};

#endif

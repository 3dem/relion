#ifndef CONVERT_PROJECTIONS_PROGRAM_H
#define CONVERT_PROJECTIONS_PROGRAM_H

#include <string>
#include <vector>
#include <src/jaz/gravis/t4Matrix.h>
#include <src/jaz/gravis/t2Vector.h>

class ConvertProjectionsProgram
{
	public:
		
		ConvertProjectionsProgram(){}
		
			int n_threads;
			std::string newstComFn, tltComFn, outFn, outFnCrop, tsFn;
			double thicknessCmd;
			bool flipYZ, flipZ, flipAngles, ali, aliSize;
			double offset3Dx, offset3Dy, offset3Dz;
			
		void run();
		
		
	protected:
		
		void readNewstCom(
			std::string filename, 
			std::string& origStackFn, std::string& aliFn, std::string& xfFn, 
			int& w_ali, int& h_ali, 
			double& binning_ali);
		
		void readTiltCom(
			std::string filename,
			std::string& aliFn,
			std::string& tltFn,
			double& w_out, double& h_out, 
			double& thickness_out, double& binning_out, 
                        std::vector<bool>& exclude,
                        double& shiftX,
                        double& shiftZ);
		
		std::vector<gravis::d4Matrix> loadInvAffineTransforms(
			std::string xformFile, 
			gravis::d2Vector centerOrig, 
			gravis::d2Vector centerAli, 
			double binning,
			bool forAli);
		
		std::vector<gravis::d4Matrix> loadTiltProjections(
			std::string tiltFile, 
			double centerX, double centerY, 
			bool flipAngles);
};

#endif

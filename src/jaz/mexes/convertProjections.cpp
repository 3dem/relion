#include "include/mex_parser.h"
#include <programs/convert_projections.h>


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	MexParser mp(nlhs, plhs, nrhs, prhs);
	
	ConvertProjectionsProgram cpp;
	
	cpp.newstComFn = mp.getString(MexParser::Right, 0, "Input command to IMOD's newstack (usually newst.com)");	
	cpp.tltComFn = mp.getString(MexParser::Right, 1, "Input command to IMOD's tilt (usually tilt.com)");	
	cpp.outFn = mp.getString(MexParser::Right, 2, "Output filename");	
	cpp.outFnCrop = mp.getString(MexParser::Right, 3, "Output filename for culled stack (with frames excluded)");
	
	cpp.tsFn = mp.getString(MexParser::Right, 4, "Original tilt series (only size needed, unless views have to be excluded)");
	
	cpp.thicknessCmd = mp.getDouble(MexParser::Right, 5, "Thickness of original IMOD tomogram (overrides the value in tilt.com)", -1.0);
	
	cpp.offset3Dx = mp.getDouble(MexParser::Right, 6, "3D offset, X", 0.0);
	cpp.offset3Dy = mp.getDouble(MexParser::Right, 7, "3D offset, Y", 0.0);
	cpp.offset3Dz = mp.getDouble(MexParser::Right, 8, "3D offset, Z", 0.0);
		
	cpp.flipYZ = mp.checkOption(MexParser::Right, 9, "Interchange the Y and Z coordinates", false);
	cpp.flipZ = mp.checkOption(MexParser::Right, 10, "Change the sign of the Z coordinate", false);
	cpp.flipAngles = mp.checkOption(MexParser::Right, 11, "Change the sign of all tilt angles", false);
			
	cpp.ali = mp.checkOption(MexParser::Right, 12, "Map to aligned stack (.ali)", false);
	cpp.aliSize = mp.checkOption(MexParser::Right, 13, "Use the size indicated in newst.com", true);
		
	cpp.n_threads = 1;
	
	
	if (!mp.finalize()) return;
				
	cpp.run();
}

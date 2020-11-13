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

int main(int argc, char *argv[])
{
	IOParser parser;
	
	double fractionalDose, voltage, Cs, Q0, hand, pixelSize;
	std::string outFn, outFnCrop, name, tsFn, orderFn, tltFn, ctfFindFn, ctfPlotterFn;
	CtfSource ctfSource;
	ImodImport ii;

	try
	{
		parser.setCommandLine(argc, argv);
		
		
		int gen_section = parser.addSection("General options");
		
		tsFn = parser.getOption("--ts", "Tilt series file name");
		
		name = parser.getOption("--name", "Tomogram name (default: ts_<N>)", "");
		hand = textToDouble(parser.getOption("--hand", "Handedness of the tilt geometry", "-1"));
		
		outFn = parser.getOption("--io", "Input and output tomogram set");
		outFnCrop = parser.getOption("--oc", 
			"Output filename for culled stack (if frames have been excluded)", "");
		
		
		int optics_section = parser.addSection("Optics options");
		
		ctfFindFn = parser.getOption("--ctffind", "Output file written by CTFFind 4", "");
		ctfPlotterFn = parser.getOption("--ctfplotter", "Output file written by ctfplotter", "");
		pixelSize = textToDouble(parser.getOption("--angpix", "Pixel size in Å"));
		voltage = textToDouble(parser.getOption("--voltage", "Voltage in kV", "300"));
		Cs = textToDouble(parser.getOption("--Cs", "Spherical aberration in mm", "2.7"));
		Q0 = textToDouble(parser.getOption("--Q0", "Amplitude contrast", "0.1"));
		
		
		int dose_section = parser.addSection("Electron dose options");
				
		orderFn = parser.getOption("--ol", "Frame-order list (*.csv)", "order_list.csv");
		tltFn = parser.getOption("--tlt", "Tilt angles list (*.tlt, the one from IMOD is used by default)", "");
		fractionalDose = textToDouble(parser.getOption("--fd", "Fractional dose in electrons per Å²"));
		
		
		int imod_section = parser.addSection("IMOD import options");
		
		ii.inDir = parser.getOption("--id", "Directory in which IMOD was run (assumed to contain all other IMOD files)");
		
		ii.newstComFn = parser.getOption("--nc", "Input command to IMOD's newstack", "newst.com");
		ii.tltComFn = parser.getOption("--tc", "Input command to IMOD's tilt", "tilt.com");
		
		ii.tsFn = tsFn;
		
		ii.thicknessCmd = textToDouble(parser.getOption(
			"--thick", 
			"Thickness of original IMOD tomogram (overrides the value in tilt.com)", 
			"-1.0"));
		
		ii.offset3Dx = textToDouble(parser.getOption("--offx", "3D offset, X", "0.0"));
		ii.offset3Dy = textToDouble(parser.getOption("--offy", "3D offset, Y", "0.0"));
		ii.offset3Dz = textToDouble(parser.getOption("--offz", "3D offset, Z", "0.0"));
		
		ii.flipYZ = parser.checkOption("--flipYZ", "Interchange the Y and Z coordinates");
		ii.flipZ = parser.checkOption("--flipZ", "Change the sign of the Z coordinate");
		ii.flipAngles = parser.checkOption("--flipAng", "Change the sign of all tilt angles");
				
		ii.ali = parser.checkOption("--ali", "Map to aligned stack (.ali)");
		ii.aliSize = parser.checkOption("--ali_size", "Use the size indicated in newst.com");
		
		Log::readParams(parser);
				
		if (parser.checkForErrors())
		{
			exit(1);
		}
	}
	catch (RelionError XE)
	{
		parser.writeUsage(std::cout);
		std::cerr << XE;
		exit(1);
	}
	
	if (ctfFindFn == "" && ctfPlotterFn == "")
	{
		ctfSource = None;
		Log::print("Not considering the CTF");
	}
	else if (ctfPlotterFn == "")
	{
		ctfSource = CtfFind;
		Log::print("The CTF will be read from the CTFFind output file "+ctfFindFn);
	}
	else if (ctfFindFn == "")
	{
		ctfSource = CtfPlotter;
		Log::print("The CTF will be read from the ctfplotter output file "+ctfPlotterFn);
	}
	else
	{
		REPORT_ERROR_STR("The CTF can only be imported from CTFFind (--ctffind)"
		              << " or ctfplotter (--ctfplotter), not both.");
	}
	
	
	Log::print("Importing frame alignment from "+ii.inDir);
	
	ImodImport::Mapping mapping = ii.import();
	
	if (mapping.framesMissing)
	{
		if (outFnCrop == "")
		{
			REPORT_ERROR_STR("According to " << ii.inDir+ii.tltComFn << ", frames have been excluded, "
							 << "but no output filename has been specified for the culled stack (--oc).");
		}
				
		Log::print("Writing culled frame stack to "+outFnCrop);
		ii.writeCulledStack(mapping.oldFrameIndex, outFnCrop);
		
		tsFn = outFnCrop;
	}
	
	if (tltFn == "")
	{
		tltFn = ii.lastTltFn;
	}
	
	Log::print("Deducing frame sequence from "+orderFn+" and "+tltFn);
		
	std::vector<std::vector<double>> order = ZIO::readFixedDoublesTable(orderFn, 2, ',');
	std::vector<double> tilts = ZIO::readDoubles(tltFn);
	
	int fc = tilts.size();
	
	std::vector<double> cumulativeDose(fc);
	
	for (int i = 0; i < fc; i++)
	{
		double minDist = 720.0;
		int bestJ = -1;
		
		// find the tilt angle in *.tlt closest to the angle in the order list
		for (int j = 0; j < order.size(); j++)
		{
			double d = order[j][1] - tilts[i];
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
	
	if (ctfSource == CtfFind)
	{
		ctfs = TomoCtfHelper::loadCtffind4(ctfFindFn, fc, voltage, Cs, Q0);
	}
	else
	{
		ctfs = TomoCtfHelper::loadCtfplotter(ctfPlotterFn, fc, voltage, Cs, Q0);
	}
	
	TomogramSet tomograms(outFn);
	
	int my_index = tomograms.size();
	
	if (name == "")
	{
		name = "ts_"+ZIO::itoa(my_index);
	}
	
	tomograms.addTomogram(
		name, tsFn,
		mapping.projections, 
		mapping.w, mapping.h, mapping.d, 
		cumulativeDose, fractionalDose,
		ctfs, hand, pixelSize);
	
	tomograms.write(outFn);
	
	return 0;	
}

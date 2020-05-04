#ifndef TOMOLIST_H
#define TOMOLIST_H

#include <string>
#include <vector>

/*
 will contain: 

	-	tomo name
	-	original tilt series
	-	projections file
	-	ctf file
	-	dose file
	-	fiducial-erased tilt series (optional, can be generated)
	-	file with fiducial positions in 3D (optional)
*/

class TomoList
{
	public:
		
		TomoList();
		TomoList(std::string filename);
		
		
			std::vector<std::vector<std::string>> table;
			
			
		int addTomogram(
				std::string ts, 
				std::string proj, 
				std::string ctf = "noctf", 
				std::string dose = "nodose");
		
		std::string getTiltSeriesFilename(long int t) const;
		std::string getProjectionsFilename(long int t) const;
		std::string getCtfFilename(long int t) const;
		std::string getDoseFilename(long int t) const;
		
		void setTiltSeriesFilename(long int t, std::string s);
		void setProjectionsFilename(long int t, std::string s);
		void setCtfFilename(long int t, std::string s);
		void setDoseFilename(long int t, std::string s);
		
		bool isKnown(long int t);		
		std::vector<long int> getKnownIndices();
		
		void write(std::string filename);
};

#endif

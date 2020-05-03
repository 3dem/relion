#include <src/jaz/util/log.h>
#include <unistd.h>
#include <src/args.h>


int main(int argc, char *argv[])
{
	
	IOParser parser;	
	parser.setCommandLine(argc, argv);
	
	Log::readParams(parser);
	
	const int N = 70;
	
	Log::beginSection("initialising");
	
		Log::beginSection("component X");
		
			Log::print("loading X");
			usleep(920000u);
			
			Log::beginProgress("transmogrifying", N);
			
				for (int i = 0; i < N; i++)
				{
					Log::updateProgress(i);
					
					usleep(50000u);
				}
					
			Log::endProgress(); 
			
			Log::print("component X ready");
			
		Log::endSection();
		
	Log::endSection();
	
	Log::beginSection("processing");
		
		Log::beginSection("component X");
		
			Log::beginSection("preparing");
			usleep(230000u);
			
				Log::print("step 1");
				usleep(400000u);
				Log::print("step 2");
				usleep(270000u);
				
			Log::endSection();
			
		Log::endSection();
	
	Log::endSection();
	
	Log::beginProgress("performing global optimisation", N);
	
	for (int i = 0; i < N; i++)
	{
		Log::updateProgress(i);
		
		usleep(70000u);
	}
			
	Log::endProgress(); 
	
	Log::print("writing");
	
	usleep(370000u);
	
	Log::print("all done");
	
	return 0;
}

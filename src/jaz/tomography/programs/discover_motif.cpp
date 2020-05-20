#include "discover_motif.h"

#include <src/jaz/image/filter.h>
#include <src/jaz/image/normalization.h>
#include <src/jaz/image/structure_tensor.h>
#include <src/jaz/image/local_extrema.h>
#include <src/jaz/tomography/projection/projection.h>
#include <src/jaz/tomography/lattice/point_cloud.h>
#include <src/jaz/tomography/lattice/motif.h>
#include <src/jaz/tomography/lattice/motif_refinement.h>
#include <src/jaz/single_particle/vtk_helper.h>
#include <src/jaz/mesh/mesh.h>
#include <src/jaz/mesh/mesh_builder.h>

#include <omp.h>


void DiscoverMotifProgram::run()
{
	std::istringstream sts(pointsFn);
	std::string token;
	
	std::vector<std::string> pointsFns;
	
	while(std::getline(sts, token, ','))
	{
		pointsFns.push_back(token);
	}
	
	const int cloudNum = pointsFns.size();
	
	std::vector<PointCloud> clouds(cloudNum);
	
	for (int cl = 0; cl < cloudNum; cl++)
	{
		clouds[cl].readCsv(pointsFns[cl]);
		
		if (doDistHist)
		{
			std::vector<size_t> distHist = clouds[cl].distanceHist(dhBins, maxDistHist);
			
			std::ofstream dhf(outFn+"_dist_density.dat");
			
			double sum = 0.0;
			
			for (int b = 1; b < dhBins; b++)
			{
				sum += distHist[b];
			}
					
			const double scale = 1.0 / sum;
			
			for (int b = 1; b < dhBins; b++)
			{
				const double x0 = ((double)b       ) * maxDistHist / (double)dhBins;
				const double x1 = ((double)b + 1.0 ) * maxDistHist / (double)dhBins;
				
				const double vol = x1*x1*x1 - x0*x0*x0;
				const double nrm = scale * distHist[b] / vol;
				
				dhf << x0 << " " << nrm << "\n";
				dhf << x1 << " " << nrm << "\n";
			}
			
			return;
		}
		
		clouds[cl].findNeighbors(spacingMin, spacingMax);
		
		double avgN(0.0);
		
		for (int p = 0; p < clouds[cl].vertices.size(); p++)
		{
			avgN += clouds[cl].neighbors[p].size();
		}
		
		avgN /= clouds[cl].vertices.size();
		
		std::cout << "    cloud #" << cl << ", average neighbor count: " << avgN << std::endl;
		
		if (doNcHist)
		{
			std::vector<int> ncHist(maxNcHist);
			
			for (int p = 0; p < clouds[cl].vertices.size(); p++)
			{
				const int n = clouds[cl].neighbors[p].size();
				
				if (n < maxNcHist)
				{
					ncHist[n]++;
				}
			}
			
			std::ofstream nhf(outFn+"_neighbors_hist.dat");
			
			nhf << 0 << " " << ncHist[0] << "\n";
			nhf << 0.5 << " " << ncHist[0] << "\n";
			
			for (int i = 1; i < ncHist.size(); i++)
			{
				nhf << (i - 0.5) << " " << ncHist[i] << "\n";
				nhf << (i + 0.5) << " " << ncHist[i] << "\n";			
			}
			
			return;
		}
	}
	
	
	double t0 = omp_get_wtime();
		
	std::vector<std::pair<Motif, double>> motifs = Motif::discover(
                clouds, neighborCount, motIters, evalIters, radius, motifCount, 
				refineDiscovery, num_threads);
	
	double t1 = omp_get_wtime();
	
	std::cout << "discovery completed in " << (t1 - t0) << " sec." << std::endl;
	
	for (int i = 0; i < motifs.size(); i++)
	{
		Motif motif = motifs[i].first;
		double score = motifs[i].second;
		
		Mesh mesh0 = motif.toMesh();
		
		std::stringstream sts;
		sts << i;
		
		mesh0.writeObj(outFn + "_motif_" + sts.str() + ".obj");
				
		std::vector<MotifDetection> occurences 
				= motif.detectIn(clouds, minSc, radius, refineDetection);
		
		std::cout << i << ": " << occurences.size() 
				<< " occurences found with a score above " << minSc << "\n";
				
		std::vector<Mesh> occMesh = motif.visualize(occurences);
		
		for (int i = 0; i < occMesh.size(); i++)
		{
			std::stringstream sts2;
			sts2 << i;
			occMesh[i].writePly(outFn + "_detections_" + sts.str() + "_tomo_" + sts2.str() + ".ply");
		}
		
		
		std::pair<Motif, std::vector<MotifDetection>> refdets 
				= motif.refine(clouds, occurences, radius, reg);
		
		
		std::vector<Mesh> occMeshRef = refdets.first.visualize(refdets.second);
		
		for (int i = 0; i < occMeshRef.size(); i++)
		{
			std::stringstream sts2;
			sts2 << i;
			occMeshRef[i].writePly(outFn + "_detections-refined_" + sts.str() + "_tomo_" + sts2.str() + ".ply");
		}
		
		Mesh mesh1 = refdets.first.toMesh();
		mesh1.writeObj(outFn + "_motif-ref_" + sts.str() + ".obj");
		
		
		Catalogue cat = MotifDetection::toCatalogue(refdets.second, cloudBinning);
		
		cat.write(outFn + "_particles_" + sts.str() + ".cat");
					
	}
}

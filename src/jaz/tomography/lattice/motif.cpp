#include "motif.h"
#include "motif_alignment.h"
#include "motif_refinement.h"
#include <src/jaz/mesh/mesh_builder.h>
#include <src/jaz/optimization/lbfgs.h>
#include <src/jaz/optimization/gradient_descent.h>
#include <src/jaz/optimization/nelder_mead.h>
#include <omp.h>

using namespace gravis;

Motif::Motif()
{	
}

Motif::Motif(const std::vector<d3Vector>& points)
:	points(points)
{	
}

std::vector<MotifDetection> Motif::detectIn(
        const std::vector<PointCloud>& clouds,
        double minScore,
        double tolerance,
        bool refine) const
{
	const size_t tc = clouds.size();
	
	std::vector<MotifDetection> out;
	out.reserve(tc * clouds[0].vertices.size());
	
	for (int t = 0; t < tc; t++)
	{
		const size_t cpc = clouds[t].vertices.size();
		const size_t mpc = points.size() - 1;
				
		std::vector<d3Vector> motifPts(mpc);
		
		for (size_t i = 0; i < mpc; i++)
		{
			motifPts[i] = points[i+1] - points[0];
		}
		
		for (size_t j = 0; j < cpc; j++)
		{
			const size_t nc = clouds[t].neighbors[j].size();
			
			if (nc < 2) continue;
			
			std::vector<d3Vector> testPts = clouds[t].getNeighborhood(j);
			
			const double offset = 1.0;
			std::pair<d4Matrix, double> alignment = alignPointsRansac(motifPts, testPts, tolerance, offset);
	
			if (refine)
			{
				alignment = refineAlignment(motifPts, testPts, tolerance, alignment.first);
			}
	
			const double score = evaluateAlignment(motifPts, testPts, tolerance, alignment.first) - 1.0;
			
			if (score > minScore) 
			{
				d4Matrix T = alignment.first;
				T(0,3) += clouds[t].vertices[j].x;
				T(1,3) += clouds[t].vertices[j].y;
				T(2,3) += clouds[t].vertices[j].z;
				
				out.push_back(MotifDetection(T, j, t, score));
			}
		}
	}
	
	return out;
}

std::vector<Mesh> Motif::visualize(const std::vector<MotifDetection>& occurences) const
{
	Mesh mesh0 = toMesh(0.5);
	Mesh mesh1 = mesh0;
	
	int tc = 0;
	
	for (int i = 0; i < occurences.size(); i++)
	{
		const MotifDetection& md = occurences[i];
		const int minTc = md.tomoIndex + 1;
		if (minTc > tc) tc = minTc;
	}
	
	std::vector<Mesh> out(tc);
	
	for (int i = 0; i < occurences.size(); i++)
	{
		const MotifDetection& md = occurences[i];
		
		for (int v = 0; v < mesh0.vertices.size(); v++)
		{
			mesh1.vertices[v] = (md.alignment * d4Vector(mesh0.vertices[v])).toVector3();
		}
		
		const int t = md.tomoIndex;
		
		MeshBuilder::insert(mesh1, out[t], dRGB(md.score, 0.0, 0.0));
	}
	
	return out;
}

Mesh Motif::toMesh(double length, double pointedness) const
{
	Mesh out;
	
	const d3Vector p = points[0];
	
	for (int i = 1; i < points.size(); i++)
	{
		const d3Vector q = points[i];
		
		if (pointedness == 0.0)
		{
			MeshBuilder::addBar(p, p + length*(q - p), 0.1, 5, out);
		}
		else
		{
			MeshBuilder::addPointyBar(p, p + length*(q - p), 17, pointedness, 5, out);
		}
	}
	
	return out;
}

std::pair<Motif, std::vector<MotifDetection>> Motif::refine(
	const std::vector<PointCloud>& clouds, 
	const std::vector<MotifDetection>& detections, 
	double tolerance,
	double reg)
{
	MotifRefinement mr(*this, clouds, detections, tolerance, reg);
	
	std::vector<double> x0 = mr.getInitialParams();
	
	std::vector<double> x1 = LBFGS::optimize(x0, mr, 1, 20000, 1e-4, 1e-16);
	//std::vector<double> x1 = GradientDescent::optimize(x0, mr, 0.01, 0.0001, 0.00001, 20000, 0.0, true);
	//std::vector<double> x1 = NelderMead::optimize(x0, mr, 0.01, 1e-5, 20000, 1.0, 2.0, 0.5, 0.5, true);
	
	return mr.getDetections(x1);
}

std::vector<std::pair<Motif, double>> Motif::discover(
	const std::vector<PointCloud>& clouds, 
	int neighbors, 
	int motifsNum, 
	int evalsNum,
	double tolerance,
    int number,
    bool refine,
	int num_threads)
{
	//const int pc = cloud.vertices.size();
	const int tc = clouds.size();
	int pc = 0;
	
	for (int t = 0; t < tc; t++)
	{
		const size_t cpc = clouds[t].vertices.size();
		pc += cpc;
	}
		
	std::vector<int> pt2cloud(pc), pt2index(pc);
	size_t curr_index = 0;
	
	std::vector<int> goodPts(0);
	goodPts.reserve(pc);
	
	for (int t = 0; t < tc; t++)
	{
		const size_t cpc = clouds[t].vertices.size();
	
		for (int p = 0; p < cpc; p++)
		{
			pt2cloud[curr_index + p] = t;
			pt2index[curr_index + p] = p;
					
			if (clouds[t].neighbors[p].size() == neighbors)
			{
				goodPts.push_back(curr_index + p);
			}
		}
		
		curr_index += cpc;
	}
		
	const int gpc = goodPts.size();

	std::cout << "Found " << gpc << " points with " << neighbors << " neighbors." << std::endl;
	
	int mc;
	bool tryAllMotifs;
	
	if (gpc < motifsNum)
	{
		std::cerr << "Warning: unable to try " << motifsNum << " motifs - only " << gpc
				  << " points with " << neighbors << " neighbors present. Trying all of them.\n";
		
		mc = gpc;
		tryAllMotifs = true;
	}
	else
	{
		mc = motifsNum;
		tryAllMotifs = false;
	}
	
	int ec;
	bool tryAllPts;
	
	if (pc < evalsNum)
	{
		std::cerr << "Warning: unable to evaluate motifs on " << evalsNum << " points - only " << pc
				  << " points present. Trying all of them.\n";
		
		tryAllPts = true;
		ec = pc;
	}
	else
	{
		tryAllPts = false;
		ec = evalsNum;
	}
	
	std::vector<int> bestMotifs(number, -1);
	std::vector<int> bestMotifsCloud(number, -1);
	std::vector<double> bestScores(number, 0.0);
	
	std::vector<bool> tried(gpc, false);
	
	#pragma omp parallel for num_threads(num_threads)
	for (int m = 0; m < mc; m++)
	{
		const int ig = tryAllMotifs? m : rand() % gpc;
		
		if (tried[ig]) continue;
		else tried[ig] = true;
		
		const int glob_i = goodPts[ig];		
		const int cloud_i = pt2cloud[glob_i];
		const int index_i = pt2index[glob_i];
		
		std::vector<d3Vector> motifPts = clouds[cloud_i].getNeighborhood(index_i);
		
		double score = 0.0;
		
		for (int e = 0; e < ec; e++)
		{
			const int glob_j = tryAllPts? e : rand() % pc;
			
			if (glob_j == glob_i) continue;
			
			const int cloud_j = pt2cloud[glob_j];
			const int index_j = pt2index[glob_j];
			
			const int qnc = clouds[cloud_j].neighbors[index_j].size();
			
			if (qnc < 2) continue;
			
			std::vector<d3Vector> testPts = clouds[cloud_j].getNeighborhood(index_j);
			
			std::pair<d4Matrix, double> alignment = alignPointsRansac(motifPts, testPts, tolerance, 1.0);
	
			if (refine)
			{
				std::pair<d4Matrix, double> alignment2 = refineAlignment(motifPts, testPts, tolerance, alignment.first);
				score += alignment2.second;
			}
			else
			{
				score += alignment.second;
			}
		}
		
		if (score < bestScores[number-1]) continue;
		
		#pragma omp critical
		{
			for (int n = 0; n < number; n++)
			{
				if (bestScores[n] < score)
				{
					for (int nn = number-1; nn > n; nn--)
					{
						bestMotifs[nn] = bestMotifs[nn-1];
						bestMotifsCloud[nn] = bestMotifsCloud[nn-1];
						bestScores[nn] = bestScores[nn-1];
					}
					
					bestMotifs[n] = index_i;
					bestMotifsCloud[n] = cloud_i;
					bestScores[n] = score;
					
					std::cout << "  ";
					
					for (int nn = 0; nn < number; nn++)
					{
						std::cout << bestScores[nn] << " ";
					}
					
					std::cout << std::endl;
					
					break;
				}
			}
		}
	}

	int triedNum(0);

	for (int i = 0; i < tried.size(); i++)
	{
		if (tried[i]) triedNum++;
	}

	std::cout << triedNum << " points tried." << std::endl;
	
	std::vector<std::pair<Motif, double>> out;
	
	for (int i = 0; i < number; i++)
	{
		const int index = bestMotifs[i];
		const int tomo = bestMotifsCloud[i];
        const double score = bestScores[i];
		
		const PointCloud& cl = clouds[tomo];
		
		const int nc = cl.neighbors[index].size();
		
		std::cout << nc << " neighbors: ";
		
		Motif m;
		m.points.resize(nc + 1);
		
		m.points[0] = d3Vector(0.0, 0.0, 0.0);
		
		for (int n = 0; n < nc; n++)
		{
			m.points[n+1] = cl.vertices[cl.neighbors[index][n]] - cl.vertices[index];
			std::cout << m.points[n+1] << " ";
		}
		std::cout << "\n";
		
		out.push_back(std::make_pair(m, score));
	}
	
	return out;
}

std::pair<d4Matrix, double> Motif::alignPointsRansac(
		const std::vector<d3Vector>& pts_a,
		const std::vector<d3Vector>& pts_b,
		double tolerance, double offset)
{
	const int pc_a = pts_a.size();
	const int pc_b = pts_b.size();
	
	const double r2 = tolerance * tolerance;
	
	d4Matrix bestP;
	double maxScore = 0.0;
	
	for (int a0 = 0; a0 < pc_a; a0++)
	{
		const double l_a0 = pts_a[a0].length();
		const d3Vector va0 = pts_a[a0] / l_a0;
		
		for (int a1 = 0; a1 < pc_a; a1++)
		{
			if (a1 == a0) continue;
			
			const double l_a1 = pts_a[a1].length();
			const double a1_on_a0 = pts_a[a1].dot(va0);
						
			const d3Vector va1 = pts_a[a1].cross(pts_a[a0]).normalize();
			const d3Vector va2 = pts_a[a0].cross(va1).normalize();
					
			for (int b0 = 0; b0 < pc_b; b0++)
			{
				const double l_b0 = pts_b[b0].length();
				const d3Vector vb0 = pts_b[b0] / l_b0;
				
				const double dl0 = l_b0 - l_a0;
				
				if (std::abs(dl0) > tolerance) continue;
					
				for (int b1 = 0; b1 < pc_b; b1++)
				{
					if (b1 == b0) continue;
					
					const double l_b1 = pts_b[b1].length();
					
					const double dl1 = l_b1 - l_a1;
					
					if (std::abs(dl1) > tolerance) continue;
					
					const double b1_on_b0 = pts_b[b1].dot(vb0);
					
					const double d1on0 = b1_on_b0 - a1_on_a0;
					
					if (std::abs(d1on0) > tolerance) continue;
					
					const d3Vector vb1 = pts_b[b1].cross(pts_b[b0]).normalize();
					const d3Vector vb2 = pts_b[b0].cross(vb1).normalize();
					
					const d3Matrix At(
							va0.x, va0.y, va0.z,
							va1.x, va1.y, va1.z,
							va2.x, va2.y, va2.z);
					
					const d3Matrix B(
							vb0.x, vb1.x, vb2.x,
							vb0.y, vb1.y, vb2.y,
							vb0.z, vb1.z, vb2.z);
					
					const d3Matrix P = B * At;
					
					double score = -offset;
					
					for (int aa = 0; aa < pc_a; aa++)
					{
						const d3Vector aab = P * pts_a[aa];
						
						for (int bb = 0; bb < pc_b; bb++)
						{
							const double d2 = (pts_b[bb] - aab).norm2();
							score += r2 / (r2 + d2);
						}
					}
					
					if (score > maxScore)
					{
						maxScore = score;
						bestP = P;
					}
				}
			}
		}
	}
			
	return std::make_pair(bestP, maxScore);
}

std::pair<d4Matrix, double> Motif::refineAlignment(
		const std::vector<d3Vector> &pts_a, 
		const std::vector<d3Vector> &pts_b, 
		double tolerance, 
		const d4Matrix initialRot)
{
	std::vector<d3Vector> motifPts(pts_a.size() + 1);
	motifPts[0] = d3Vector(0.0);
	
	for (int i = 0; i < pts_a.size(); i++)
	{
        motifPts[i+1] = pts_a[i];
	}
			
	std::vector<d3Vector> targetPts(pts_b.size() + 1);
	targetPts[0] = d3Vector(0.0);
	
	for (int i = 0; i < pts_b.size(); i++)
	{
        targetPts[i+1] = pts_b[i];
	}
	
	
	MotifAlignment ma(motifPts, targetPts, tolerance);
	std::vector<double> par0 = MotifAlignment::getParams(initialRot);
	
	std::vector<double> par1 = LBFGS::optimize(par0, ma, 1, 20000, 1e-5, 1e-20);
	//std::vector<double> par1 = GradientDescent::optimize(par0, ma, 0.01, 0.0001, 0.00001, 20000);
	//std::vector<double> par1 = NelderMead::optimize(par0, ma, 0.01, 1e-5, 20000, 1.0, 2.0, 0.5, 0.5, false);
	
	d4Matrix A = MotifAlignment::getMatrix(par1);
	
    return std::make_pair(A, -ma.f(par1, 0));
}

double Motif::evaluateAlignment(const std::vector<d3Vector> &pts_a, const std::vector<d3Vector> &pts_b, double tolerance, const d4Matrix A)
{
    std::vector<d3Vector> motifPts(pts_a.size() + 1);
    motifPts[0] = d3Vector(0.0);

    for (int i = 0; i < pts_a.size(); i++)
    {
        motifPts[i+1] = pts_a[i];
    }

    std::vector<d3Vector> targetPts(pts_b.size() + 1);
    targetPts[0] = d3Vector(0.0);

    for (int i = 0; i < pts_b.size(); i++)
    {
        targetPts[i+1] = pts_b[i];
    }


    MotifAlignment ma(motifPts, targetPts, tolerance);
    std::vector<double> par0 = MotifAlignment::getParams(A);

    return -ma.f(par0, 0);
}

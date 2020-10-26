#include "motif_refinement.h"
#include "motif.h"
#include "motif_detection.h"
#include "motif_alignment.h"
#include "point_cloud.h"
#include <src/jaz/math/Euler_angles_dynamo.h>


using namespace gravis;


MotifRefinement::MotifRefinement(
		const Motif& motif,
		const std::vector<PointCloud>& clouds, 
		const std::vector<MotifDetection>& detections, 
		double tolerance, 
		double reg)

:	mc(motif.points.size()),
	dc(detections.size()),
	tolerance(tolerance),
	reg(reg)
{
	int tpc(0);
	
	initialMotifPts = motif.points;
	
	initialCenters.resize(dc);
	initialAngles.resize(dc);
	tomoInds.resize(dc);
	
	for (int d = 0; d < dc; d++)
	{
		const MotifDetection& dd = detections[d];
		const d4Matrix& A = dd.alignment;
		
		const d3Vector ang = EulerDynamo::matrixToAngles(A);
		
		initialCenters[d] = d3Vector(A(0,3), A(1,3), A(2,3));
		initialAngles[d] = d3Vector(ang[0], ang[1], ang[2]);
		tomoInds[d] = dd.tomoIndex;
		
		tpc += clouds[dd.tomoIndex].neighbors[dd.centerIndex].size();
	}
	
	targetPts.resize(tpc);
	
	int ind(0);
	
	neighborInd.resize(dc);
	neighborNum.resize(dc);
	detectionInds.resize(dc);
	
	for (int d = 0; d < dc; d++)
	{
		const MotifDetection& dd = detections[d];
		const int nc = clouds[dd.tomoIndex].neighbors[dd.centerIndex].size();
		
		for (int n = 0; n < nc; n++)
		{
			const int ni = clouds[dd.tomoIndex].neighbors[dd.centerIndex][n];
			targetPts[ind + n] = clouds[dd.tomoIndex].vertices[ni];
		}
		
		neighborInd[d] = ind;
		neighborNum[d] = nc;
		detectionInds[d] = dd.centerIndex;
		
		ind += nc;
	}
}

double MotifRefinement::f(
		const std::vector<double> &x, 
		void *tempStorage) const
{
	std::vector<d3Vector> mpos(mc), t(dc);
	std::vector<d3Matrix> R(dc);
	
	mpos[0] = d3Vector(0.0);
	
	for (int m = 1; m < mc; m++)
	for (int i = 0; i < 3; i++)
	{
		mpos[m][i] = x[3*(m-1) + i];
	}
	
	for (int d = 0; d < dc; d++)
	{
		d3Vector ang;
		
		for (int i = 0; i < 3; i++)
		{
			t[d][i] = x[3*(mc-1) + 6*d + i];
			ang[i] = x[3*(mc-1) + 6*d + 3 + i];
		}
		
		R[d] = EulerDynamo::anglesToMatrix3(ang[0], ang[1], ang[2]);
	}
	
	const double r2 = tolerance * tolerance;
		
	double score(0.0);
	
	for (int d = 0; d < dc; d++)
	{
		const size_t n0 = neighborInd[d];
		const size_t nc = neighborNum[d];
		
		for (int m = 0; m < mc; m++)
		{
			const d3Vector sp = R[d] * mpos[m] + t[d];
			
			for (size_t n = 0; n < nc; n++)
			{
				const d3Vector sq = targetPts[n0 + n];
				
				const double d2 = (sp - sq).norm2();
				score += r2 / (r2 + d2);
			}
		}
	}
	
	for (size_t m = 1; m < mc; m++)
	{
		double d2 = (mpos[m] - initialMotifPts[m]).norm2();		
		score -= dc * reg * d2;
	}
	
	return -score;
}

void MotifRefinement::grad(
		const std::vector<double> &x, 
		std::vector<double> &gradDest, 
		void *tempStorage) const
{
	std::vector<d3Vector> mpos(mc), t(dc);
	std::vector<t4Vector<d3Matrix>> RdR(dc, t4Vector<d3Matrix>(d3Matrix(),d3Matrix(),d3Matrix(),d3Matrix()));
	
	mpos[0] = d3Vector(0.0);
	
	for (int m = 1; m < mc; m++)
	for (int i = 0; i < 3; i++)
	{
		mpos[m][i] = x[3*(m-1) + i];
	}
	
	const size_t det0 = 3 * (mc-1);
	
	for (int d = 0; d < dc; d++)
	{
		d3Vector ang;
		
		for (int i = 0; i < 3; i++)
		{
			t[d][i] = x[det0 + 6*d + i];
			ang[i] = x[det0 + 6*d + 3 + i];
		}
		
		RdR[d] = EulerDynamo::anglesToMatrixAndDerivatives(ang[0], ang[1], ang[2]);
	}
	
	const double r2 = tolerance * tolerance;
	
	for (int i = 0; i < gradDest.size(); i++)
	{
		gradDest[i] = 0.0;
	}
	
	for (int d = 0; d < dc; d++)
	{
		const size_t n0 = neighborInd[d];
		const size_t nc = neighborNum[d];
		
		for (int m = 0; m < mc; m++)
		{
			const d3Matrix& dR_dPhi   = RdR[d].x;
			const d3Matrix& dR_dTheta = RdR[d].y;
			const d3Matrix& dR_dChi   = RdR[d].z;
			
			const d3Matrix& R = RdR[d].w;
			
			
			const d3Vector sp = R * mpos[m] + t[d];
			
			const d3Vector dSp_dPhi   = dR_dPhi * mpos[m];
			const d3Vector dSp_dTheta = dR_dTheta * mpos[m];
			const d3Vector dSp_dChi   = dR_dChi * mpos[m];
			
			const d3Vector dSp_dSx(R(0,0), R(1,0), R(2,0));
			const d3Vector dSp_dSy(R(0,1), R(1,1), R(2,1));
			const d3Vector dSp_dSz(R(0,2), R(1,2), R(2,2));
			
			
			for (size_t n = 0; n < nc; n++)
			{
				const d3Vector sq = targetPts[n0 + n];
				
				const double d2 = (sp - sq).norm2();
				
				//score += r2 / (r2 + d2);
				
				const d3Vector dD2_dSp = 2.0 * (sp - sq);
				const double dScore_dD2 = -r2 / ((r2 + d2)*(r2 + d2));			
				
				const d3Vector dScore_dSp = dScore_dD2 * dD2_dSp;
				
				gradDest[det0 + 6*d + 0] -= dScore_dSp.x;
				gradDest[det0 + 6*d + 1] -= dScore_dSp.y;
				gradDest[det0 + 6*d + 2] -= dScore_dSp.z;
				
				gradDest[det0 + 6*d + 3] -= dScore_dSp.dot(dSp_dPhi);
				gradDest[det0 + 6*d + 4] -= dScore_dSp.dot(dSp_dTheta);
				gradDest[det0 + 6*d + 5] -= dScore_dSp.dot(dSp_dChi);
				
				if (m > 0)
				{
					gradDest[3*(m-1) + 0] -= dScore_dSp.dot(dSp_dSx);
					gradDest[3*(m-1) + 1] -= dScore_dSp.dot(dSp_dSy);
					gradDest[3*(m-1) + 2] -= dScore_dSp.dot(dSp_dSz);
				}
			}
		}
	}
	
	for (size_t m = 1; m < mc; m++)
	{
		const d3Vector del = mpos[m] - initialMotifPts[m];		
		
		//score -= dc * reg * d2;
		
		const d3Vector dD2_dm = 2.0 * del;
		const double dScore_dD2 = -dc * reg;
		
		const d3Vector dScore_dm = dScore_dD2 * dD2_dm;
		
		gradDest[3*(m-1) + 0] -= dScore_dm.x;
		gradDest[3*(m-1) + 1] -= dScore_dm.y;
		gradDest[3*(m-1) + 2] -= dScore_dm.z;
	}
}

std::vector<double> MotifRefinement::getInitialParams()
{
	std::vector<double> out(3*(mc-1) + 6*dc);
	const int det0 = 3*(mc-1);
	
	for (int m = 1; m < mc; m++)
	for (int i = 0; i < 3; i++)
	{
		out[3*(m-1) + i] = initialMotifPts[m][i];
	}
	
	for (int d = 0; d < dc; d++)
	for (int i = 0; i < 3; i++)
	{
		out[det0 + 6*d + i] = initialCenters[d][i];
		out[det0 + 6*d + 3 + i] = initialAngles[d][i];
	}
	
	return out;
}

std::pair<Motif, std::vector<MotifDetection>> MotifRefinement::getDetections(const std::vector<double> &x)
{
	std::vector<MotifDetection> out(dc);
	
	std::vector<d3Vector> mpos(mc);
	
	mpos[0] = d3Vector(0.0);
	
	for (int m = 1; m < mc; m++)
	for (int i = 0; i < 3; i++)
	{
		mpos[m][i] = x[3*(m-1) + i];
	}
	
	Motif motif(mpos);
	
	for (int d = 0; d < dc; d++)
	{
		d3Vector ang;
		
		for (int i = 0; i < 3; i++)
		{
			ang[i] = x[3*(mc-1) + 6*d + 3 + i];
		}
		
		d4Matrix R = EulerDynamo::anglesToMatrix4(ang[0], ang[1], ang[2]);
		
		for (int i = 0; i < 3; i++)
		{
			R(i,3) = x[3*(mc-1) + 6*d + i];
		}
		
		std::vector<d3Vector> dpos(neighborNum[d]);
		const size_t n0 = neighborInd[d];

		for (int n = 0; n < neighborNum[d]; n++)
		{
			dpos[n] = targetPts[n0 + n];
		}
					
		MotifAlignment ma(mpos, dpos, tolerance);
		std::vector<double> par0 = MotifAlignment::getParams(R);
				
		const double score = -ma.f(par0, 0);
					
		out[d] = MotifDetection(R, detectionInds[d], tomoInds[d], score);
	}
	
	return std::make_pair(motif, out);
}

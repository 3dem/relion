#include "motif_alignment.h"
#include <src/jaz/math/Euler_angles_dynamo.h>

using namespace gravis;

#define TX 0
#define TY 1
#define TZ 2
#define PHI 3
#define THETA 4
#define CHI 5


MotifAlignment::MotifAlignment(
		const std::vector<gravis::d3Vector> &motifPts, 
		const std::vector<gravis::d3Vector> &targetPts, 
		double tolerance)
:	
	motifPts(motifPts),
	targetPts(targetPts),
	tolerance(tolerance)
{	
}

double MotifAlignment::f(
		const std::vector<double> &x, void *tempStorage) const
{
	const d3Vector t(x[TX],x[TY],x[TZ]);
	const d3Matrix R = EulerDynamo::anglesToMatrix3(x[PHI], x[THETA], x[CHI]);
	
	const double r2 = tolerance * tolerance;
		
	double score(0.0);
	
	for (int p = 0; p < motifPts.size(); p++)
	{
		const d3Vector sp = R * motifPts[p] + t;
		
		for (int q = 0; q < targetPts.size(); q++)
		{
			const d3Vector sq = targetPts[q];
			
			const double d2 = (sp - sq).norm2();
			score += r2 / (r2 + d2);
		}
	}
	
	return -score;
}

void MotifAlignment::grad(
		const std::vector<double> &x, 
		std::vector<double> &gradDest, 
		void *tempStorage) const
{
	const d3Vector t(x[TX],x[TY],x[TZ]);
	
	const t4Vector<d3Matrix> RdR = EulerDynamo::anglesToMatrixAndDerivatives(x[PHI], x[THETA], x[CHI]);
	
	const d3Matrix R = RdR.w;
	
	const d3Matrix dR_dPhi = RdR.x;
	const d3Matrix dR_dTheta = RdR.y;
	const d3Matrix dR_dChi = RdR.z;
	
	const double r2 = tolerance * tolerance;
		
	for (int i = 0; i < gradDest.size(); i++)
	{
		gradDest[i] = 0.0;
	}
	
	for (int p = 0; p < motifPts.size(); p++)
	{
		const d3Vector sp = R * motifPts[p] + t;
		
		const d3Vector dSp_dPhi   = dR_dPhi * motifPts[p];
		const d3Vector dSp_dTheta = dR_dTheta * motifPts[p];
		const d3Vector dSp_dChi   = dR_dChi * motifPts[p];
		
		for (int q = 0; q < targetPts.size(); q++)
		{
			const d3Vector sq = targetPts[q];
			
			const double d2 = (sp - sq).norm2();
			//score += r2 / (r2 + d2);
			
			const d3Vector dD2_dSp = 2.0 * (sp - sq);
			const double dScore_dD2 = -r2 / ((r2 + d2)*(r2 + d2));			
			
			const d3Vector dScore_dSp = dScore_dD2 * dD2_dSp;
			
			gradDest[TX] -= dScore_dSp.x;
			gradDest[TY] -= dScore_dSp.y;
			gradDest[TZ] -= dScore_dSp.z;
			
			gradDest[PHI]   -= dScore_dSp.dot(dSp_dPhi);
			gradDest[THETA] -= dScore_dSp.dot(dSp_dTheta);
			gradDest[CHI]   -= dScore_dSp.dot(dSp_dChi);
		}
	}
}

d4Matrix MotifAlignment::getMatrix(const std::vector<double> &x)
{
	d4Matrix out = EulerDynamo::anglesToMatrix4(x[PHI], x[THETA], x[CHI]);
	
	out(0,3) = x[TX];
	out(1,3) = x[TY];
	out(2,3) = x[TZ];
	
	return out;
}

std::vector<double> MotifAlignment::getParams(d4Matrix A)
{
	std::vector<double> out(6);
	
	out[TX] = A(0,3);
	out[TY] = A(1,3);
	out[TZ] = A(2,3);
	
	d3Vector ang = EulerDynamo::matrixToAngles(A);
	
	out[PHI]   = ang[0];
	out[THETA] = ang[1];
	out[CHI]   = ang[2];
	
	return out;
}

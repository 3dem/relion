#ifndef MODULAR_ALIGNMENT_H
#define MODULAR_ALIGNMENT_H

#include <src/jaz/optimization/optimization.h>
#include <src/jaz/image/buffered_image.h>
#include <src/jaz/tomography/particle_set.h>
#include <src/jaz/mesh/mesh_builder.h>
#include <src/jaz/util/log.h>
#include <src/jaz/tomography/motion/trajectory.h>
#include <src/jaz/tomography/particle_set.h>
#include <src/jaz/tomography/extraction.h>
#include <src/jaz/tomography/prediction.h>
#include <src/jaz/tomography/tomogram.h>
#include <src/ctf.h>
#include <src/jaz/image/centering.h>
#include <src/jaz/math/Tait_Bryan_angles.h>
#include <src/jaz/math/Gaussian_process.h>
#include <src/jaz/util/zio.h>
#include <omp.h>

class CTF;

struct ModularAlignmentSettings
{
	bool constParticles, constAngles, constShifts, perFrame2DDeformation;
};

template<class MotionModel, class DeformationModel2D>
class ModularAlignment : public FastDifferentiableOptimization
{
	public:

		ModularAlignment(
				const std::vector<BufferedImage<double>>& CCs,
				const std::vector<gravis::d4Matrix>& frameProj,
				ParticleSet& dataSet,
				const std::vector<ParticleIndex>& partIndices,
				const MotionModel& motionModel,
				const DeformationModel2D& deformationModel2D,
				ModularAlignmentSettings settings,
				const Tomogram& tomogram,
				double paddingFactor,
				int progressBarOffset,
				int num_threads,
				bool verbose);
		

			const MotionModel& motionModel;
			const DeformationModel2D& deformationModel2D;

			const std::vector<BufferedImage<double>>& CCs;  // one frame stack for each particle
			const std::vector<gravis::d4Matrix>& frameProj; // initial projection matrices
			ParticleSet& dataSet;
			const std::vector<ParticleIndex>& partIndices;

			ModularAlignmentSettings settings;

			gravis::d3Vector tomoCentre;
			double frameDose, pixelSize, paddingFactor;
			int progressBarOffset, num_threads;

			bool verbose, devMode;

			int fc, pc, mpc, dc, maxRange;


			std::vector<gravis::d3Vector> initialPos;

			gravis::d4Matrix minusCentre, plusCentre;
			
			

		double gradAndValue(const std::vector<double>& x, std::vector<double>& gradDest) const;
			
		std::vector<gravis::d4Matrix> getProjections(
				const std::vector<double>& x,
				const std::vector<int>& frameSequence) const;

		std::vector<gravis::d3Vector> getParticlePositions(
				const std::vector<double>& x) const;

		Trajectory getTrajectory(
				const std::vector<double>& x,
				int p,
				const std::vector<int>& frameSequence) const;

		std::vector<Trajectory> exportTrajectories(
				const std::vector<double>& x,
				const ParticleSet& dataSet,
				const std::vector<int>& frameSequence) const;
		
		void visualiseTrajectories2D(
				const std::vector<double>& x,
				double scale,
				const std::string& tomo_name,
				const std::string& file_name_root);
		
		void report(int iteration, double cost, const std::vector<double>& x) const;
		

		int getParamCount() const
		{
			const int fs = getFrameStride();
			const int ds = settings.perFrame2DDeformation? dc * fc : dc;
			
			return fs * fc + 3 * pc + mpc * (fc - 1) + ds;
		}


	protected:

		inline int getFrameStride() const
		{
			int out = 0;

			if (!settings.constAngles) out += 3;
			if (!settings.constShifts) out += 2;

			return out;
		}
		
		inline int getAlignmentBlockOffset() const
		{
			return 0;
		}
		
		inline int getPositionsBlockOffset(int fs) const
		{
			return fs * fc;
		}
		
		inline int getMotionBlockOffset(int fs) const
		{
			return fs * fc + 3 * pc;
		}
		
		inline int get2DDeformationsBlockOffset(int fs) const
		{
			return fs * fc + 3 * pc + mpc * (fc - 1);
		}

		inline void readViewParams(
			const std::vector<double>& x, int f,
			double& phi, double& theta, double& psi,
			double& dx, double& dy) const;
};

template<class MotionModel, class DeformationModel2D>
ModularAlignment<MotionModel, DeformationModel2D>::ModularAlignment(
		const std::vector<BufferedImage<double>>& CCs,
		const std::vector<gravis::d4Matrix>& frameProj,
		ParticleSet& dataSet,
		const std::vector<ParticleIndex>& partIndices,
		const MotionModel& motionModel,
		const DeformationModel2D& deformationModel2D,
		ModularAlignmentSettings settings,
		const Tomogram& tomogram,
		double paddingFactor,
		int progressBarOffset,
		int num_threads,
		bool verbose)
:	
	motionModel(motionModel),
	deformationModel2D(deformationModel2D),
	CCs(CCs),
	frameProj(frameProj),
	dataSet(dataSet),
	partIndices(partIndices),
	settings(settings),
	tomoCentre(tomogram.centre),
	frameDose(tomogram.getFrameDose()),
	pixelSize(tomogram.optics.pixelSize),
	paddingFactor(paddingFactor),
	progressBarOffset(progressBarOffset),
	num_threads(num_threads),
	verbose(verbose),
	devMode(false),
	fc(frameProj.size()),
	pc(partIndices.size()),
    mpc(motionModel.getParameterCount()),
    dc(deformationModel2D.getParameterCount()),
	maxRange(CCs[0].xdim / (2 * paddingFactor))
{	
	initialPos.resize(pc);
	
	for (int p = 0; p < pc; p++)
	{
		initialPos[p] = dataSet.getPosition(partIndices[p]);
	}
	
	minusCentre = gravis::d4Matrix(
			1, 0, 0, -tomoCentre.x, 
			0, 1, 0, -tomoCentre.y, 
			0, 0, 1, -tomoCentre.z, 
			0, 0, 0, 1 );
	
	plusCentre = gravis::d4Matrix(
			1, 0, 0, tomoCentre.x, 
			0, 1, 0, tomoCentre.y, 
			0, 0, 1, tomoCentre.z, 
			0, 0, 0, 1 );
}

/*
	Parameter Layout:
	
	0:                  ([phi, theta, psi], [dx, dy]) * fc           frame alignment: fs * fc
	fs * fc:             [dx0, dy0, dz0]  *  pc;                     static part. shifts: 3 * pc
	fs * fc + 3 * pc:    [b0x, b0y, b0z][b1x, ..., bBz] * (fc-1);    motion: mpc * (fc - 1)
	
	fs * fc  +  3 * pc  +  mpc * (fc - 1)   total
*/

template<class MotionModel, class DeformationModel2D>
double ModularAlignment<MotionModel, DeformationModel2D>::gradAndValue(
        const std::vector<double>& x, std::vector<double>& gradDest) const
{
	const int fs = getFrameStride();
	const int xs = x.size();
	const int data_pad = 512;
	const int step_grad = xs + data_pad;
	const int step_frame = fc + data_pad;
	const int pos_block = getPositionsBlockOffset(fs);
	const int mot_block = getMotionBlockOffset(fs);
	const int def_block = get2DDeformationsBlockOffset(fs);

	std::vector<gravis::d4Matrix> P(fc), P_phi(fc), P_theta(fc), P_psi(fc);

	for (int f = 0; f < fc; f++)
	{
		double phi, theta, psi, dx, dy;
		readViewParams(x, f, phi, theta, psi, dx, dy);

		const gravis::d4Matrix Q = 
				TaitBryan::anglesToMatrix4(phi, theta, psi);
		
		gravis::t4Vector<gravis::d3Matrix> dQ =
				TaitBryan::anglesToMatrixAndDerivatives(phi, theta, psi);

		gravis::d4Matrix Q_phi(dQ[0]);
		gravis::d4Matrix Q_theta(dQ[1]);
		gravis::d4Matrix Q_psi(dQ[2]);

		Q_phi(3,3) = 0.0;
		Q_theta(3,3) = 0.0;
		Q_psi(3,3) = 0.0;

		const gravis::d4Matrix centProj = minusCentre * frameProj[f];

		P[f] = plusCentre * Q * centProj;

		P_phi[f]   = plusCentre * Q_phi   * centProj;
		P_theta[f] = plusCentre * Q_theta * centProj;
		P_psi[f]   = plusCentre * Q_psi   * centProj;

		P[f](0,3) += dx;
		P[f](1,3) += dy;
	}


	std::vector<double> grad_par(step_grad * num_threads, 0.0);
	std::vector<double> val_par(data_pad * num_threads, 0.0);
	std::vector<gravis::d3Vector> dC_dPos(step_frame * num_threads, gravis::d3Vector(0.0, 0.0, 0.0));


	#pragma omp parallel for num_threads(num_threads)
	for (int p = 0; p < pc; p++)
	{
		const int th = omp_get_thread_num();

		gravis::d3Vector shift = settings.constParticles?
			gravis::d3Vector(0.0, 0.0, 0.0) :
			gravis::d3Vector(x[pos_block + 3*p], x[pos_block + 3*p+1], x[pos_block + 3*p+2]);

		for (int f = 0; f < fc; f++)
		{
			const gravis::d4Vector pos4(initialPos[p] + shift);

			const gravis::d2Vector p0 = (frameProj[f] * gravis::d4Vector(initialPos[p])).xy();
			const gravis::d2Vector pl = (P[f] * pos4).xy();
						
			const int def_block_f = def_block + (settings.perFrame2DDeformation? f * dc : 0);
			
			gravis::d2Vector def, def_x, def_y;			
			deformationModel2D.computeShiftAndGradient(pl, &x[def_block_f], def, def_x, def_y);
			        
			const gravis::d2Vector p1 = pl + def;
			
			const gravis::d2Vector dp = p1 - p0;

			const double dx_img = (dp.x + maxRange) * paddingFactor;
			const double dy_img = (dp.y + maxRange) * paddingFactor;

			const gravis::d2Vector pl_phi   = (P_phi[f]   * pos4).xy();
			const gravis::d2Vector pl_theta = (P_theta[f] * pos4).xy();
			const gravis::d2Vector pl_psi   = (P_psi[f]   * pos4).xy();
			

			const gravis::d3Vector g0 = -((double)paddingFactor) *
				Interpolation::cubicXYGradAndValue_clip(CCs[p], dx_img, dy_img, f);

			val_par[th*data_pad] += g0.z;
			
			/*
				dp_phi = [U,V] pl_phi
				
				g0^T dp_phi = g0^T [U,V] pl_phi = g^T pl_phi
				
				=> 
				
				g^T = g0^T [U,V] = [<g0,U>, <g0,V>]
			*/
			// alternatively: deformationModel2D.transformImageGradient
			const gravis::d2Vector g(
			            (def_x.x + 1.0) * g0.x + def_x.y * g0.y,
			            def_y.x * g0.x + (def_y.y + 1.0) * g0.y);
			
			deformationModel2D.updateCostGradient(pl, g0.xy(), &x[def_block_f], &gradDest[def_block_f]);


			if (settings.constAngles)
			{
				if (!settings.constShifts)
				{
					grad_par[th*step_grad + fs*f    ]  +=  g.x;
					grad_par[th*step_grad + fs*f + 1]  +=  g.y;
				}
			}
			else
			{
				if (settings.constShifts)
				{
					grad_par[th*step_grad + fs*f    ]  +=  pl_phi.x   * g.x  +  pl_phi.y   * g.y;
					grad_par[th*step_grad + fs*f + 1]  +=  pl_theta.x * g.x  +  pl_theta.y * g.y;
					grad_par[th*step_grad + fs*f + 2]  +=  pl_psi.x   * g.x  +  pl_psi.y   * g.y;
				}
				else
				{
					grad_par[th*step_grad + fs*f    ]  +=  pl_phi.x   * g.x  +  pl_phi.y   * g.y;
					grad_par[th*step_grad + fs*f + 1]  +=  pl_theta.x * g.x  +  pl_theta.y * g.y;
					grad_par[th*step_grad + fs*f + 2]  +=  pl_psi.x   * g.x  +  pl_psi.y   * g.y;
					grad_par[th*step_grad + fs*f + 3]  +=  g.x;
					grad_par[th*step_grad + fs*f + 4]  +=  g.y;
				}
			}

			if (!settings.constParticles)
			{
				const gravis::d3Vector dC_dPos_f(
					P[f](0,0) * g.x  +  P[f](1,0) * g.y,
					P[f](0,1) * g.x  +  P[f](1,1) * g.y,
					P[f](0,2) * g.x  +  P[f](1,2) * g.y);

				dC_dPos[th*step_frame + f] = dC_dPos_f;

				grad_par[th*step_grad + pos_block + 3*p    ]  +=  dC_dPos_f.x;
				grad_par[th*step_grad + pos_block + 3*p + 1]  +=  dC_dPos_f.y;
				grad_par[th*step_grad + pos_block + 3*p + 2]  +=  dC_dPos_f.z;

				if (f < fc-1)
				{
					motionModel.updatePosition(&x[mot_block + f*mpc], p, shift);
				}
			}
		}

		if (!settings.constParticles)
		{
			motionModel.updateCostGradient(
				&dC_dPos[th*step_frame], p, fc, 
				&grad_par[th*step_grad + mot_block]);
		}
	}

	for (int i = 0; i < xs; i++)
	{
		gradDest[i] = 0.0;
	}

	double cost = 0.0;

	for (int th = 0; th < num_threads; th++)
	{
		cost += val_par[th*data_pad];
	}

	for (int th = 0; th < num_threads; th++)
	for (int i = 0; i < xs; i++)
	{
		gradDest[i] += grad_par[th*step_grad + i];
	}

	if (!settings.constParticles)
	{
		cost += motionModel.computePriorCostAndGradient(
			&x[mot_block], fc, &gradDest[mot_block]);
	}

	return cost;
}

template<class MotionModel, class DeformationModel2D>
std::vector<gravis::d4Matrix> ModularAlignment<MotionModel, DeformationModel2D>::getProjections(
		const std::vector<double> &x,
		const std::vector<int>& frameSequence) const
{
	std::vector<gravis::d4Matrix> out(fc);
		
	for (int f = 0; f < fc; f++)
	{	
		double phi, theta, psi, dx, dy;		
		readViewParams(x, f, phi, theta, psi, dx, dy);
		
		const gravis::d4Matrix Q = TaitBryan::anglesToMatrix4(phi, theta, psi);		
		
		const int fa = frameSequence[f];
		
		out[fa] = plusCentre * Q * minusCentre * frameProj[f];
		
		out[fa](0,3) += dx;
		out[fa](1,3) += dy;
	}
		
	return out;
}

template<class MotionModel, class DeformationModel2D>
std::vector<gravis::d3Vector> ModularAlignment<MotionModel, DeformationModel2D>::getParticlePositions(
		const std::vector<double>& x) const
{
	std::vector<gravis::d3Vector> out(pc);
	
	const int fs = getFrameStride();
	const int pos_block = getPositionsBlockOffset(fs);
	
	for (int p = 0; p < pc; p++)
	{
		out[p] = initialPos[p] + gravis::d3Vector(
				x[pos_block + 3*p    ],
				x[pos_block + 3*p + 1],
				x[pos_block + 3*p + 2]);
	}

	return out;
}

template<class MotionModel, class DeformationModel2D>
Trajectory ModularAlignment<MotionModel, DeformationModel2D>::getTrajectory(
		const std::vector<double> &x, int p,
		const std::vector<int>& frameSequence) const
{
	Trajectory out(fc);
	
	if (settings.constParticles) return out;
	
	const int fs = getFrameStride();
	const int mot_block = getMotionBlockOffset(fs);

	gravis::d3Vector shift(0.0, 0.0, 0.0);
	
	for (int f = 0; f < fc; f++)
	{
		const int fa = frameSequence[f];
		
		out.shifts_Ang[fa] = pixelSize * shift;
		
		if (f < fc-1)
		{
			motionModel.updatePosition(&x[mot_block + f*mpc], p, shift);
		}
	}
	
	return out;
}

template<class MotionModel, class DeformationModel2D>
std::vector<Trajectory> ModularAlignment<MotionModel, DeformationModel2D>::exportTrajectories(
		const std::vector<double>& x, 
		const ParticleSet& dataSet,
		const std::vector<int>& frameSequence) const
{
	std::vector<Trajectory> out(pc);

	for (int p = 0; p < pc; p++)
	{
		const int pp = partIndices[p].value;

		out[p] = dataSet.motionTrajectories[pp] + getTrajectory(x, p, frameSequence);
	}

	return out;
}

template<class MotionModel, class DeformationModel2D>
void ModularAlignment<MotionModel, DeformationModel2D>::visualiseTrajectories2D(
        const std::vector<double> &x, 
        double scale, 
        const std::string& tomo_name, 
        const std::string& file_name_root)
{
	std::vector<int> timeSeq(fc);

	for (int f = 0; f < fc; f++)
	{
		timeSeq[f] = f;
	}

	std::vector<std::string> plot_names {"XY", "XZ", "YZ"};
	std::vector<gravis::i2Vector> dim_indices{
		gravis::i2Vector(0,1), 
		gravis::i2Vector(0,2), 
		gravis::i2Vector(1,2)};

	for (int dim = 0; dim < 3; dim++)
	{
		CPlot2D plot2D(tomo_name + " Motion " + plot_names[dim]);
		plot2D.SetXAxisSize(600);
		plot2D.SetYAxisSize(600);
		plot2D.SetDrawLegend(false);
		plot2D.SetFlipY(true);

		for (int p = 0; p < pc; p++)
		{
			Trajectory track = getTrajectory(x, p, timeSeq);

			// Mark start of each track
			CDataSet start;
			start.SetDrawMarker(true);
			start.SetMarkerSize(8);
			start.SetDatasetColor(0.2,0.5,1.0);

			const gravis::d3Vector a = initialPos[p] + scale * track.shifts_Ang[0] / pixelSize;
			const gravis::i2Vector di = dim_indices[dim];

			CDataPoint point3(a[di[0]],a[di[1]]);
			start.AddDataPoint(point3);
			plot2D.AddDataSet(start);
		}

		for (int p = 0; p < pc; p++)
		{
			Trajectory track = getTrajectory(x, p, timeSeq);

			CDataSet curve;
			curve.SetDrawMarker(false);
			curve.SetDatasetColor(0.0,0.0,0.0);
			curve.SetLineWidth(0.5);

			for (int f = 0; f < fc; f++)
			{
				const gravis::d3Vector a = initialPos[p] + scale * track.shifts_Ang[f] / pixelSize;
				const gravis::i2Vector di = dim_indices[dim];

				CDataPoint point(a[di[0]],a[di[1]]);

				curve.AddDataPoint(point);
			}

			plot2D.AddDataSet(curve);
		}

		std::string label_x = plot_names[dim].substr(0,1) +
				" (in pixels; trajectory scaled by " + ZIO::itoa(scale) + ")";

		std::string label_y = plot_names[dim].substr(1,1);

		plot2D.SetXAxisTitle(label_x);
		plot2D.SetYAxisTitle(label_y);

		FileName fn_eps = file_name_root + "_" + plot_names[dim] + ".eps";

		plot2D.OutputPostScriptPlot(fn_eps);
	}
}

template<class MotionModel, class DeformationModel2D>
void ModularAlignment<MotionModel, DeformationModel2D>::report(
		int iteration, double cost, const std::vector<double> &x) const
{
	if (devMode)
	{
		int prec = (int)(log(iteration) / log(10));
		int step = pow(10, prec);

		if (prec == 0 || iteration % step == 0)
		{
			std::cout << iteration << " \t " << std::setw(11) << std::setprecision(10) << cost << std::endl;
		}
	}
	else
	{
		Log::updateProgress(progressBarOffset + iteration);
	}
}



template<class MotionModel, class DeformationModel2D>
inline void ModularAlignment<MotionModel, DeformationModel2D>::readViewParams(
		const std::vector<double>& x, int f,
		double& phi, double& theta, double& psi,
		double& dx, double& dy) const
{
	const int fs = getFrameStride();
	int offset = f * fs;
	
	if (settings.constAngles)
	{
		phi   = 0.0;
		theta = 0.0;
		psi   = 0.0;
	}
	else
	{
		phi   = x[offset  ];
		theta = x[offset+1];
		psi   = x[offset+2];

		offset += 3;
	}

	if (settings.constShifts)
	{
		dx    = 0.0;
		dy    = 0.0;
	}
	else
	{
		dx    = x[offset  ];
		dy    = x[offset+1];

		offset += 2;
	}
}





#endif

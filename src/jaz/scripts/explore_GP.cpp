#include <iostream>
#include <src/Eigen/Dense>
#include <src/Eigen/SVD>
#include <time.h>
#include <src/jaz/math/svd_helper.h>
#include <src/jaz/math/Gaussian_process.h>
#include <src/jaz/tomography/data_set.h>
#include <src/jaz/util/zio.h>
#include <src/jaz/util/log.h>


using namespace Eigen;
using namespace gravis;
 
double diff(timespec start, timespec end);

int main()
{
	{
		std::string catFn = "refined_particles_2.star";
		
		DataSet* dataSet = DataSet::load(catFn, "");	
		std::vector<std::vector<int>> particles = dataSet->splitByTomogram();
		
		const int pc = particles[0].size();
		std::vector<d3Vector> initialPos(pc);
		
		for (int p = 0; p < pc; p++)
		{
			initialPos[p] = dataSet->getPosition(particles[0][p]);
		}
		
		bool sqExpKernel = true;
		double sig_vel_px = 1.111;
		double sig_div_px = 2222.22;
		
		
		std::cout << "    sig_vel_px = " << sig_vel_px << std::endl;
		std::cout << "    sig_div_px = " << sig_div_px << std::endl;
		std::cout << "    #particles = " << pc << std::endl;
		
		GpKernel* kernel(0);
		
		if (sqExpKernel)
		{
			kernel = new SquareExponentialKernel(sig_vel_px, sig_div_px);
		}
		else
		{
			kernel = new ExponentialKernel(sig_vel_px, sig_div_px);
		}
		
		Matrix2D<double> C = GaussianProcess::computeCovariance(initialPos, kernel);
		
		delete kernel;
		
		
		std::cout << "\n    Decomposing covariance matrix...\n";
		
		
		timespec time1, time2;	
		
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
		
		GaussianProcess::Basis defBasis = GaussianProcess::getBasis(C, pc);
		
		clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
		
		const int bc = defBasis.eigenvalues.size();
		
		if (bc == pc)
		{
			Log::print("Keeping all " + ZIO::itoa(bc) + " eigendeformations");
		}
		else
		{
			Log::print("Keeping " + ZIO::itoa(bc) + " out of " + ZIO::itoa(pc) + " eigendeformations");
		}
		
		std::cout << "\n    time: " << diff(time1,time2) << " sec. " << std::endl;
		
		std::cout << "\n  Eigenvalues:" << std::endl;
		
		std::cout
			<< '\n'
			<< "    " << defBasis.eigenvalues[0] << ", \n"
			<< "    " << defBasis.eigenvalues[1] << ", \n"
			<< "    " << defBasis.eigenvalues[2] << ", \n"
			<< "    " << "... , \n"
			<< "    " << defBasis.eigenvalues[pc-1] << std::endl;
		
		std::cout << "\n  Eigenvectors:" << std::endl;
		
		for (int j = 0; j < 3; j++)
		{
			std::cout
				<< '\n'
				<< "    " << defBasis.eigenvectors[j*pc + 0] << ", \n"
				<< "    " << defBasis.eigenvectors[j*pc + 1] << ", \n"
				<< "    " << defBasis.eigenvectors[j*pc + 2] << ", \n"
				<< "    " << "... , \n"
				<< "    " << defBasis.eigenvectors[j*pc + pc-1] << std::endl;
		}
		
		std::ofstream evDat("debug/deformation_eigenvalues.dat");
		
		for (int i = 0; i < defBasis.eigenvalues.size(); i++)
		{
			evDat << i << ' ' << defBasis.eigenvalues[i] << '\n';
		}
		
		return 0;	
		
	}
	
	{
		const int s = 3000;
		
		std::cout << "size = " << s << std::endl;
		
		const int nr = s;
		const int nc = s;
		
		MatrixXd A0 = MatrixXd::Random(nr,nc);
		
		Matrix2D<double> A1(nr,nc);
		
		for (int c = 0; c < nc; c++)
		for (int r = 0; r < nr; r++)
		{
			A1(r,c) = A0(r,c);
		}
		
		/*{
			std::cout << "\nEigen::Jacobi: \n" << std::endl;
			
			timespec time1, time2;		
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
			
			JacobiSVD<MatrixXd> svd(A0, ComputeThinU | ComputeThinV);
			
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
			
			std::cout << "    eigenvalue 1: "<< svd.singularValues()[0] << std::endl;		
			std::cout << "    time: " << diff(time1,time2) << std::endl;
		}*/
		
		{
			std::cout << "\nEigen::BDCSVD: \n" << std::endl;
			
			timespec time1, time2;		
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
			
			BDCSVD<MatrixXd> svd(A0, ComputeFullV);
			
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
			
			std::cout << "    eigenvalue 1: "<< svd.singularValues()[0] << std::endl;		
			std::cout << "    time: " << diff(time1,time2) << std::endl;
		}
		
		{
			std::cout << "\nXMIPP: \n" << std::endl;
			
			Matrix2D<double> U, Vt;
			Matrix1D<double> S;
			
			timespec time1, time2;	
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time1);
			SvdHelper::decompose(A1, U, S, Vt);
			clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &time2);
			
			std::cout << "    eigenvalue 1: " << S(0) << std::endl;		
			std::cout << "    time: " << diff(time1,time2) << std::endl;
		}
		
		std::cout << std::endl;	
	}
	
	return 0;
}

double diff(timespec start, timespec end)
{
    timespec temp;
	
    if (end.tv_nsec < start.tv_nsec) 
	{
        temp.tv_sec = end.tv_sec - start.tv_sec - 1;
        temp.tv_nsec = 1000000000 + end.tv_nsec - start.tv_nsec;
    } 
	else 
	{
        temp.tv_sec = end.tv_sec - start.tv_sec;
        temp.tv_nsec = end.tv_nsec - start.tv_nsec;
    }
	
    return (double)temp.tv_sec + (double)temp.tv_nsec / 1e9;
}

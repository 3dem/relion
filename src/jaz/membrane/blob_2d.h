#ifndef BLOB_2D_H
#define BLOB_2D_H

#include <src/jaz/gravis/t3Vector.h>
#include <src/jaz/tomography/tomogram.h>

#define FIRST_BLOB_FREQUENCY 2

class Blob2D
{
	public:
		
		Blob2D();		
		Blob2D(gravis::d2Vector center, double smoothingRadius = 0.0);
		Blob2D(const std::vector<double>& params, double smoothingRadius = 0.0);
		
		
			gravis::d2Vector center;			
			std::vector<dComplex> amplitudes;
			double smoothingRadius;
			
			std::vector<dComplex> sin_cos_table;
			
			
		std::vector<double> radialAverage(
				const RawImage<float>& frame,
				const RawImage<float>& weight,
				double relevantRadius = -1.0) const;
		
		std::pair<std::vector<double>,std::vector<double>> radialAverageAndWeight(
				const RawImage<float>& frame,
				const RawImage<float>& weight,
				double relevantRadius = -1.0) const;
		
		std::pair<std::vector<double>,std::vector<double>> radialAverageAndWeightInSectors(
				const RawImage<float>& frame,
				const RawImage<float>& weight,
		        int sectors,
				double relevantRadius = -1.0) const;

		double radialAverageError(
				const RawImage<float>& frame,
				const RawImage<float>& weight,
				const std::vector<double>& radAvg) const;

		BufferedImage<float> drawError(
				const RawImage<float>& frame,
				const RawImage<float>& weight,
				const std::vector<double>& radAvg) const;

		BufferedImage<float> radialAverageProjection(
				const RawImage<float>& frame,
				const std::vector<double>& radAvg) const;
		
		void erase(
		        const RawImage<float>& micrographs,
				RawImage<float>& erased_out,
				RawImage<float>& blob_out,
				const RawImage<float>& weight,
				double radius, double taper) const;
		
		void eraseLocally(
		        const RawImage<float>& micrographs,
				RawImage<float>& erased_out,
				RawImage<float>& blob_out,
				const RawImage<float>& weight,
				double radius, 
		        double taper,
		        double smoothness) const;
		
		std::pair<BufferedImage<float>, BufferedImage<float>> 
			transformToPolar(
		        const RawImage<float>& frame,
		        const RawImage<float>& mask,
				double maxRadius) const;
		
		
		inline std::vector<double> toVector() const;
		
		inline int findMaxRadius(gravis::i2Vector imageSize) const;
		
		inline double getOffset(gravis::d2Vector v) const;
		inline double getOffset(double phi) const;

		inline double smoothOrigin(double r) const;
		
		inline double getDistance(gravis::d2Vector imgPos) const;

		
		double scanForMinimalRadius(int samples) const;
		double scanForMaximalRadius(int samples) const;
		
		std::pair<gravis::d2Vector,gravis::d2Vector> scanForBoundingBox(
		        double radius, double padding, int samples) const;

		static std::vector<double> rotate(
				const std::vector<double>& params,
				double angle,
				gravis::d2Vector axis);
		
};

class DelineatedBlob2D : public Blob2D
{
	public: 
		
		DelineatedBlob2D();
		DelineatedBlob2D(const Blob2D& blob, double radius);
		DelineatedBlob2D(gravis::d2Vector center, double radius, double smoothingRadius = 0.0);
		DelineatedBlob2D(const std::vector<double>& params);
		
			double radius;
			
			
		double getRadius(double phi) const;
		double getSignedDistance(gravis::d2Vector imgPos) const;
		double getRelativeSignedDistance(gravis::d2Vector imgPos) const;

		gravis::d2Vector getOutlinePoint(double phi) const;
		gravis::d2Vector estimateNormal(double phi, double scale) const;
		
		double perimeter() const;
		Blob2D getBlob2D() const;
		
		static std::vector<DelineatedBlob2D> read(const std::string& filename);
		static std::vector<double> stripRadius(const std::vector<double>& params);
		static std::vector<double> addRadius(double radius, const std::vector<double>& params);
};


inline std::vector<double> Blob2D::toVector() const
{
	std::vector<double> out(2 * amplitudes.size() + 2);
	
	for (int i = 0; i < 2; i++)
	{
		out[i] = center[i];
	}
	
	for (int i = 0; i < amplitudes.size(); i++)
	{
		out[2*i+2] = amplitudes[i].real;
		out[2*i+3] = amplitudes[i].imag;
	}
	
	return out;
}

inline int Blob2D::findMaxRadius(gravis::i2Vector imageSize) const
{
	double maxRad = 0;
	
	for (int y = 0; y < imageSize.y; y++)
	for (int x = 0; x < imageSize.x; x++)
	{
		const double dx = x + center.x;
		const double dy = y + center.y;

		const double r = sqrt(dx*dx + dy*dy);
		
		if (r > maxRad) maxRad = r;
	}
		
	return (int) std::ceil(maxRad);
}

inline double Blob2D::getOffset(gravis::d2Vector v) const
{
	const int cc = amplitudes.size();
	
	if (cc < 1) return 0.0;
	
	const double phi = atan2(v.y, v.x);

	return getOffset(phi);
}

inline double Blob2D::getOffset(double phi) const
{
	const int cc = amplitudes.size();
	
	double out = 0.0;
	
	for (int i = 0; i < cc; i++)
	{
		out += amplitudes[i].real * cos((i+FIRST_BLOB_FREQUENCY) * phi);
		out += amplitudes[i].imag * sin((i+FIRST_BLOB_FREQUENCY) * phi);
	}
	
	return out;
}

inline double Blob2D::smoothOrigin(double r) const
{
	if (r < smoothingRadius)
	{
		/* f(x)  = x²/2R + R/2   =>  f(R) = R,  f(0) = R/2		   
		   f'(x) = x/R           =>  f'(R) = 1             */
		
		return 0.5 * r * r / smoothingRadius + smoothingRadius / 2;
	}
	else
	{
		return r;
	}
}

inline double Blob2D::getDistance(gravis::d2Vector imgPos) const
{
	const gravis::d2Vector d = imgPos - center;
	return d.length() + getOffset(d);
}

#endif

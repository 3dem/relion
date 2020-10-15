#ifndef POINT_INSERTION_H
#define POINT_INSERTION_H

#include <src/jaz/gravis/t3Vector.h>
#include <src/complex.h>
#include <src/jaz/image/raw_image.h>
#include <src/jaz/optics/dual_contrast/dual_contrast_voxel.h>


template <typename SrcType, typename DestType>
class ClippedPointInsertion
{
	public:
		
		inline void insert(
		        const tComplex<SrcType>& value, 
		        const SrcType& weight, 
		        const gravis::d3Vector& pos,
		        RawImage<tComplex<DestType>>& value_out,
		        RawImage<DestType>& weight_out) const;
		
		inline void insert_dualContrast(
		        const gravis::t2Vector<tComplex<SrcType>>& value, 
		        const gravis::t3Vector<SrcType>& weight, 
		        const gravis::d3Vector& pos,
		        RawImage<DualContrastVoxel<DestType>>& dest) const;
};

class WrappedPointInsertion
{};

template <typename SrcType, typename DestType>
inline void ClippedPointInsertion<SrcType, DestType>::insert(
        const tComplex<SrcType>& value, 
        const SrcType& weight, 
        const gravis::d3Vector& pos,
        RawImage<tComplex<DestType>>& value_out,
        RawImage<DestType>& weight_out) const
{
	const int wh3 = value_out.xdim;
	const int h3  = value_out.ydim;
	const int d3  = value_out.zdim;
	
	const int x0 = std::floor(pos.x);
	const int y0 = std::floor(pos.y);
	const int z0 = std::floor(pos.z);
	
	for (int dz = 0; dz < 2; dz++)
	for (int dy = 0; dy < 2; dy++)
	for (int dx = 0; dx < 2; dx++)
	{
		const int xg = x0 + dx;			
		const int yg = y0 + dy;
		const int zg = z0 + dz;
		
		if ( xg < wh3
		  && yg >= -h3/2 && yg < h3/2
		  && zg >= -d3/2 && zg < d3/2)
		{
			const int xi = xg;			
			const int yi = yg >= 0? yg : yg + h3;
			const int zi = zg >= 0? zg : zg + d3;
			
			const double fx = 1.0 - std::abs(pos.x - xg);
			const double fy = 1.0 - std::abs(pos.y - yg);
			const double fz = 1.0 - std::abs(pos.z - zg);
			
			const double m = fx * fy * fz;
			
			value_out( xi,yi,zi) += m * value;
			weight_out(xi,yi,zi) += m * weight;
			
			if (xi == 0)
			{
				const int yim = (h3 - yi) % h3;
				const int zim = (d3 - zi) % d3;
				
				value_out( 0,yim,zim) += m * value.conj();
				weight_out(0,yim,zim) += m * weight;
			}
		}
	}
}

template <typename SrcType, typename DestType>
inline void ClippedPointInsertion<SrcType, DestType>::insert_dualContrast(
        const gravis::t2Vector<tComplex<SrcType>>& value, 
        const gravis::t3Vector<SrcType>& weight, 
        const gravis::d3Vector& pos,
        RawImage<DualContrastVoxel<DestType>>& dest) const
{
	const int wh3 = dest.xdim;
	const int h3  = dest.ydim;
	const int d3  = dest.zdim;
	
	const int x0 = std::floor(pos.x);
	const int y0 = std::floor(pos.y);
	const int z0 = std::floor(pos.z);
	
	for (int dz = 0; dz < 2; dz++)
	for (int dy = 0; dy < 2; dy++)
	for (int dx = 0; dx < 2; dx++)
	{
		const int xg = x0 + dx;			
		const int yg = y0 + dy;
		const int zg = z0 + dz;
		
		if ( xg < wh3
		  && yg >= -h3/2 && yg < h3/2
		  && zg >= -d3/2 && zg < d3/2)
		{
			const int xi = xg;			
			const int yi = yg >= 0? yg : yg + h3;
			const int zi = zg >= 0? zg : zg + d3;
			
			const double fx = 1.0 - std::abs(pos.x - xg);
			const double fy = 1.0 - std::abs(pos.y - yg);
			const double fz = 1.0 - std::abs(pos.z - zg);
			
			const double m = fx * fy * fz;
			
			dest(xi,yi,zi).data_sin += m * value[0];
			dest(xi,yi,zi).data_cos += m * value[1];
			
			dest(xi,yi,zi).weight_sin2    += m * weight[0];
			dest(xi,yi,zi).weight_sin_cos += m * weight[1];
			dest(xi,yi,zi).weight_cos2    += m * weight[2];
			
			if (xi == 0)
			{
				const int yim = (h3 - yi) % h3;
				const int zim = (d3 - zi) % d3;
				
				dest(0,yim,zim).data_sin += m * value[0];
				dest(0,yim,zim).data_cos += m * value[1];
				
				dest(0,yim,zim).weight_sin2    += m * weight[0];
				dest(0,yim,zim).weight_sin_cos += m * weight[1];
				dest(0,yim,zim).weight_cos2    += m * weight[2];
			}
		}
	}
}

#endif

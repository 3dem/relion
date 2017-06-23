// WARNING - the following is dicey - CPU may use different values for some
// constants in some situations - hopefully none of the kernels here.
//#include "src/acc/cuda/cuda_settings.h"


namespace CpuKernels
{

template <typename T>
void multiply(
	T *A,
	T *OUT,
	T S,
	int image_size)
{
	for (int i = 0; i < image_size; i ++)
		OUT[i] = A[i]*S;
}

template <typename T>
void multiply(
	T *A,
	T S,
	int image_size)
{
	for (int i = 0; i < image_size; i ++)
		A[i] *= S;
}

template <typename T>
void multiply(
	T *A,
	T *B,
	T *OUT,
	T S,
	int image_size)
{
	for (int i = 0; i < image_size; i ++)
		OUT[i] = A[i]*B[i]*S;
}


template <typename T>
void translate2D(
	T * g_image_in,
	T * g_image_out,
	int image_size,
	int xdim,
	int ydim,
	int dx,
	int dy)
{
	int x,y,xp,yp;
	int new_pixel;

	for(int pixel=0; pixel<image_size; pixel++)
	{
		x = pixel % xdim;
		y = (pixel-x) / (xdim);

		xp = x + dx;
		yp = y + dy;

		if( yp>=0 && xp>=0 && yp<ydim && xp<xdim)
		{
			new_pixel = yp*xdim + xp;
			if(new_pixel>=0 && new_pixel<image_size) // if displacement is negative, new_pixel could be less than 0
				g_image_out[new_pixel] = g_image_in[pixel];
		}
	}
}

template <typename T>
void translate3D(
	T * g_image_in,
	T * g_image_out,
	int image_size,
	int xdim,
	int ydim,
	int zdim,
	int dx,
	int dy,
	int dz)
{
	int x,y,z,xp,yp,zp,xy;
	int new_voxel;

	for(int voxel=0; voxel<image_size; voxel++)
	{
		int xydim = xdim*ydim;

		if(voxel<image_size)
		{
			z =  voxel / xydim;
			zp = z + dz;

			xy = voxel % xydim;
			y =  xy / xdim;
			yp = y + dy;

			x =  xy % xdim;
			xp = x + dx;

			if( zp>=0 && yp>=0 && xp>=0 && zp<zdim && yp<ydim && xp<xdim)
			{
				new_voxel = zp*xydim +  yp*xdim + xp;
				if(new_voxel>=0 && new_voxel<image_size) // if displacement is negative, new_pixel could be less than 0
					g_image_out[new_voxel] = g_image_in[voxel];
			}
		}
	}
}

} //namespace CpuKernels


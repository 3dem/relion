#include <src/jaz/image/buffered_image.h>
#include <src/jaz/math/fft.h>
#include <src/jaz/tomography/projection/Fourier_backprojection.h>
#include <src/jaz/util/zio.h>

using namespace gravis;


int main(int argc, char** argv)
{
	const int box = 512;

	BufferedImage<float> box2D(box, box);

	box2D.fill(1.f);

	const double dot_rad = box/16;
	const double dot_edge = 3;

	for (int y = 0; y < box; y++)
	for (int x = 0; x < box; x++)
	{
		const double xx = x < box/2? x : x - box;
		const double yy = y < box/2? y : y - box;
		const double r = sqrt(xx*xx + yy*yy);

		if (r < dot_rad)
		{
			box2D(x,y) = 1.f;
		}
		else if (r < dot_rad + dot_edge)
		{
			box2D(x,y) = 1.0 - (r - dot_rad) / dot_edge;
		}
		else
		{
			box2D(x,y) = 0.f;
		}
	}

	box2D.write("DEBUG_box2D_initial.mrc");


	BufferedImage<fComplex> box2D_FS(box/2 + 1, box);
	FFT::FourierTransform(box2D, box2D_FS, FFT::Both);



	std::vector<int> crops {512, 384, 256, 192, 128, 64};

	for (int ic = 0; ic < crops.size(); ic++)
	{
		const int crop = crops[ic];

		BufferedImage<float> box3D(crop, crop, crop);

		BufferedImage<fComplex> box3D_FS(crop/2 + 1, crop, crop);

		box3D_FS.fill(fComplex(0.f, 0.f));

		BufferedImage<float> box3D_div(crop, crop, crop);

		BufferedImage<float> ctf(box/2 + 1, box, box);
		ctf.fill(1.f);

		{
			BufferedImage<float>
					ctfImgFS(crop/2 + 1, crop, crop),
					multiImageFS(crop/2 + 1, crop, crop);

			ctfImgFS.fill(0.f);
			multiImageFS.fill(0.f);

			d4Matrix P;
			P.loadIdentity();

			P *= crop / (double) box;

			FourierBackprojection::backprojectSlice_forward_with_multiplicity(
					box2D_FS,
					ctf,
					P,
					box3D_FS,
					ctfImgFS,
					multiImageFS);

			box3D_FS *= (float) sqrt(box / (double) crop);

			BufferedImage<fComplex> div3D_FS(crop/2 + 1, crop, crop);

			for (size_t i = 0; i < div3D_FS.getSize(); i++)
			{
				div3D_FS[i] = box3D_FS[i] / (ctfImgFS[i] + 0.001);
			}

			FFT::inverseFourierTransform(div3D_FS, box3D_div, FFT::Both);
		}

		box3D_div.write("DEBUG_box3D_div_box_"+ZIO::itoa(box)+"_crop_"+ZIO::itoa(crop)+".mrc");

		FFT::inverseFourierTransform(box3D_FS, box3D, FFT::Both);

		box3D.write("DEBUG_box3D_data_box_"+ZIO::itoa(box)+"_crop_"+ZIO::itoa(crop)+".mrc");
	}


	return 0;
}

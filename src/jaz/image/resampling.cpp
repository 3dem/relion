#include "resampling.h"

gravis::i2Vector Resampling::getFourierCroppedSize2D(int w, int h, double factor, bool keepSizeEven)
{
	int w_out = (int)(w / factor);
	int h_out = (int)(h / factor);

	if (keepSizeEven)
	{
		w_out -= w_out % 2;
		h_out -= h_out % 2;
	}
	
	return gravis::i2Vector(w_out, h_out);
}

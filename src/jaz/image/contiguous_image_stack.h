#ifndef JAZ_CONTIGUOUS_IMAGE_STACK_H
#define JAZ_CONTIGUOUS_IMAGE_STACK_H

#include <src/image.h>
#include <src/jaz/image/raw_image.h>

#include <vector>


template <typename T>
class ContiguousImageStack
{
	public:

		ContiguousImageStack()
			: particle_count(0), frame_count(0), width(0), height(0)
		{}

		ContiguousImageStack(const ContiguousImageStack&) = delete;
		ContiguousImageStack& operator = (const ContiguousImageStack&) = delete;

		void resize(long int particleCount, long int frameCount,
				long int width, long int height)
		{
			images.clear();

			particle_count = particleCount;
			frame_count = frameCount;
			this->width = width;
			this->height = height;

			if (particleCount <= 0 || frameCount <= 0 || width <= 0 || height <= 0)
			{
				storage.clear();
				particle_count = frame_count = this->width = this->height = 0;
				return;
			}

			const long int imageCount = particleCount * frameCount;
			storage.data.setDimensions(width, height, 1, imageCount);
			storage.data.coreAllocateReuse();

			images.resize(imageCount);

			const long int imageSize = width * height;
			for (long int i = 0; i < imageCount; i++)
			{
				images[i].data.setDimensions(width, height, 1, 1);
				images[i].data.data = storage.data.data + i * imageSize;
				images[i].data.destroyData = false;
			}
		}

		long int particleCount() const { return particle_count; }
		long int frameCount() const { return frame_count; }
		long int xdim() const { return width; }
		long int ydim() const { return height; }

		Image<T>& operator () (long int particle, long int frame)
		{
			return images[particle * frame_count + frame];
		}

		const Image<T>& operator () (long int particle, long int frame) const
		{
			return images[particle * frame_count + frame];
		}

		RawImage<T> getParticleRef(long int particle)
		{
			return RawImage<T>(
				width, height, frame_count,
				storage.data.data + particle * frame_count * width * height);
		}

		const RawImage<T> getParticleRef(long int particle) const
		{
			return RawImage<T>(
				width, height, frame_count,
				storage.data.data + particle * frame_count * width * height);
		}

		const T* data() const { return storage.data.data; }
		T* data() { return storage.data.data; }

	private:

		long int particle_count, frame_count, width, height;
		Image<T> storage;
		std::vector<Image<T>> images;
};

#endif

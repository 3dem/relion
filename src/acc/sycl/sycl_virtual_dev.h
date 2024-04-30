#ifndef _SYCL_VIRTUAL_DEVICE_RELION_H
#define _SYCL_VIRTUAL_DEVICE_RELION_H

#include <cstdint>
#include <string>
#include "src/acc/settings.h"

enum class syclMallocType {device, host};
enum class syclBackendType {openCL, levelZero, CUDA, HIP, host};
enum class syclDeviceType {gpu, cpu, fpga, host};
enum class syclDeviceID {cpuAVX2, cpuAVX512, intelATS, intelARC, intelPVC};
enum class syclQueueType {outOfOrder, inOrder, enableProfiling};

class virtualSYCL
{
public:
	virtual ~virtualSYCL() = default;

	virtual void* syclMalloc(const size_t bytes, const syclMallocType type, const char *name = "") = 0;

	virtual void syclFree(void *ptr) = 0;

	virtual void syclMemcpy(void *dest, const void *src, const size_t bytes) = 0;
	virtual void syclMemcpyAfterWaitAll(void *dest, const void *src, const size_t bytes) = 0;

	virtual void syclMemset(void *ptr, const int value, const size_t bytes) = 0;
	virtual void syclMemsetAfterWaitAll(void *ptr, const int value, const size_t bytes) = 0;

	virtual void syclFillInt32(void *ptr, const std::int32_t value, const size_t count) = 0;
	virtual void syclFillInt32AfterWaitAll(void *ptr, const std::int32_t value, const size_t count) = 0;
	virtual void syclFillUint32(void *ptr, const std::uint32_t value, const size_t count) = 0;
	virtual void syclFillUint32AfterWaitAll(void *ptr, const std::uint32_t value, const size_t count) = 0;

	virtual void syclFillInt64(void *ptr, const std::int64_t value, const size_t count) = 0;
	virtual void syclFillInt64AfterWaitAll(void *ptr, const std::int64_t value, const size_t count) = 0;
	virtual void syclFillUint64(void *ptr, const std::uint64_t value, const size_t count) = 0;
	virtual void syclFillUint64AfterWaitAll(void *ptr, const std::uint64_t value, const size_t count) = 0;

	virtual void syclFillXfloat(void *ptr, const XFLOAT value, const size_t count) = 0;
	virtual void syclFillXfloatAfterWaitAll(void *ptr, const XFLOAT value, const size_t count) = 0;
	virtual void syclFillRfloat(void *ptr, const RFLOAT value, const size_t count) = 0;
	virtual void syclFillRfloatAfterWaitAll(void *ptr, const RFLOAT value, const size_t count) = 0;

	virtual void syclFillFloat(void *ptr, const float value, const size_t count) = 0;
	virtual void syclFillFloatAfterWaitAll(void *ptr, const float value, const size_t count) = 0;
	virtual void syclFillDouble(void *ptr, const double value, const size_t count) = 0;
	virtual void syclFillDoubleAfterWaitAll(void *ptr, const double value, const size_t count) = 0;

	virtual void syclPrefetch(void *ptr, const size_t bytes) = 0;
	virtual void syclPrefetchAfterWaitAll(void *ptr, const size_t bytes) = 0;

	virtual bool isAsyncQueue() const = 0;
	virtual bool canSupportFP64() const = 0;
	virtual syclQueueType getSyclQueueType() const = 0;

	virtual void destroyMemoryPool() = 0;

	virtual std::string getName() = 0;
	virtual int getDeviceID() const = 0;
	virtual void setDeviceID(int id) = 0;
	virtual int getCardID() const = 0;
	virtual void setCardID(int id) = 0;
	virtual int getStackID() const = 0;
	virtual void setStackID(int id) = 0;
	virtual int getNumStack() const = 0;
	virtual void setNumStack(int id) = 0;
	virtual int getSliceID() const = 0;
	virtual void setSliceID(int id) = 0;
	virtual int getNumSlice() const = 0;
	virtual void setNumSlice(int id) = 0;

	virtual void waitAll() = 0;

	virtual void printDeviceInfo(bool printAll) = 0;
};

struct virtualSyclPtr
{
	virtualSYCL *dev;
};

#endif

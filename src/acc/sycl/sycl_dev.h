#ifndef _SYCL_DEVICE_RELION_H
#define _SYCL_DEVICE_RELION_H

#include "src/acc/sycl/sycl_virtual_dev.h"
#include "src/acc/sycl/sycl_memory_pool.h"

#include <cassert>
#include <string>
#include <array>
#include <mutex>
#include <vector>
#include <cstdint>
#include <sycl/sycl.hpp>
#include <sycl/backend.hpp>
#ifdef USE_ONEDPL
 #include <oneapi/dpl/algorithm>
 #include <oneapi/dpl/execution>
 #include <oneapi/dpl/iterator>
 #include <oneapi/dpl/numeric>
 #include <oneapi/dpl/ranges>
#endif

class devSYCL : public virtualSYCL
{
public:
	devSYCL(const devSYCL &q);
	devSYCL(const devSYCL *q);

	devSYCL(const bool isInOrder = false, const bool isAsync = false);
	devSYCL(sycl::device &d, int id, const bool isInOrder = false, const bool isAsync = false);
	devSYCL(sycl::device &d, const syclQueueType qType, int id, const bool isAsync = false);
	devSYCL(sycl::context &c, sycl::device &d, int id, const bool isInOrder = false, const bool isAsync = false);
	devSYCL(sycl::context &c, sycl::device &d, const syclQueueType qType, int id, const bool isAsync = false);
	devSYCL(const syclDeviceType dev, int id, const bool isInOrder = false, const bool isAsync = false);
	devSYCL(const syclBackendType be, const syclDeviceType dev, int id, const bool isInOrder = false, const bool isAsync = false);
#ifdef SYCL_EXT_INTEL_QUEUE_INDEX
	devSYCL(sycl::device &d, const int qIndex, int id, const bool isInOrder = false, const bool isAsync = false);
	devSYCL(sycl::device &d, const int qIndex, const syclQueueType qType, int id, const bool isAsync = false);
	devSYCL(sycl::context &c, sycl::device &d, const int qIndex, int id, const bool isInOrder = false, const bool isAsync = false);
	devSYCL(sycl::context &c, sycl::device &d, const int qIndex, const syclQueueType qType, int id, const bool isAsync = false);
#endif
	~devSYCL();

	void pushEvent(const sycl::event &evt)
	{
		if (_queueType == syclQueueType::inOrder)
			_prev_submission = evt;
		else
			_event.push_back(evt);
	}

	template <typename T>
	sycl::event syclSubmitAsync(T&& f)
	{
		_event.push_back(_devQ->submit(f));
		return getLastEvent();
	}

	template <typename T>
	sycl::event syclSubmitSync(T&& f)
	{
		_devQ->submit(f).wait_and_throw();
		return sycl::event();
	}

	template <typename T>
	sycl::event syclSubmitSeq(T&& f)
	{
		_prev_submission = _devQ->submit(f);
		return _prev_submission;
	}

	template <typename T>
	sycl::event syclSubmit(T&& f)
	{
		if (_queueType == syclQueueType::inOrder)
			return syclSubmitSeq(f);
		else
		{
			if (_isAsyncQueue)
				return syclSubmitAsync(f);
			else
				return syclSubmitSync(f);
		}
	}

	template <typename T>
	T* syclMalloc(const size_t count, const syclMallocType type = syclMallocType::device, const char *name = "")
	{ return static_cast<T*>(syclMalloc(count * sizeof(T), type, name)); }

	void* syclMalloc(const size_t bytes, const syclMallocType type, const char *name = "") override;

	void syclFree(void *ptr) override;

	void syclMemcpy(void *dest, const void *src, const size_t bytes) override;
	void syclMemcpyAfterWaitAll(void *dest, const void *src, const size_t bytes) override;

	void syclMemset(void *ptr, const int value, const size_t bytes) override;
	void syclMemsetAfterWaitAll(void *ptr, const int value, const size_t bytes) override;

	template <typename T>
	void syclFill(void *ptr, const T &pattern, const size_t count)
	{
		assert( ptr != nullptr);
		assert(count > 0);

		pushEvent(_devQ->fill<T>(ptr, pattern, count));
	}

	template <typename T>
	void syclFillAfterWaitAll(void *ptr, const T &pattern, const size_t count)
	{
		assert( ptr != nullptr);
		assert(count > 0);

		waitAll();
		pushEvent(_devQ->fill<T>(ptr, pattern, count));
	}

	void syclFillInt32(void *ptr, const std::int32_t value, const size_t count) override;
	void syclFillInt32AfterWaitAll(void *ptr, const std::int32_t value, const size_t count) override;
	void syclFillUint32(void *ptr, const std::uint32_t value, const size_t count) override;
	void syclFillUint32AfterWaitAll(void *ptr, const std::uint32_t value, const size_t count) override;

	void syclFillInt64(void *ptr, const std::int64_t value, const size_t count) override;
	void syclFillInt64AfterWaitAll(void *ptr, const std::int64_t value, const size_t count) override;
	void syclFillUint64(void *ptr, const std::uint64_t value, const size_t count) override;
	void syclFillUint64AfterWaitAll(void *ptr, const std::uint64_t value, const size_t count) override;

	void syclFillXfloat(void *ptr, const XFLOAT value, const size_t count) override;
	void syclFillXfloatAfterWaitAll(void *ptr, const XFLOAT value, const size_t count) override;
	void syclFillRfloat(void *ptr, const RFLOAT value, const size_t count) override;
	void syclFillRfloatAfterWaitAll(void *ptr, const RFLOAT value, const size_t count) override;

	void syclFillFloat(void *ptr, const float value, const size_t count) override;
	void syclFillFloatAfterWaitAll(void *ptr, const float value, const size_t count) override;
	void syclFillDouble(void *ptr, const double value, const size_t count) override;
	void syclFillDoubleAfterWaitAll(void *ptr, const double value, const size_t count) override;

	void syclPrefetch(void *ptr, const size_t bytes) override;
	void syclPrefetchAfterWaitAll(void *ptr, const size_t bytes) override;

	void printDeviceInfo(bool printAll = false) override;

	void waitAll() override;

	bool isAsyncQueue() const override		{ return _isAsyncQueue; }
	bool canSupportFP64() const override	{ return _isFP64Supported; }
	syclQueueType getSyclQueueType() const override	{ return _queueType; }
#ifdef SYCL_EXT_INTEL_QUEUE_INDEX
	int getComputeIndex() const		{ return _computeIndex; }
#else
	int getComputeIndex() const		{ return -1; }
#endif

	std::string getName() override;
	int getDeviceID() const override	{ return deviceID; }
	void setDeviceID(int id) override	{ deviceID = id; }
	int getCardID() const override	{ return cardID; }
	void setCardID(int id) override	{ cardID = id; }
	int getStackID() const override	{ return stackID; }
	void setStackID(int id) override	{ stackID = id; }
	int getNumStack() const override	{ return nStack; }
	void setNumStack(int n) override	{ nStack = n; }
	int getSliceID() const override	{ return sliceID; }
	void setSliceID(int id) override	{ sliceID = id; }
	int getNumSlice() const override	{ return nSlice; }
	void setNumSlice(int n) override	{ nSlice = n; }

	void reCalculateRange(sycl::range<3> &wg, const sycl::range<3> &wi);

	void setSyclHostPool(syclMemoryPool *sM, const size_t size)		{ _syclHostPool = sM; _host_block_size = size; }
	void setSyclDevicePool(syclMemoryPool *sM, const size_t size)	{ _syclDevicePool = sM; _device_block_size = size; }
	syclMemoryPool* getSyclHostPool() const		{ return _syclHostPool; }
	syclMemoryPool* getSyclDevicePool() const	{ return _syclDevicePool; }
	const size_t getSyclHostBlockSize() const	{ return _host_block_size; }
	const size_t getSyclDeviceBlockSize() const	{ return _device_block_size; }
	void destroyMemoryPool() override;

	sycl::queue* getQueue()		{ return _devQ; }
	sycl::device& getDevice()		{ return _devD; }
	sycl::context& getContext()		{ return _devC; }
	sycl::event& getLastEvent()			{ return _event.back(); }
	sycl::event& getPrevSubmission()	{ return _prev_submission; }
	std::vector<sycl::event>& getEventList()	{ return _event; }
#ifdef USE_ONEDPL
	oneapi::dpl::execution::device_policy<>& getDevicePolicy()	{ return _devicePolicy; }
#endif

	std::array<int, 3> maxItem;	// Maximum work item size in each 3-dimension
	size_t maxGroup;	// Maximum work item size per group
	uint64_t localMem;	// Maximum local memory size per group
	uint64_t globalMem;	// Maximum device memory size per queue
	uint32_t maxUnit;	// Maximum compute units
	std::array<long, 3> maxWorkGroup;	// Maximum number of work-groups that can be submitted in each dimension
	size_t maxGlobalWorkGroup;			// Maximum number of work-groups that can be submitted across all the dimensions

private:
	syclQueueType _queueType;
	sycl::queue *_devQ;
	sycl::device _devD;
	sycl::context _devC;
	sycl::async_handler _exHandler;
	std::vector<sycl::event> _event;
	sycl::event _prev_submission;
	std::string deviceName;
	syclMemoryPool *_syclHostPool;
	syclMemoryPool *_syclDevicePool;
	size_t _host_block_size;
	size_t _device_block_size;
	int cardID;
	int deviceID;
	int stackID;
	int nStack;
	int sliceID;
	int nSlice;
	bool _isAsyncQueue;
	bool _isFP64Supported;
	int _computeIndex;
#ifdef USE_ONEDPL
	oneapi::dpl::execution::device_policy<> _devicePolicy;
#endif
	inline static std::mutex _lock;
};

#endif

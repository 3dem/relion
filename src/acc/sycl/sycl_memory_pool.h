#pragma once

#include <mutex>
#include <atomic>
#include <cstddef>
#include <utility>
#include <list>
#include <map>
#include <sycl/sycl.hpp>

using memoryStorage = std::byte*;
using memoryAllocMap = std::map<memoryStorage,memoryStorage>;
using memoryFreeMap = std::map<memoryStorage,memoryStorage>;
using memoryFreeSizeMap = std::multimap<size_t,std::pair<memoryStorage,memoryStorage>>;
using alloc_kind = sycl::usm::alloc;

constexpr size_t memoryAlignSize = 128;	// 128 Byte (1024-bit)
constexpr size_t memoryPadSize = 128;	// 128 Byte (1024-bit)
constexpr size_t defDeviceBlockSize = 256UL*1024*1024;	// 256 MiB
constexpr size_t defHostBlockSize = 256UL*1024*1024;	// 256 MiB
constexpr size_t defSharedBlockSize = 256UL*1024*1024;	// 256 MiB
constexpr size_t numMutex = 256;

struct alignas(memoryAlignSize) alignedMutex
{
	std::mutex mutex;
};

class syclMemoryBlock
{
	public:
		syclMemoryBlock(sycl::queue* const Q, const alloc_kind type, const size_t block_bytes);
		syclMemoryBlock(const syclMemoryBlock&) = delete;
		syclMemoryBlock() = delete;
		~syclMemoryBlock();

		void* alloc(const size_t bytes);
		bool free(void* const pos);
		void clean();
		const size_t size() const		{ return _block_bytes; }
		const size_t alloc_size() const	{ return _alloc_bytes; }
		const size_t free_size() const	{ return _free_bytes; }
		const bool empty() const;
		const bool isUnused() const;
		void printStats() const;

	private:
		size_t _blockID;
		sycl::queue* const _Q;
		const alloc_kind _mem_type;
		const size_t _block_bytes;
		memoryStorage _mem;
		memoryAllocMap _alloc_list;
		memoryFreeMap _free_list;
		memoryFreeMap _rfree_list;		// Holds value-key pair of _free_list for faster backward free-list lookup
		size_t _free_bytes;
		size_t _alloc_bytes;
		size_t _max_alloc;
		inline static std::atomic<size_t> _count {0};
		inline static alignedMutex _container_lock[numMutex];

	friend class syclMemoryPool;
};

class syclMemoryPool
{
	public:
		syclMemoryPool(sycl::queue* const Q, const alloc_kind type, const size_t block_size = 0);
		syclMemoryPool(const syclMemoryPool&) = delete;
		syclMemoryPool() = delete;
		~syclMemoryPool();

		template <typename T>
		T* alloc(const size_t count)	{ return static_cast<T*>(alloc(count * sizeof(T))); }
		void* alloc(const size_t bytes);
		bool free(void* const pos);
		void clean();
		void shrink(const size_t num = 1);
		const size_t getAllocSize(const size_t size) const;
		const size_t get_alloc_size() const;
		const size_t get_free_size() const;
		const size_t get_pool_size() const;
		const alloc_kind get_memory_type() const	{ return _mem_type; }
		const bool is_host() const		{ return _mem_type == alloc_kind::host; }
		const bool is_device() const	{ return _mem_type == alloc_kind::device; }
		const bool is_shared() const	{ return _mem_type == alloc_kind::shared; }
		void printStats() const;

	private:
		std::list<syclMemoryBlock>::iterator addNewBlock();

		sycl::queue* const _Q;
		size_t _block_size;
		const alloc_kind _mem_type;
		std::list<syclMemoryBlock> _memory_block_list;
		inline static std::mutex _mem_lock;
};

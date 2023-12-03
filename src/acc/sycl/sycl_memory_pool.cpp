#include "src/acc/sycl/sycl_memory_pool.h"

#include <iostream>
#include <cstring>


syclMemoryBlock::syclMemoryBlock(sycl::queue* const Q, const alloc_kind type, const size_t block_bytes) : _Q {Q}, _mem_type {type}, _block_bytes {block_bytes}
{
	_mem = nullptr;
	if (_mem_type == alloc_kind::device)
		_mem = static_cast<memoryStorage>(sycl::aligned_alloc_device(memoryAlignSize, _block_bytes, *_Q));
	else if (_mem_type == alloc_kind::host)
		_mem = static_cast<memoryStorage>(sycl::aligned_alloc_host(memoryAlignSize, _block_bytes, *_Q));
	else if (_mem_type == alloc_kind::shared)
		_mem = static_cast<memoryStorage>(sycl::aligned_alloc_shared(memoryAlignSize, _block_bytes, *_Q));

	if (_mem == nullptr)
	{
		_blockID = numMutex-1;
		_alloc_bytes = 0;
		_max_alloc = 0;
		_free_bytes = 0;
	}
	else
	{
		_blockID = _count % numMutex;
		_count++;

		_alloc_bytes = 0;
		_max_alloc = 0;
		_free_bytes = _block_bytes;
		_free_list.emplace(_mem, _mem + _block_bytes);
		_rfree_list.emplace(_mem + _block_bytes, _mem);
	}
}

syclMemoryBlock::~syclMemoryBlock()
{
	if (_mem)
	{
		sycl::free(_mem, *_Q);
#ifdef DIAG_SYCL_MEMORY_POOL
		printStats();
#endif
	}
}

void syclMemoryBlock::clean()
{
	if (_mem)
	{
		std::lock_guard<std::mutex> locker(_container_lock[_blockID].mutex);
#ifdef DIAG_SYCL_MEMORY_POOL
		printStats();
#endif
		_alloc_bytes = 0;
		_max_alloc = 0;
		_free_bytes = _block_bytes;

		_alloc_list.clear();
		_free_list.clear();
		_free_list.emplace(_mem, _mem + _block_bytes);
		_rfree_list.clear();
		_rfree_list.emplace(_mem + _block_bytes, _mem);
	}
}

const bool syclMemoryBlock::isUnused() const
{
	if (_mem != nullptr && _block_bytes > 0 && _alloc_list.size() == 0)
		return true;
	else
		return false;
}

const bool syclMemoryBlock::empty() const
{
	if (_mem == nullptr)
		return true;
	else
		return false;
}

void syclMemoryBlock::printStats() const
{
	if (_mem)
	{
		std::cout << "[" << _blockID << "] " << (_mem_type == alloc_kind::device ? "DEVICE" : "HOST") << " :Efficiency= " << 100*_alloc_bytes/_block_bytes << "% / " << 100*_max_alloc/_block_bytes << "%\n";
		std::cout << "[" << _blockID << "] " << (_mem_type == alloc_kind::device ? "DEVICE" : "HOST") << " :Alloc(" << _alloc_list.size() << ")/Free(" << _free_list.size() << ")/Total= " << _alloc_bytes << " " << _free_bytes << " " << _block_bytes << std::endl;
	}
}

void* syclMemoryBlock::alloc(const size_t bytes)
{
	std::lock_guard<std::mutex> locker(_container_lock[_blockID].mutex);

	if (bytes > _free_bytes)
		return nullptr;

	void* mem = nullptr;
	auto free = _free_list.begin();
	while (free != _free_list.cend())
	{
		if (free->second - free->first >= bytes)
		{
			mem = static_cast<void*>(free->first);
			_alloc_list.emplace(free->first, free->first + bytes);
			if (free->second > free->first + bytes)
			{
				_rfree_list[free->second] = free->first + bytes;

				auto node = _free_list.extract(free);
				node.key() = free->first + bytes;
				_free_list.insert(std::move(node));
			}
			else
			{
				_rfree_list.erase(free->second);
				_free_list.erase(free);
			}

			_alloc_bytes += bytes;
			_free_bytes -= bytes;
			if (_alloc_bytes > _max_alloc)
				_max_alloc = _alloc_bytes;
			break;
		}
		free++;
	}
#ifdef INIT_SYCL_MEMORY_TO_ZERO
	// If needed to set zero for allocated memory
	if (mem)
	{
		if (_mem_type == alloc_kind::device || _mem_type == alloc_kind::shared)
			_Q->memset(mem, 0, bytes).wait();
		else
			memset(mem, 0, bytes);
	}
#endif
	return mem;
}

bool syclMemoryBlock::free(void* const pos)
{
	std::lock_guard<std::mutex> locker(_container_lock[_blockID].mutex);

	bool isFree = false;
	const auto alloc_found = _alloc_list.find(static_cast<memoryStorage>(pos));
	if (alloc_found != _alloc_list.cend())
	{
		// Merge with forward free-list
		bool isForwardMerged = false;
		const auto free_found = _free_list.find(alloc_found->second);
		memoryStorage end_of_forward_free;
		if (free_found != _free_list.cend())
		{
			_rfree_list[free_found->second] = alloc_found->first;

			end_of_forward_free = free_found->second;
			auto node = _free_list.extract(free_found);
			node.key() = alloc_found->first;
			_free_list.insert(std::move(node));

			isForwardMerged = true;
		}
		// Merge with backward free-list
		bool isBackwardMerged = false;
		const auto rfree_found = _rfree_list.find(alloc_found->first);
		if (rfree_found != _rfree_list.cend())
		{
			if (isForwardMerged)
			{
				_free_list[rfree_found->second] = end_of_forward_free;
				_free_list.erase(alloc_found->first);

				_rfree_list[end_of_forward_free] = rfree_found->second;
				_rfree_list.erase(rfree_found);
			}
			else
			{
				_free_list[rfree_found->second] = alloc_found->second;

				auto node = _rfree_list.extract(rfree_found);
				node.key() = alloc_found->second;
				_rfree_list.insert(std::move(node));
			}
			isBackwardMerged = true;
		}

		if (isForwardMerged == false && isBackwardMerged == false)
		{
			_free_list.emplace(alloc_found->first, alloc_found->second);
			_rfree_list.emplace(alloc_found->second, alloc_found->first);
		}

		_alloc_bytes -= alloc_found->second - alloc_found->first;
		_free_bytes += alloc_found->second - alloc_found->first;
		_alloc_list.erase(alloc_found);
		isFree = true;
	}
	return isFree;
}

///////////////////////////////////////////////////////////////////////////////
syclMemoryPool::syclMemoryPool(sycl::queue* const Q, const alloc_kind type, const size_t block_size) : _Q {Q}, _mem_type {type}, _block_size {block_size}
{
	if (_block_size == 0)
	{
		if (_mem_type == alloc_kind::device)
			_block_size = defDeviceBlockSize;
		else if (_mem_type == alloc_kind::host)
			_block_size = defHostBlockSize;
		else if (_mem_type == alloc_kind::shared)
			_block_size = defSharedBlockSize;
	}
}

syclMemoryPool::~syclMemoryPool()
{
}

void syclMemoryPool::clean()
{
	for (auto& block : _memory_block_list)
		block.clean();
}

void syclMemoryPool::shrink(const size_t num)
{
	std::lock_guard<std::mutex> locker(_mem_lock);

	size_t count = 0;
	auto it = _memory_block_list.begin();
	while (it != _memory_block_list.end())
	{
		if (it->empty())
			it = _memory_block_list.erase(it);
		else if (count >= num && it->isUnused())
			it = _memory_block_list.erase(it);
		else
		{
			if (count < num && it->isUnused())
				it->clean();
			count++;
			it++;
		}
	}
}

const size_t syclMemoryPool::get_alloc_size() const
{
	size_t bytes;
	for (const auto& m : _memory_block_list)
		bytes += m._alloc_bytes;
	return bytes;
}

const size_t syclMemoryPool::get_free_size() const
{
	size_t bytes;
	for (const auto& m : _memory_block_list)
		bytes += m._free_bytes;
	return bytes;
}

const size_t syclMemoryPool::get_pool_size() const
{
	size_t bytes;
	for (const auto& m : _memory_block_list)
		bytes += m._block_bytes;
	return bytes;
}

void syclMemoryPool::printStats() const
{
	for (const auto& m : _memory_block_list)
		m.printStats();
}

const size_t syclMemoryPool::getAllocSize(const size_t size) const
{
	return (size/memoryAlignSize + (size%memoryAlignSize == 0 ? 0 : 1)) * memoryAlignSize + (memoryPadSize > memoryAlignSize ? memoryPadSize : memoryAlignSize);
}

std::list<syclMemoryBlock>::iterator syclMemoryPool::addNewBlock()
{
	std::lock_guard<std::mutex> locker(_mem_lock);

	_memory_block_list.emplace_back(_Q, _mem_type, (_block_size/memoryAlignSize)*memoryAlignSize);
	return --_memory_block_list.end();
}

void* syclMemoryPool::alloc(const size_t bytes)
{
	void *mem = nullptr;
	const size_t allocBytes = getAllocSize(bytes);
	for (auto& block : _memory_block_list)
	{
		mem = block.alloc(allocBytes);
		if (mem != nullptr)
			break;
	}
	if (mem == nullptr)
	{
		auto add = addNewBlock();
		if (add->empty())
			_memory_block_list.pop_back();
		else
			mem = add->alloc(allocBytes);
	}
	return mem;
}

bool syclMemoryPool::free(void* const pos)
{
	for (auto& block : _memory_block_list)
	{
		if (block.free(pos))
			return true;
	}
	return false;
}

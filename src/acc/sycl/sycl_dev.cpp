#include <string>
#include <sstream>
#include <iostream>
#include <exception>
#include <chrono>
#include <cassert>
#include <exception>
#include <limits>

#include <algorithm>
#include "src/acc/sycl/sycl_dev.h"

#define PRINT_DEV_INFO(sd, x) std::cout << "  "#x": " << sd.get_info<sycl::info::device::x>() << std::endl
#define PRINT_DEV_INFO2020(sd, x) std::cout << "  "#x": " << sd.has(sycl::aspect::x) << std::endl

sycl::async_handler exceptionHandler = [](sycl::exception_list exceptions)
{
	for (const std::exception_ptr &e : exceptions)
	{
		try
		{
			std::rethrow_exception(e);
		}
		catch (const sycl::exception &e)
		{
			std::cerr << "Caught asynchronous SYCL exception:\n" << e.what() << "\n" << std::flush;
		}
	}
};

devSYCL::devSYCL(const bool isInOrder, const bool isAsync)
{
	_exHandler = exceptionHandler;

	try
	{
		if (isInOrder)
		{
			_devQ = new sycl::queue(sycl::default_selector_v, _exHandler, sycl::property::queue::in_order{});
			_queueType = syclQueueType::inOrder;
		}
		else
		{
			_devQ = new sycl::queue(sycl::default_selector_v, _exHandler);
			_queueType = syclQueueType::outOfOrder;
		}
	}
	catch (const sycl::exception &e)
	{
		std::cerr << "There is no available SYCL device.\n" << e.what() << std::endl;
		std::terminate();
	}

	_devD = _devQ->get_device();
	_devC = _devQ->get_context();
	_syclHostPool = nullptr;
	_syclDevicePool = nullptr;
	_host_block_size = 0;
	_device_block_size = 0;
#ifdef USE_ONEDPL
	_devicePolicy = oneapi::dpl::execution::make_device_policy(*_devQ);
#endif
	auto d = _devD;
	auto pf = d.get_platform();
	deviceName = d.get_info<sycl::info::device::name>() + " (" + d.get_info<sycl::info::device::driver_version>() + ") / " + pf.get_info<sycl::info::platform::name>();
	_isFP64Supported = d.has(sycl::aspect::fp64);
	_isAsyncQueue = isAsync;
	const auto isizes = d.get_info<sycl::info::device::max_work_item_sizes<3>>();
	maxItem[0] = isizes[0];
	maxItem[1] = isizes[1];
	maxItem[2] = isizes[2];
	maxGroup = d.get_info<sycl::info::device::max_work_group_size>();
#ifdef SYCL_EXT_ONEAPI_MAX_WORK_GROUP_QUERY
	auto groups = d.get_info<sycl::ext::oneapi::experimental::info::device::max_work_groups<3>>();
	maxWorkGroup[0] = groups[0];
	maxWorkGroup[1] = groups[1];
	maxWorkGroup[2] = groups[2];
	maxGlobalWorkGroup = d.get_info<sycl::ext::oneapi::experimental::info::device::max_global_work_groups>();
#else
	maxWorkGroup[0] = maxWorkGroup[1] = maxWorkGroup[2] = maxGlobalWorkGroup = std::numeric_limits<int>::max;
#endif
	globalMem = d.get_info<sycl::info::device::global_mem_size>();
	localMem = d.get_info<sycl::info::device::local_mem_size>();
	maxUnit = d.get_info<sycl::info::device::max_compute_units>();
	cardID = -1;
	deviceID = -1;
	stackID = -1;
	nStack = -1;
	sliceID = -1;
	nSlice = -1;
	_computeIndex = -1;
	_prev_submission = sycl::event();
	_event.push_back(sycl::event());
}

devSYCL::devSYCL(sycl::device &d, int id, const bool isInOrder, const bool isAsync)
{
	_exHandler = exceptionHandler;

	try
	{
		if (isInOrder)
		{
			_devQ = new sycl::queue(d, _exHandler, sycl::property::queue::in_order{});
			_queueType = syclQueueType::inOrder;
		}
		else
		{
			_devQ = new sycl::queue(d, _exHandler);
			_queueType = syclQueueType::outOfOrder;
		}
	}
	catch (const sycl::exception &e)
	{
		std::cerr << "Provided SYCL device failed\n" << e.what() << std::endl;
		std::terminate();
	}

	_devD = d;
	_devC = _devQ->get_context();
	_syclHostPool = nullptr;
	_syclDevicePool = nullptr;
	_host_block_size = 0;
	_device_block_size = 0;
#ifdef USE_ONEDPL
	_devicePolicy = oneapi::dpl::execution::make_device_policy(*_devQ);
#endif
	_isFP64Supported = d.has(sycl::aspect::fp64);
	_isAsyncQueue = isAsync;

	const auto isizes = d.get_info<sycl::info::device::max_work_item_sizes<3>>();
	maxItem[0] = isizes[0];
	maxItem[1] = isizes[1];
	maxItem[2] = isizes[2];
	maxGroup = d.get_info<sycl::info::device::max_work_group_size>();
#ifdef SYCL_EXT_ONEAPI_MAX_WORK_GROUP_QUERY
	auto groups = d.get_info<sycl::ext::oneapi::experimental::info::device::max_work_groups<3>>();
	maxWorkGroup[0] = groups[0];
	maxWorkGroup[1] = groups[1];
	maxWorkGroup[2] = groups[2];
	maxGlobalWorkGroup = d.get_info<sycl::ext::oneapi::experimental::info::device::max_global_work_groups>();
#else
	maxWorkGroup[0] = maxWorkGroup[1] = maxWorkGroup[2] = maxGlobalWorkGroup = std::numeric_limits<int>::max;
#endif
	globalMem = d.get_info<sycl::info::device::global_mem_size>();
	localMem = d.get_info<sycl::info::device::local_mem_size>();
	maxUnit = d.get_info<sycl::info::device::max_compute_units>();
	deviceID = id;
	_computeIndex = -1;
	auto pf = d.get_platform();
	deviceName = d.get_info<sycl::info::device::name>() + " (" + d.get_info<sycl::info::device::driver_version>() + ") / " + pf.get_info<sycl::info::platform::name>();
	_prev_submission = sycl::event();
	_event.push_back(sycl::event());
}

devSYCL::devSYCL(sycl::device &d, const syclQueueType qType, int id, const bool isAsync)
{
	_exHandler = exceptionHandler;

	try
	{
		switch (qType)
		{
			case syclQueueType::inOrder :
				_devQ = new sycl::queue(d, _exHandler, sycl::property::queue::in_order{});
				break;

			case syclQueueType::enableProfiling :
				_devQ = new sycl::queue(d, _exHandler, sycl::property::queue::enable_profiling{});
				break;

			case syclQueueType::outOfOrder :
			default :
				_devQ = new sycl::queue(d, _exHandler); 
				break;
		}
	}
	catch (const sycl::exception &e)
	{
		std::cerr << "Provided SYCL device failed\n" << e.what() << std::endl;
		std::terminate();
	}

	_devD = d;
	_devC = _devQ->get_context();
	_syclHostPool = nullptr;
	_syclDevicePool = nullptr;
	_host_block_size = 0;
	_device_block_size = 0;
	_queueType = qType;
#ifdef USE_ONEDPL
	_devicePolicy = oneapi::dpl::execution::make_device_policy(*_devQ);
#endif
	_isFP64Supported = d.has(sycl::aspect::fp64);
	_isAsyncQueue = isAsync;

	const auto isizes = d.get_info<sycl::info::device::max_work_item_sizes<3>>();
	maxItem[0] = isizes[0];
	maxItem[1] = isizes[1];
	maxItem[2] = isizes[2];
	maxGroup = d.get_info<sycl::info::device::max_work_group_size>();
#ifdef SYCL_EXT_ONEAPI_MAX_WORK_GROUP_QUERY
	auto groups = d.get_info<sycl::ext::oneapi::experimental::info::device::max_work_groups<3>>();
	maxWorkGroup[0] = groups[0];
	maxWorkGroup[1] = groups[1];
	maxWorkGroup[2] = groups[2];
	maxGlobalWorkGroup = d.get_info<sycl::ext::oneapi::experimental::info::device::max_global_work_groups>();
#else
	maxWorkGroup[0] = maxWorkGroup[1] = maxWorkGroup[2] = maxGlobalWorkGroup = std::numeric_limits<int>::max;
#endif
	globalMem = d.get_info<sycl::info::device::global_mem_size>();
	localMem = d.get_info<sycl::info::device::local_mem_size>();
	maxUnit = d.get_info<sycl::info::device::max_compute_units>();
	deviceID = id;
	_computeIndex = -1;
	auto pf = d.get_platform();
	deviceName = d.get_info<sycl::info::device::name>() + " (" + d.get_info<sycl::info::device::driver_version>() + ") / " + pf.get_info<sycl::info::platform::name>();
	_prev_submission = sycl::event();
	_event.push_back(sycl::event());
}

devSYCL::devSYCL(sycl::context &c, sycl::device &d, int id, const bool isInOrder, const bool isAsync)
{
	_exHandler = exceptionHandler;

	try
	{
		if (isInOrder)
		{
			_devQ = new sycl::queue(c, d, _exHandler, sycl::property::queue::in_order{});
			_queueType = syclQueueType::inOrder;
		}
		else
		{
			_devQ = new sycl::queue(c, d, _exHandler);
			_queueType = syclQueueType::outOfOrder;
		}
	}
	catch (const sycl::exception &e)
	{
		std::cerr << "Provided SYCL device failed\n" << e.what() << std::endl;
		std::terminate();
	}
	_devD = d;
	_devC = c;
	_syclHostPool = nullptr;
	_syclDevicePool = nullptr;
	_host_block_size = 0;
	_device_block_size = 0;
#ifdef USE_ONEDPL
	_devicePolicy = oneapi::dpl::execution::make_device_policy(*_devQ);
#endif
	_isFP64Supported = d.has(sycl::aspect::fp64);
	_isAsyncQueue = isAsync;

	const auto isizes = d.get_info<sycl::info::device::max_work_item_sizes<3>>();
	maxItem[0] = isizes[0];
	maxItem[1] = isizes[1];
	maxItem[2] = isizes[2];
	maxGroup = d.get_info<sycl::info::device::max_work_group_size>();
#ifdef SYCL_EXT_ONEAPI_MAX_WORK_GROUP_QUERY
	auto groups = d.get_info<sycl::ext::oneapi::experimental::info::device::max_work_groups<3>>();
	maxWorkGroup[0] = groups[0];
	maxWorkGroup[1] = groups[1];
	maxWorkGroup[2] = groups[2];
	maxGlobalWorkGroup = d.get_info<sycl::ext::oneapi::experimental::info::device::max_global_work_groups>();
#else
	maxWorkGroup[0] = maxWorkGroup[1] = maxWorkGroup[2] = maxGlobalWorkGroup = std::numeric_limits<int>::max;
#endif
	globalMem = d.get_info<sycl::info::device::global_mem_size>();
	localMem = d.get_info<sycl::info::device::local_mem_size>();
	maxUnit = d.get_info<sycl::info::device::max_compute_units>();
	deviceID = id;
	_computeIndex = -1;
	auto pf = d.get_platform();
	deviceName = d.get_info<sycl::info::device::name>() + " (" + d.get_info<sycl::info::device::driver_version>() + ") / " + pf.get_info<sycl::info::platform::name>();
	_prev_submission = sycl::event();
	_event.push_back(sycl::event());
}

devSYCL::devSYCL(sycl::context &c, sycl::device &d, const syclQueueType qType, int id, const bool isAsync)
{
	_exHandler = exceptionHandler;

	try
	{
		switch (qType)
		{
			case syclQueueType::inOrder :
				_devQ = new sycl::queue(c, d, _exHandler, sycl::property::queue::in_order{});
				break;

			case syclQueueType::enableProfiling :
				_devQ = new sycl::queue(c, d, _exHandler, sycl::property::queue::enable_profiling{});
				break;

			case syclQueueType::outOfOrder :
			default :
				_devQ = new sycl::queue(c, d, _exHandler);
				break;
		}
	}
	catch (const sycl::exception &e)
	{
		std::cerr << "Provided SYCL device failed\n" << e.what() << std::endl;
		std::terminate();
	}
	_devD = d;
	_devC = c;
	_syclHostPool = nullptr;
	_syclDevicePool = nullptr;
	_host_block_size = 0;
	_device_block_size = 0;
	_queueType = qType;
#ifdef USE_ONEDPL
	_devicePolicy = oneapi::dpl::execution::make_device_policy(*_devQ);
#endif
	_isFP64Supported = d.has(sycl::aspect::fp64);
	_isAsyncQueue = isAsync;

	const auto isizes = d.get_info<sycl::info::device::max_work_item_sizes<3>>();
	maxItem[0] = isizes[0];
	maxItem[1] = isizes[1];
	maxItem[2] = isizes[2];
	maxGroup = d.get_info<sycl::info::device::max_work_group_size>();
#ifdef SYCL_EXT_ONEAPI_MAX_WORK_GROUP_QUERY
	auto groups = d.get_info<sycl::ext::oneapi::experimental::info::device::max_work_groups<3>>();
	maxWorkGroup[0] = groups[0];
	maxWorkGroup[1] = groups[1];
	maxWorkGroup[2] = groups[2];
	maxGlobalWorkGroup = d.get_info<sycl::ext::oneapi::experimental::info::device::max_global_work_groups>();
#else
	maxWorkGroup[0] = maxWorkGroup[1] = maxWorkGroup[2] = maxGlobalWorkGroup = std::numeric_limits<int>::max;
#endif
	globalMem = d.get_info<sycl::info::device::global_mem_size>();
	localMem = d.get_info<sycl::info::device::local_mem_size>();
	maxUnit = d.get_info<sycl::info::device::max_compute_units>();
	deviceID = id;
	_computeIndex = -1;
	auto pf = d.get_platform();
	deviceName = d.get_info<sycl::info::device::name>() + " (" + d.get_info<sycl::info::device::driver_version>() + ") / " + pf.get_info<sycl::info::platform::name>();
	_prev_submission = sycl::event();
	_event.push_back(sycl::event());
}

#ifdef SYCL_EXT_INTEL_QUEUE_INDEX
devSYCL::devSYCL(sycl::device &d, const int qIndex, int id, const bool isInOrder, const bool isAsync)
{
	_exHandler = exceptionHandler;

	try
	{
		if (isInOrder)
		{
			sycl::property_list q_prop {sycl::property::queue::in_order{}, sycl::ext::intel::property::queue::compute_index{qIndex}};
			_devQ = new sycl::queue(d, _exHandler, q_prop);
			_queueType = syclQueueType::inOrder;
		}
		else
		{
			sycl::property_list q_prop {sycl::ext::intel::property::queue::compute_index{qIndex}};
			_devQ = new sycl::queue(d, _exHandler, q_prop);
			_queueType = syclQueueType::outOfOrder;
		}
	}
	catch (const sycl::exception &e)
	{
		std::cerr << "Provided SYCL device failed\n" << e.what() << std::endl;
		std::terminate();
	}

	_devD = d;
	_devC = _devQ->get_context();
	_syclHostPool = nullptr;
	_syclDevicePool = nullptr;
	_host_block_size = 0;
	_device_block_size = 0;
#ifdef USE_ONEDPL
	_devicePolicy = oneapi::dpl::execution::make_device_policy(*_devQ);
#endif
	_isFP64Supported = d.has(sycl::aspect::fp64);
	_isAsyncQueue = isAsync;
	_computeIndex = qIndex;

	const auto isizes = d.get_info<sycl::info::device::max_work_item_sizes<3>>();
	maxItem[0] = isizes[0];
	maxItem[1] = isizes[1];
	maxItem[2] = isizes[2];
	maxGroup = d.get_info<sycl::info::device::max_work_group_size>();
#ifdef SYCL_EXT_ONEAPI_MAX_WORK_GROUP_QUERY
	auto groups = d.get_info<sycl::ext::oneapi::experimental::info::device::max_work_groups<3>>();
	maxWorkGroup[0] = groups[0];
	maxWorkGroup[1] = groups[1];
	maxWorkGroup[2] = groups[2];
	maxGlobalWorkGroup = d.get_info<sycl::ext::oneapi::experimental::info::device::max_global_work_groups>();
#else
	maxWorkGroup[0] = maxWorkGroup[1] = maxWorkGroup[2] = maxGlobalWorkGroup = std::numeric_limits<int>::max;
#endif
	globalMem = d.get_info<sycl::info::device::global_mem_size>();
	localMem = d.get_info<sycl::info::device::local_mem_size>();
	maxUnit = d.get_info<sycl::info::device::max_compute_units>();
	deviceID = id;
	auto pf = d.get_platform();
	deviceName = d.get_info<sycl::info::device::name>() + " (" + d.get_info<sycl::info::device::driver_version>() + ") / " + pf.get_info<sycl::info::platform::name>();
	_prev_submission = sycl::event();
	_event.push_back(sycl::event());
}

devSYCL::devSYCL(sycl::device &d, const int qIndex, const syclQueueType qType, int id, const bool isAsync)
{
	_exHandler = exceptionHandler;

	try
	{
		switch (qType)
		{
			case syclQueueType::inOrder :
			{
				sycl::property_list q_prop {sycl::property::queue::in_order{}, sycl::ext::intel::property::queue::compute_index{qIndex}};
				_devQ = new sycl::queue(d, _exHandler, q_prop);
			}
				break;

			case syclQueueType::enableProfiling :
			{
				sycl::property_list q_prop {sycl::property::queue::enable_profiling{}, sycl::ext::intel::property::queue::compute_index{qIndex}};
				_devQ = new sycl::queue(d, _exHandler, q_prop);
			}
				break;

			case syclQueueType::outOfOrder :
			default :
			{
				sycl::property_list q_prop {sycl::ext::intel::property::queue::compute_index{qIndex}};
				_devQ = new sycl::queue(d, _exHandler, q_prop);
			}
				break;
		}
	}
	catch (const sycl::exception &e)
	{
		std::cerr << "Provided SYCL device failed\n" << e.what() << std::endl;
		std::terminate();
	}
	_devD = d;
	_devC = _devQ->get_context();
	_syclHostPool = nullptr;
	_syclDevicePool = nullptr;
	_host_block_size = 0;
	_device_block_size = 0;
	_queueType = qType;
#ifdef USE_ONEDPL
	_devicePolicy = oneapi::dpl::execution::make_device_policy(*_devQ);
#endif
	_isFP64Supported = d.has(sycl::aspect::fp64);
	_isAsyncQueue = isAsync;
	_computeIndex = qIndex;

	const auto isizes = d.get_info<sycl::info::device::max_work_item_sizes<3>>();
	maxItem[0] = isizes[0];
	maxItem[1] = isizes[1];
	maxItem[2] = isizes[2];
	maxGroup = d.get_info<sycl::info::device::max_work_group_size>();
#ifdef SYCL_EXT_ONEAPI_MAX_WORK_GROUP_QUERY
	auto groups = d.get_info<sycl::ext::oneapi::experimental::info::device::max_work_groups<3>>();
	maxWorkGroup[0] = groups[0];
	maxWorkGroup[1] = groups[1];
	maxWorkGroup[2] = groups[2];
	maxGlobalWorkGroup = d.get_info<sycl::ext::oneapi::experimental::info::device::max_global_work_groups>();
#else
	maxWorkGroup[0] = maxWorkGroup[1] = maxWorkGroup[2] = maxGlobalWorkGroup = std::numeric_limits<int>::max;
#endif
	globalMem = d.get_info<sycl::info::device::global_mem_size>();
	localMem = d.get_info<sycl::info::device::local_mem_size>();
	maxUnit = d.get_info<sycl::info::device::max_compute_units>();
	deviceID = id;
	auto pf = d.get_platform();
	deviceName = d.get_info<sycl::info::device::name>() + " (" + d.get_info<sycl::info::device::driver_version>() + ") / " + pf.get_info<sycl::info::platform::name>();
	_prev_submission = sycl::event();
	_event.push_back(sycl::event());
}

devSYCL::devSYCL(sycl::context &c, sycl::device &d, const int qIndex, int id, const bool isInOrder, const bool isAsync)
{
	_exHandler = exceptionHandler;

	try
	{
		if (isInOrder)
		{
			sycl::property_list q_prop {sycl::property::queue::in_order{}, sycl::ext::intel::property::queue::compute_index{qIndex}};
			_devQ = new sycl::queue(c, d, _exHandler, q_prop);
			_queueType = syclQueueType::inOrder;
		}
		else
		{
			sycl::property_list q_prop {sycl::ext::intel::property::queue::compute_index{qIndex}};
			_devQ = new sycl::queue(c, d, _exHandler, q_prop);
			_queueType = syclQueueType::outOfOrder;
		}
	}
	catch (const sycl::exception &e)
	{
		std::cerr << "Provided SYCL device failed\n" << e.what() << std::endl;
		std::terminate();
	}
	_devD = d;
	_devC = c;
	_syclHostPool = nullptr;
	_syclDevicePool = nullptr;
	_host_block_size = 0;
	_device_block_size = 0;
#ifdef USE_ONEDPL
	_devicePolicy = oneapi::dpl::execution::make_device_policy(*_devQ);
#endif
	_isFP64Supported = d.has(sycl::aspect::fp64);
	_isAsyncQueue = isAsync;
	_computeIndex = qIndex;

	const auto isizes = d.get_info<sycl::info::device::max_work_item_sizes<3>>();
	maxItem[0] = isizes[0];
	maxItem[1] = isizes[1];
	maxItem[2] = isizes[2];
	maxGroup = d.get_info<sycl::info::device::max_work_group_size>();
#ifdef SYCL_EXT_ONEAPI_MAX_WORK_GROUP_QUERY
	auto groups = d.get_info<sycl::ext::oneapi::experimental::info::device::max_work_groups<3>>();
	maxWorkGroup[0] = groups[0];
	maxWorkGroup[1] = groups[1];
	maxWorkGroup[2] = groups[2];
	maxGlobalWorkGroup = d.get_info<sycl::ext::oneapi::experimental::info::device::max_global_work_groups>();
#else
	maxWorkGroup[0] = maxWorkGroup[1] = maxWorkGroup[2] = maxGlobalWorkGroup = std::numeric_limits<int>::max;
#endif
	globalMem = d.get_info<sycl::info::device::global_mem_size>();
	localMem = d.get_info<sycl::info::device::local_mem_size>();
	maxUnit = d.get_info<sycl::info::device::max_compute_units>();
	deviceID = id;
	auto pf = d.get_platform();
	deviceName = d.get_info<sycl::info::device::name>() + " (" + d.get_info<sycl::info::device::driver_version>() + ") / " + pf.get_info<sycl::info::platform::name>();
	_prev_submission = sycl::event();
	_event.push_back(sycl::event());
}

devSYCL::devSYCL(sycl::context &c, sycl::device &d, const int qIndex, const syclQueueType qType, int id, const bool isAsync)
{
	_exHandler = exceptionHandler;

	try
	{
		switch (qType)
		{
			case syclQueueType::inOrder :
			{
				sycl::property_list q_prop {sycl::property::queue::in_order{}, sycl::ext::intel::property::queue::compute_index{qIndex}};
				_devQ = new sycl::queue(c, d, _exHandler, q_prop);
			}
				break;

			case syclQueueType::enableProfiling :
			{
				sycl::property_list q_prop {sycl::property::queue::enable_profiling{}, sycl::ext::intel::property::queue::compute_index{qIndex}};
				_devQ = new sycl::queue(c, d, _exHandler, q_prop);
			}
				break;

			case syclQueueType::outOfOrder :
			default :
			{
				sycl::property_list q_prop {sycl::ext::intel::property::queue::compute_index{qIndex}};
				_devQ = new sycl::queue(c, d, _exHandler, q_prop);
			}
				break;
		}
	}
	catch (const sycl::exception &e)
	{
		std::cerr << "Provided SYCL device failed\n" << e.what() << std::endl;
		std::terminate();
	}
	_devD = d;
	_devC = c;
	_syclHostPool = nullptr;
	_syclDevicePool = nullptr;
	_host_block_size = 0;
	_device_block_size = 0;
	_queueType = qType;
#ifdef USE_ONEDPL
	_devicePolicy = oneapi::dpl::execution::make_device_policy(*_devQ);
#endif
	_isFP64Supported = d.has(sycl::aspect::fp64);
	_isAsyncQueue = isAsync;
	_computeIndex = qIndex;

	const auto isizes = d.get_info<sycl::info::device::max_work_item_sizes<3>>();
	maxItem[0] = isizes[0];
	maxItem[1] = isizes[1];
	maxItem[2] = isizes[2];
	maxGroup = d.get_info<sycl::info::device::max_work_group_size>();
#ifdef SYCL_EXT_ONEAPI_MAX_WORK_GROUP_QUERY
	auto groups = d.get_info<sycl::ext::oneapi::experimental::info::device::max_work_groups<3>>();
	maxWorkGroup[0] = groups[0];
	maxWorkGroup[1] = groups[1];
	maxWorkGroup[2] = groups[2];
	maxGlobalWorkGroup = d.get_info<sycl::ext::oneapi::experimental::info::device::max_global_work_groups>();
#else
	maxWorkGroup[0] = maxWorkGroup[1] = maxWorkGroup[2] = maxGlobalWorkGroup = std::numeric_limits<int>::max;
#endif
	globalMem = d.get_info<sycl::info::device::global_mem_size>();
	localMem = d.get_info<sycl::info::device::local_mem_size>();
	maxUnit = d.get_info<sycl::info::device::max_compute_units>();
	deviceID = id;
	auto pf = d.get_platform();
	deviceName = d.get_info<sycl::info::device::name>() + " (" + d.get_info<sycl::info::device::driver_version>() + ") / " + pf.get_info<sycl::info::platform::name>();
	_prev_submission = sycl::event();
	_event.push_back(sycl::event());
}
#endif

devSYCL::devSYCL(const syclDeviceType dev, int id, const bool isInOrder, const bool isAsync)
{
	_exHandler = exceptionHandler;

	try
	{
		switch (dev)
		{
			case syclDeviceType::gpu :
				if (isInOrder)
				{
					_devQ = new sycl::queue(sycl::gpu_selector_v, _exHandler, sycl::property::queue::in_order{});
					_queueType = syclQueueType::inOrder;
				}
				else
				{
					_devQ = new sycl::queue(sycl::gpu_selector_v, _exHandler);
					_queueType = syclQueueType::outOfOrder;
				}
				break;

			case syclDeviceType::cpu :
				_devQ = new sycl::queue(sycl::cpu_selector_v, _exHandler);
				_queueType = syclQueueType::outOfOrder;
				break;

			case syclDeviceType::host :
			case syclDeviceType::fpga :
			default:
				_devQ = new sycl::queue(sycl::default_selector_v, _exHandler);
				_queueType = syclQueueType::outOfOrder;
				break;
		}
	}
	catch (const sycl::exception &e)
	{
		std::cerr << "There is no available SYCL device.\n" << e.what() << std::endl;
		std::terminate();
	}

	_devD = _devQ->get_device();
	_devC = _devQ->get_context();
	_syclHostPool = nullptr;
	_syclDevicePool = nullptr;
	_host_block_size = 0;
	_device_block_size = 0;
#ifdef USE_ONEDPL
	_devicePolicy = oneapi::dpl::execution::make_device_policy(*_devQ);
#endif
	_isFP64Supported = _devD.has(sycl::aspect::fp64);
	_isAsyncQueue = isAsync;
	_computeIndex = -1;

	auto d = _devD;
	auto pf = d.get_platform();
	const auto isizes = d.get_info<sycl::info::device::max_work_item_sizes<3>>();
	maxItem[0] = isizes[0];
	maxItem[1] = isizes[1];
	maxItem[2] = isizes[2];
	maxGroup = d.get_info<sycl::info::device::max_work_group_size>();
#ifdef SYCL_EXT_ONEAPI_MAX_WORK_GROUP_QUERY
	auto groups = d.get_info<sycl::ext::oneapi::experimental::info::device::max_work_groups<3>>();
	maxWorkGroup[0] = groups[0];
	maxWorkGroup[1] = groups[1];
	maxWorkGroup[2] = groups[2];
	maxGlobalWorkGroup = d.get_info<sycl::ext::oneapi::experimental::info::device::max_global_work_groups>();
#else
	maxWorkGroup[0] = maxWorkGroup[1] = maxWorkGroup[2] = maxGlobalWorkGroup = std::numeric_limits<int>::max;
#endif
	globalMem = d.get_info<sycl::info::device::global_mem_size>();
	localMem = d.get_info<sycl::info::device::local_mem_size>();
	maxUnit = d.get_info<sycl::info::device::max_compute_units>();
	deviceID = id;
	deviceName = d.get_info<sycl::info::device::name>() + " (" + d.get_info<sycl::info::device::driver_version>() + ") / " + pf.get_info<sycl::info::platform::name>();
	_prev_submission = sycl::event();
	_event.push_back(sycl::event());
}

devSYCL::devSYCL(const syclBackendType be, const syclDeviceType dev, int id, const bool isInOrder, const bool isAsync)
{
	_exHandler = exceptionHandler;

	try
	{
		switch (dev)
		{
			case syclDeviceType::gpu :
				if (isInOrder)
				{
					_devQ = new sycl::queue(sycl::gpu_selector_v, _exHandler, sycl::property::queue::in_order{});
					_queueType = syclQueueType::inOrder;
				}
				else
				{
					_devQ = new sycl::queue(sycl::gpu_selector_v, _exHandler);
					_queueType = syclQueueType::outOfOrder;
				}
				break;

			case syclDeviceType::cpu :
				_devQ = new sycl::queue(sycl::cpu_selector_v, _exHandler);
				_queueType = syclQueueType::outOfOrder;
				break;

			case syclDeviceType::host :
			case syclDeviceType::fpga :
			default:
				_devQ = new sycl::queue(sycl::default_selector_v, _exHandler);
				_queueType = syclQueueType::outOfOrder;
				break;
			}
	}
	catch (const sycl::exception &e)
	{
		std::cerr << "There is no available SYCL device.\n" << e.what() << std::endl;
		std::terminate();
	}

	_devD = _devQ->get_device();
	_devC = _devQ->get_context();
	_syclHostPool = nullptr;
	_syclDevicePool = nullptr;
	_host_block_size = 0;
	_device_block_size = 0;
#ifdef USE_ONEDPL
	_devicePolicy = oneapi::dpl::execution::make_device_policy(*_devQ);
#endif
	_isFP64Supported = _devD.has(sycl::aspect::fp64);
	_isAsyncQueue = isAsync;
	_computeIndex = -1;

	auto d = _devD;
	auto pf = d.get_platform();
	const auto isizes = d.get_info<sycl::info::device::max_work_item_sizes<3>>();
	maxItem[0] = isizes[0];
	maxItem[1] = isizes[1];
	maxItem[2] = isizes[2];
	maxGroup = d.get_info<sycl::info::device::max_work_group_size>();
#ifdef SYCL_EXT_ONEAPI_MAX_WORK_GROUP_QUERY
	auto groups = d.get_info<sycl::ext::oneapi::experimental::info::device::max_work_groups<3>>();
	maxWorkGroup[0] = groups[0];
	maxWorkGroup[1] = groups[1];
	maxWorkGroup[2] = groups[2];
	maxGlobalWorkGroup = d.get_info<sycl::ext::oneapi::experimental::info::device::max_global_work_groups>();
#else
	maxWorkGroup[0] = maxWorkGroup[1] = maxWorkGroup[2] = maxGlobalWorkGroup = std::numeric_limits<int>::max;
#endif
	globalMem = d.get_info<sycl::info::device::global_mem_size>();
	localMem = d.get_info<sycl::info::device::local_mem_size>();
	maxUnit = d.get_info<sycl::info::device::max_compute_units>();
	deviceID = id;
	deviceName = d.get_info<sycl::info::device::name>() + " (" + d.get_info<sycl::info::device::driver_version>() + ") / " + pf.get_info<sycl::info::platform::name>();
	_prev_submission = sycl::event();
	_event.push_back(sycl::event());
}

devSYCL::devSYCL(const devSYCL &q)
:	_devQ {q._devQ}, _devD {q._devD}, _devC {q._devC}, deviceName {q.deviceName}, _exHandler {q._exHandler},
	_syclHostPool {q._syclHostPool}, _syclDevicePool {q._syclDevicePool}, _host_block_size {q._host_block_size}, _device_block_size {q._device_block_size},
	_isAsyncQueue {q._isAsyncQueue}, _isFP64Supported {q._isFP64Supported}, _computeIndex {q._computeIndex},
	cardID {q.cardID}, deviceID {q.deviceID}, stackID {q.stackID}, nStack {q.nStack}, sliceID {q.sliceID}, nSlice {q.nSlice},
	_queueType {q._queueType}, maxItem {q.maxItem}, maxGroup {q.maxGroup}, maxWorkGroup {q.maxWorkGroup}, maxGlobalWorkGroup {q.maxGlobalWorkGroup}, maxUnit {q.maxUnit},
	globalMem {q.globalMem}, localMem {q.localMem}, _prev_submission {sycl::event()}, _event {sycl::event()}
#ifdef USE_ONEDPL
	, _devicePolicy {q._devicePolicy}
#endif
{}

devSYCL::devSYCL(const devSYCL *q)
:	_devQ {q->_devQ}, _devD {q->_devD}, _devC {q->_devC}, deviceName {q->deviceName}, _exHandler {q->_exHandler},
	_syclHostPool {q->_syclHostPool}, _syclDevicePool {q->_syclDevicePool}, _host_block_size {q->_host_block_size}, _device_block_size {q->_device_block_size},
	_isAsyncQueue {q->_isAsyncQueue}, _isFP64Supported {q->_isFP64Supported}, _computeIndex {q->_computeIndex},
	cardID {q->cardID}, deviceID {q->deviceID}, stackID {q->stackID}, nStack {q->nStack}, sliceID {q->sliceID}, nSlice {q->nSlice},
	_queueType {q->_queueType}, maxItem {q->maxItem}, maxGroup {q->maxGroup}, maxWorkGroup {q->maxWorkGroup}, maxGlobalWorkGroup {q->maxGlobalWorkGroup}, maxUnit {q->maxUnit},
	globalMem {q->globalMem}, localMem {q->localMem}, _prev_submission {sycl::event()}, _event {sycl::event()}
#ifdef USE_ONEDPL
	, _devicePolicy {q->_devicePolicy}
#endif
{}

devSYCL::~devSYCL()
{
	waitAll();
	delete _devQ;
}

void devSYCL::destroyMemoryPool()
{
	if (_syclHostPool)
		delete _syclHostPool;
	if (_syclDevicePool)
		delete _syclDevicePool;
}

void* devSYCL::syclMalloc(const size_t bytes, const syclMallocType type, const char *name)
{
	assert(bytes > 0);
	void* ptr = nullptr;
	if (type == syclMallocType::device)
	{
		if (_syclDevicePool != nullptr && bytes <= _device_block_size)
			ptr = _syclDevicePool->alloc(bytes);
		if (_syclDevicePool == nullptr || ptr == nullptr)
		{
			std::lock_guard<std::mutex> locker(_lock);
			ptr = sycl::malloc_device(bytes, *_devQ);
		}
	}
	else if (type == syclMallocType::host)
	{
		if (_syclHostPool != nullptr && bytes <= _host_block_size)
			ptr = _syclHostPool->alloc(bytes);
		if (_syclHostPool == nullptr || ptr == nullptr)
		{
			std::lock_guard<std::mutex> locker(_lock);
			ptr = sycl::malloc_host(bytes, *_devQ);
		}
	}
	assert(ptr != nullptr);

	return ptr;
}

void devSYCL::syclFree(void *ptr)
{
	if (ptr)
	{
		bool isFree = false;
		auto kind = sycl::get_pointer_type(ptr, _devC);
		if (_syclDevicePool && kind == alloc_kind::device)
			isFree = _syclDevicePool->free(ptr);
		else if (_syclHostPool && kind == alloc_kind::host)
			isFree = _syclHostPool->free(ptr);
		if (isFree == false && (kind == alloc_kind::device || kind == alloc_kind::host))
		{
			std::lock_guard<std::mutex> locker(_lock);
			sycl::free(ptr, *_devQ);
		}
		else if (isFree == false)	// This should not happen
			free(ptr);
	}
}

void devSYCL::syclMemcpy(void *dest, const void *src, const size_t bytes)
{
	assert(dest != nullptr);
	assert( src != nullptr);
	assert(bytes > 0);

	pushEvent(_devQ->memcpy(dest, src, bytes));
}

void devSYCL::syclMemcpyAfterWaitAll(void *dest, const void *src, const size_t bytes)
{
	assert(dest != nullptr);
	assert( src != nullptr);
	assert(bytes > 0);

	waitAll();
	pushEvent(_devQ->memcpy(dest, src, bytes));
}

void devSYCL::syclMemset(void *ptr, const int value, const size_t bytes)
{
	assert(ptr != nullptr);
	assert(bytes > 0);

	pushEvent(_devQ->memset(ptr, value, bytes));
}

void devSYCL::syclMemsetAfterWaitAll(void *ptr, const int value, const size_t bytes)
{
	assert(ptr != nullptr);
	assert(bytes > 0);

	waitAll();
	pushEvent(_devQ->memset(ptr, value, bytes));
}

void devSYCL::syclPrefetch(void *ptr, const size_t bytes)
{
	assert(ptr != nullptr);
	assert(bytes > 0);

	pushEvent(_devQ->prefetch(ptr, bytes));
}

void devSYCL::syclPrefetchAfterWaitAll(void *ptr, const size_t bytes)
{
	assert(ptr != nullptr);
	assert(bytes > 0);

	waitAll();
	pushEvent(_devQ->prefetch(ptr, bytes));
}

void devSYCL::syclFillInt32(void *ptr, const std::int32_t value, const size_t count)
{
	return syclFill<std::int32_t>(ptr, value, count);
}

void devSYCL::syclFillInt32AfterWaitAll(void *ptr, const std::int32_t value, const size_t count)
{
	return syclFillAfterWaitAll<std::int32_t>(ptr, value, count);
}

void devSYCL::syclFillUint32(void *ptr, const std::uint32_t value, const size_t count)
{
	return syclFill<std::uint32_t>(ptr, value, count);
}

void devSYCL::syclFillUint32AfterWaitAll(void *ptr, const std::uint32_t value, const size_t count)
{
	return syclFillAfterWaitAll<std::uint32_t>(ptr, value, count);
}

void devSYCL::syclFillInt64(void *ptr, const std::int64_t value, const size_t count)
{
	return syclFill<std::int64_t>(ptr, value, count);
}

void devSYCL::syclFillInt64AfterWaitAll(void *ptr, const std::int64_t value, const size_t count)
{
	return syclFillAfterWaitAll<std::int64_t>(ptr, value, count);
}

void devSYCL::syclFillUint64(void *ptr, const std::uint64_t value, const size_t count)
{
	return syclFill<std::uint64_t>(ptr, value, count);
}

void devSYCL::syclFillUint64AfterWaitAll(void *ptr, const std::uint64_t value, const size_t count)
{
	return syclFillAfterWaitAll<std::uint64_t>(ptr, value, count);
}

void devSYCL::syclFillXfloat(void *ptr, const XFLOAT value, const size_t count)
{
	return syclFill<XFLOAT>(ptr, value, count);
}

void devSYCL::syclFillXfloatAfterWaitAll(void *ptr, const XFLOAT value, const size_t count)
{
	return syclFillAfterWaitAll<XFLOAT>(ptr, value, count);
}

void devSYCL::syclFillRfloat(void *ptr, const RFLOAT value, const size_t count)
{
	return syclFill<RFLOAT>(ptr, value, count);
}

void devSYCL::syclFillRfloatAfterWaitAll(void *ptr, const RFLOAT value, const size_t count)
{
	return syclFillAfterWaitAll<RFLOAT>(ptr, value, count);
}

void devSYCL::syclFillFloat(void *ptr, const float value, const size_t count)
{
	return syclFill<float>(ptr, value, count);
}

void devSYCL::syclFillFloatAfterWaitAll(void *ptr, const float value, const size_t count)
{
	return syclFillAfterWaitAll<float>(ptr, value, count);
}

void devSYCL::syclFillDouble(void *ptr, const double value, const size_t count)
{
	return syclFill<double>(ptr, value, count);
}

void devSYCL::syclFillDoubleAfterWaitAll(void *ptr, const double value, const size_t count)
{
	return syclFillAfterWaitAll<double>(ptr, value, count);
}

void devSYCL::waitAll()
{
	if (_queueType == syclQueueType::inOrder)
	{
		_prev_submission.wait_and_throw();
		_prev_submission = sycl::event();
	}
	else
	{
		sycl::event::wait_and_throw(_event);
		_event.clear();
		_event.push_back(sycl::event());
	}
}

void devSYCL::reCalculateRange(sycl::range<3> &wg, const sycl::range<3> &wi)
{
	if (wg[2] > wi[2]*maxWorkGroup[2])
	{
		auto wg12 = wg[1]*wg[2];
		wg[2] = wi[2] * maxWorkGroup[2];
		wg[1] = wg12 % wg[2] == 0 ? wg12 / wg[2] : wg12 / wg[2] + 1;
	}
	if (wg[1] > wi[1]*maxWorkGroup[1])
	{
		auto wg01 = wg[0]*wg[1];
		wg[1] = wi[1] * maxWorkGroup[1];
		wg[0] = wg01 % wg[1] == 0 ? wg01 / wg[1] : wg01 / wg[1] + 1;
	}
#if 0
	if (wg[0] > wi[0]*maxWorkGroup[0])
	{
		std::cerr << "The number of work-groups are too big. Please increase work-item.\n";
	}
#endif
}

std::string devSYCL::getName()
{
	std::stringstream ss;

	if (cardID < 0)
		ss << "[" << deviceID << "] " << deviceName;
	else if (nStack > 1)
	{
		if (nSlice < 0)
			ss << "[" << deviceID << "] Card[" << cardID << "]Stack[" << stackID << "] / " << deviceName;
		else if (sliceID >= 0)
			ss << "[" << deviceID << "] Card[" << cardID << "]Stack[" << stackID << "]{" << sliceID << "} / " << deviceName;
		else if (sliceID < 0)
			ss << "[" << deviceID << "] Card[" << cardID << "]Stack[" << stackID << "]{0-" << nSlice-1 << "} / " << deviceName;
	}
	else
	{
		if (nSlice < 0)
			ss << "[" << deviceID << "] Card[" << cardID << "] / " << deviceName;
		else if (sliceID >= 0)
			ss << "[" << deviceID << "] Card[" << cardID << "]{" << sliceID << "} / " << deviceName;
		else if (sliceID < 0)
			ss << "[" << deviceID << "] Card[" << cardID << "]{0-" << nSlice-1 << "} / " << deviceName;
	}

	return ss.str();
}

void devSYCL::printDeviceInfo(bool printAll)
{
	auto d = _devQ->get_device();
	PRINT_DEV_INFO(d, name);
	PRINT_DEV_INFO(d, version);
	PRINT_DEV_INFO(d, driver_version);
	PRINT_DEV_INFO(d, partition_max_sub_devices);

	PRINT_DEV_INFO(d, max_compute_units);
#ifdef SYCL_EXT_ONEAPI_MAX_WORK_GROUP_QUERY
	std::cout << "  Max global number of work-groups: " << maxGlobalWorkGroup << std::endl;
	std::cout << "  Max number of work-groups: " << maxWorkGroup[0] << " x " << maxWorkGroup[1] << " x " << maxWorkGroup[2] << std::endl;
#endif
	PRINT_DEV_INFO(d, max_work_group_size);
	PRINT_DEV_INFO(d, max_work_item_dimensions);
	std::cout << "  max_work_item_sizes: " << maxItem[0] << " x " << maxItem[1] << " x " << maxItem[2] << std::endl;
	PRINT_DEV_INFO(d, max_num_sub_groups);
	const auto subGroups = d.get_info<sycl::info::device::sub_group_sizes>();
	std::cout << "  sub_group_sizes: ";
	for (auto &subg : subGroups)
		std::cout << subg << " ";
	std::cout << std::endl;

	PRINT_DEV_INFO(d, max_mem_alloc_size);
	PRINT_DEV_INFO(d, global_mem_cache_line_size);
	PRINT_DEV_INFO(d, global_mem_cache_size);
	PRINT_DEV_INFO(d, global_mem_size);
	PRINT_DEV_INFO(d, local_mem_size);

	if (printAll)
	{
		PRINT_DEV_INFO(d, profiling_timer_resolution);
		PRINT_DEV_INFO2020(d, queue_profiling);

		PRINT_DEV_INFO(d, preferred_vector_width_int);
		PRINT_DEV_INFO(d, preferred_vector_width_long);
		PRINT_DEV_INFO(d, preferred_vector_width_half);
		PRINT_DEV_INFO(d, preferred_vector_width_float);
		PRINT_DEV_INFO(d, preferred_vector_width_double);
		PRINT_DEV_INFO(d, native_vector_width_int);
		PRINT_DEV_INFO(d, native_vector_width_long);
		PRINT_DEV_INFO(d, native_vector_width_half);
		PRINT_DEV_INFO(d, native_vector_width_float);
		PRINT_DEV_INFO(d, native_vector_width_double);

		PRINT_DEV_INFO2020(d, usm_device_allocations);
		PRINT_DEV_INFO2020(d, usm_host_allocations);

		const auto domains = d.get_info<sycl::info::device::partition_affinity_domains>();
		std::cout << "  partition_affinity_domain:";
		for (auto& domain : domains)
		{
			switch(domain)
			{
				case sycl::info::partition_affinity_domain::numa :
					std::cout << " numa";
					break;
				case sycl::info::partition_affinity_domain::L1_cache :
					std::cout << " L1_cache";
					break;
				case sycl::info::partition_affinity_domain::L2_cache :
					std::cout << " L2_cache";
					break;
				case sycl::info::partition_affinity_domain::L3_cache :
					std::cout << " L3_cache";
					break;
				default :
					break;
			}
		}
		std::cout << std::endl;
	}

	const auto hconfigs = d.get_info<sycl::info::device::half_fp_config>();
	std::cout << "  half_fp_config:";
	for (auto& hconfig : hconfigs)
	{
		switch(hconfig)
		{
			case sycl::info::fp_config::fma :
				std::cout << " fma";
				break;
			case sycl::info::fp_config::denorm :
				std::cout << " denorm";
				break;
			case sycl::info::fp_config::inf_nan :
				std::cout << " inf_nan";
				break;
			case sycl::info::fp_config::round_to_nearest :
				std::cout << " round_to_nearest";
				break;
			case sycl::info::fp_config::round_to_zero :
				std::cout << " round_to_zero";
				break;
			case sycl::info::fp_config::round_to_inf :
				std::cout << " round_to_inf";
				break;
			case sycl::info::fp_config::correctly_rounded_divide_sqrt :
				std::cout << " correctly_rounded_divide_sqrt";
				break;
			case sycl::info::fp_config::soft_float :
				std::cout << " soft_float";
				break;
			default :
				break;
		}
	}
	std::cout << std::endl;

	const auto sconfigs = d.get_info<sycl::info::device::single_fp_config>();
	std::cout << "  single_fp_config:";
	for (auto& sconfig : sconfigs)
	{
		switch(sconfig)
		{
			case sycl::info::fp_config::fma :
				std::cout << " fma";
				break;
			case sycl::info::fp_config::denorm :
				std::cout << " denorm";
				break;
			case sycl::info::fp_config::inf_nan :
				std::cout << " inf_nan";
				break;
			case sycl::info::fp_config::round_to_nearest :
				std::cout << " round_to_nearest";
				break;
			case sycl::info::fp_config::round_to_zero :
				std::cout << " round_to_zero";
				break;
			case sycl::info::fp_config::round_to_inf :
				std::cout << " round_to_inf";
				break;
			case sycl::info::fp_config::correctly_rounded_divide_sqrt :
				std::cout << " correctly_rounded_divide_sqrt";
				break;
			case sycl::info::fp_config::soft_float :
				std::cout << " soft_float";
				break;
			default :
				break;
		}
	}
	std::cout << std::endl;

	const auto dconfigs = d.get_info<sycl::info::device::double_fp_config>();
	std::cout << "  double_fp_config:";
	for (auto& dconfig : dconfigs)
	{
		switch(dconfig)
		{
			case sycl::info::fp_config::fma :
				std::cout << " fma";
				break;
			case sycl::info::fp_config::denorm :
				std::cout << " denorm";
				break;
			case sycl::info::fp_config::inf_nan :
				std::cout << " inf_nan";
				break;
			case sycl::info::fp_config::round_to_nearest :
				std::cout << " round_to_nearest";
				break;
			case sycl::info::fp_config::round_to_zero :
				std::cout << " round_to_zero";
				break;
			case sycl::info::fp_config::round_to_inf :
				std::cout << " round_to_inf";
				break;
			case sycl::info::fp_config::correctly_rounded_divide_sqrt :
				std::cout << " correctly_rounded_divide_sqrt";
				break;
			case sycl::info::fp_config::soft_float :
				std::cout << " soft_float";
				break;
			default :
				break;
		}
	}
	std::cout << "\n";
}

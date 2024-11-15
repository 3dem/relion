#ifdef _SYCL_ENABLED
// A large amount of this code is direct from cuda_ml_optimizer and so could
// be shared (but possibly with difficulty since it is enough different that
// we either need a lot of #ifdefs, or a lot of macros/some other mechanism to
// abstract the differences). The biggest differences are the type of memory
// objects used (std::vector vs. CudaGlobalPtr and CudaCustomAllocator), the
// lack of transfers to/from the device, and on-device operations (which are
// replaced by loops/function calls).
//
// CudaFFT has been replaced with lib FFTW, if RELION is configured with mix
// precision, both single and double precision FFTW are linked into RELION.
// Install fftw-static.x86_64 and fftw-static.i686 to get the libraries without
// having to pull them at build time. Over time we hope to replace FFTW with
// MKL.
//
// NOTE: Since the GPU code was ported back to CPU there may be additional
// changes made in the CUDA code which may not have made it here.

#include <cstdio>
#include <ctime>
#include <cmath>
#include <csignal>
#include <ctime>
#include <vector>
#include <iostream>
#include <fstream>
#include <map>
#include <sycl/sycl.hpp>
#include <sycl/backend.hpp>

#include "src/acc/sycl/device_stubs.h"
#include "src/ml_optimiser.h"
#include "src/acc/acc_ptr.h"
#include "src/acc/acc_projector.h"
#include "src/acc/acc_projector_plan.h"
#include "src/acc/sycl/sycl_benchmark_utils.h"
#include "src/acc/sycl/sycl_helper_functions.h"
#include "src/acc/sycl/mkl_fft.h"
#include "src/acc/data_types.h"
#include "src/parallel.h"

#include "src/acc/utilities.h"
#include "src/acc/utilities_impl.h"

#include "src/acc/acc_helper_functions.h"
#include "src/acc/acc_ml_optimiser.h"
#include "src/acc/acc_ml_optimiser_impl.h"
#include "src/acc/sycl/sycl_ml_optimiser.h"
#include "src/acc/sycl/sycl_memory_pool.h"

#include "src/acc/sycl/sycl_dev.h"

#define PRINT_SYCL_INFO(x) std::cout << "  "#x": " << device.get_info<sycl::info::device::x>() << std::endl

#ifdef SYCL_EXT_INTEL_CSLICE
template <typename RangeTy, typename ElemTy>
bool contains(RangeTy &&Range, const ElemTy &Elem) {
	return std::find(Range.begin(), Range.end(), Elem) != Range.end();
}

bool isPartitionableByCSlice(sycl::device &Dev) {
	return contains(Dev.get_info<sycl::info::device::partition_properties>(), sycl::info::partition_property::ext_intel_partition_by_cslice);
}
#else
bool isPartitionableByCSlice(sycl::device &Dev) {
	return false;
}
#endif

#ifdef SYCL_EXT_INTEL_QUEUE_INDEX
int numPartitionableByQueueIndex(sycl::device &Dev) {
	return Dev.get_info<sycl::ext::intel::info::device::max_compute_queue_indices>();
}
#else
int numPartitionableByQueueIndex(sycl::device &Dev) {
	return 1;
}
#endif

std::mutex	fft_mutex;

MlSyclDataBundle::MlSyclDataBundle(virtualSYCL *dev) : generateProjectionPlanOnTheFly {false}, _devAcc {dev}
{
}

MlSyclDataBundle::~MlSyclDataBundle()
{
	auto dev = dynamic_cast<devSYCL*>(_devAcc);
	if (dev->getSyclHostPool())
		dev->getSyclHostPool()->shrink();
	if (dev->getSyclDevicePool())
		dev->getSyclDevicePool()->shrink();

	projectors.clear();
	backprojectors.clear();
	coarseProjectionPlans.clear();
}

void MlSyclDataBundle::setup(MlOptimiser *baseMLO)
{
	auto dev = dynamic_cast<devSYCL*>(_devAcc);
	if (dev->getSyclHostPool() == nullptr || dev->getSyclDevicePool() == nullptr)
	{
		char* pEnvBlockSize = std::getenv("relionSyclBlockSize");
		const size_t def_block_size = (pEnvBlockSize == nullptr) ? 0 : std::strtoul(pEnvBlockSize, nullptr, 10);
		size_t block_size;
		if (dev->getSyclHostPool() == nullptr)
		{
			if (def_block_size > 0)
				block_size = def_block_size;
			else
			{
				char* pEnvHostBlockSize = std::getenv("relionSyclHostBlockSize");
				block_size = (pEnvHostBlockSize == nullptr) ? defHostBlockSize : std::strtoul(pEnvHostBlockSize, nullptr, 10);
			}
			dev->setSyclHostPool(new syclMemoryPool(dev->getQueue(), alloc_kind::host, block_size), block_size);
		}
		if (dev->getSyclDevicePool() == nullptr)
		{
			if (def_block_size > 0)
				block_size = def_block_size;
			else
			{
				char* pEnvDeviceBlockSize = std::getenv("relionSyclDeviceBlockSize");
				block_size = (pEnvDeviceBlockSize == nullptr) ? defDeviceBlockSize : std::strtoul(pEnvDeviceBlockSize, nullptr, 10);
			}
			dev->setSyclDevicePool(new syclMemoryPool(dev->getQueue(), alloc_kind::device, block_size), block_size);
		}
	}

	/*======================================================
				PROJECTOR AND BACKPROJECTOR
	======================================================*/

	unsigned nr_proj = baseMLO->mymodel.PPref.size();
	unsigned nr_bproj = baseMLO->wsum_model.BPref.size();

	projectors.resize(nr_proj);
	backprojectors.resize(nr_bproj);

	//Loop over classes
	for (int imodel = 0; imodel < nr_proj; imodel++)
	{
		projectors[imodel].setMdlDim(
				_devAcc,
				baseMLO->mymodel.PPref[imodel].data.xdim,
				baseMLO->mymodel.PPref[imodel].data.ydim,
				baseMLO->mymodel.PPref[imodel].data.zdim,
				baseMLO->mymodel.PPref[imodel].data.yinit,
				baseMLO->mymodel.PPref[imodel].data.zinit,
				baseMLO->mymodel.PPref[imodel].r_max,
				baseMLO->mymodel.PPref[imodel].padding_factor);

		projectors[imodel].initMdl(baseMLO->mymodel.PPref[imodel].data.data);
	}

	for (int imodel = 0; imodel < nr_bproj; imodel++)
	{
		backprojectors[imodel].setMdlDim(
				_devAcc,
				baseMLO->wsum_model.BPref[imodel].data.xdim,
				baseMLO->wsum_model.BPref[imodel].data.ydim,
				baseMLO->wsum_model.BPref[imodel].data.zdim,
				baseMLO->wsum_model.BPref[imodel].data.yinit,
				baseMLO->wsum_model.BPref[imodel].data.zinit,
				baseMLO->wsum_model.BPref[imodel].r_max,
				baseMLO->wsum_model.BPref[imodel].padding_factor);

		backprojectors[imodel].initMdl();
	}

	/*======================================================
						PROJECTION PLAN
	======================================================*/

	unsigned nr_classes = baseMLO->mymodel.nr_classes;
	coarseProjectionPlans.resize(nr_classes);

	//Can we pre-generate projector plan and corresponding euler matrices for all particles
	if (baseMLO->do_skip_align || baseMLO->do_skip_rotate || baseMLO->do_auto_refine || baseMLO->mymodel.orientational_prior_mode != NOPRIOR || baseMLO->mydata.is_tomo)
		generateProjectionPlanOnTheFly = true;
	else
		generateProjectionPlanOnTheFly = false;

	for (int iclass = 0; iclass < nr_classes; iclass++)
	{
		//If doing predefined projector plan at all and is this class significant
		if (!generateProjectionPlanOnTheFly && baseMLO->mymodel.pdf_class[iclass] > 0.)
		{
			std::vector<int> exp_pointer_dir_nonzeroprior;
			std::vector<int> exp_pointer_psi_nonzeroprior;
			std::vector<RFLOAT> exp_directions_prior;
			std::vector<RFLOAT> exp_psi_prior;

			long unsigned itrans_max = baseMLO->sampling.NrTranslationalSamplings() - 1;
			long unsigned nr_idir = baseMLO->sampling.NrDirections(0, &exp_pointer_dir_nonzeroprior);
			long unsigned nr_ipsi = baseMLO->sampling.NrPsiSamplings(0, &exp_pointer_psi_nonzeroprior );

#ifdef _SYCL_ENABLED
			coarseProjectionPlans[iclass].setSyclDevice(_devAcc);
#endif
			coarseProjectionPlans[iclass].setup(
					baseMLO->sampling,
					exp_directions_prior,
					exp_psi_prior,
					exp_pointer_dir_nonzeroprior,
					exp_pointer_psi_nonzeroprior,
					NULL, //Mcoarse_significant
					baseMLO->mymodel.pdf_class,
					baseMLO->mymodel.pdf_direction,
					nr_idir,
					nr_ipsi,
					0, //idir_min
					nr_idir - 1, //idir_max
					0, //ipsi_min
					nr_ipsi - 1, //ipsi_max
					0, //itrans_min
					itrans_max,
					0, //current_oversampling
					1, //nr_oversampled_rot
					iclass,
					true, //coarse
					!IS_NOT_INV,
					baseMLO->do_skip_align,
					baseMLO->do_skip_rotate,
					baseMLO->mymodel.orientational_prior_mode
					);
		}
	}
}

MlOptimiserSYCL::MlOptimiserSYCL(MlOptimiser *baseMLOptimiser, MlSyclDataBundle *b, const bool isStream, const char *timing_fnm) :
	baseMLO {baseMLOptimiser},
	bundle {b},
	transformer1 {baseMLOptimiser->mymodel.data_dim},
	transformer2 {baseMLOptimiser->mymodel.data_dim},
	refIs3D {baseMLO->mymodel.ref_dim == 3},
	dataIs3D {baseMLO->mymodel.data_dim == 3},
	shiftsIs3D {baseMLO->mymodel.data_dim == 3 || baseMLO->mydata.is_tomo},
	generateProjectionPlanOnTheFly {bundle->generateProjectionPlanOnTheFly},
	_useStream {isStream}
{
	if (baseMLOptimiser == nullptr)
		_devAcc = nullptr;
	else
		setupDevice();
}

MlOptimiserSYCL::~MlOptimiserSYCL()
{
	if (_useStream)
		for (auto q : classStreams)
			delete q;

	classStreams.clear();
#ifndef USE_EXISTING_SYCL_DEVICE
	if (_devAcc != nullptr)
		delete _devAcc;
#endif
}

std::vector<virtualSYCL*> MlOptimiserSYCL::getDevices(const syclDeviceType select, const std::tuple<bool,bool,bool> syclOpt, const syclBackendType BE, const bool verbose)
{
	bool isSubSub, isInOrderQueue, isAsyncSubmission;
	std::tie(isSubSub, isInOrderQueue, isAsyncSubmission) = syclOpt;
	std::vector<virtualSYCL*> selectedDevices;

	if (select == syclDeviceType::gpu)
	{
		int nDevice = 0, nCard = 0;
		auto devices = sycl::device::get_devices(sycl::info::device_type::gpu);
		for (auto &device : devices)
		{
			auto plName = device.get_platform().get_info<sycl::info::platform::name>();
			auto beType = device.get_backend();
			if ( (BE == syclBackendType::levelZero && beType == sycl::backend::ext_oneapi_level_zero)
				|| (BE == syclBackendType::CUDA && beType == sycl::backend::ext_oneapi_cuda)
				|| (BE == syclBackendType::HIP && beType == sycl::backend::ext_oneapi_hip)
				|| (BE == syclBackendType::openCL && beType == sycl::backend::opencl) )
			{
				auto ndev = device.get_info<sycl::info::device::partition_max_sub_devices>();
				if (ndev > 1)
				{
					auto subs = device.create_sub_devices<sycl::info::partition_property::partition_by_affinity_domain>(sycl::info::partition_affinity_domain::numa);
					auto ctx = sycl::context(subs);
					int nStack = 0;
					for (auto &sub : subs)
					{
#if defined(SYCL_EXT_INTEL_CSLICE) || defined(SYCL_EXT_INTEL_QUEUE_INDEX)
 #ifdef SYCL_EXT_INTEL_CSLICE
						if (isSubSub && isPartitionableByCSlice(sub))
						{
							auto subsubs = sub.create_sub_devices<sycl::info::partition_property::ext_intel_partition_by_cslice>();
							auto ctxctx = sycl::context(subsubs);
							int nSlice = 0;
							for (auto &subsub : subsubs)
							{
								selectedDevices.push_back(new devSYCL(ctxctx, subsub, nDevice++, isInOrderQueue, isAsyncSubmission));
								auto addedDevice = dynamic_cast<devSYCL*>(selectedDevices.back());
								addedDevice->setCardID(nCard);
								addedDevice->setStackID(nStack);
								addedDevice->setNumStack(subs.size());
								addedDevice->setSliceID(nSlice++);
								addedDevice->setNumSlice(1);
								if (verbose)
								{
									std::cout << std::string(80, '*') << std::endl;
									std::cout << "Created SYCL device is " << addedDevice->getName() << std::endl;
									std::cout << "maxComputeUnit= " << addedDevice->maxUnit << ", maxWorkGroupSize= " << addedDevice->maxGroup << ", globalMemSize= " << addedDevice->globalMem << std::endl;
									addedDevice->printDeviceInfo();
									std::cout << "\n";
								}
							}
						}
 #else
						int maxQindex = numPartitionableByQueueIndex(sub);
						if (isSubSub && maxQindex > 1)
						{
							int nSlice = 0;
							for (int i = 0; i < maxQindex; i++)
							{
								selectedDevices.push_back(new devSYCL(ctx, sub, i, nDevice++, isInOrderQueue, isAsyncSubmission));
								auto addedDevice = dynamic_cast<devSYCL*>(selectedDevices.back());
								addedDevice->setCardID(nCard);
								addedDevice->setStackID(nStack);
								addedDevice->setNumStack(subs.size());
								addedDevice->setSliceID(nSlice++);
								addedDevice->setNumSlice(1);
								if (verbose)
								{
									std::cout << std::string(80, '*') << std::endl;
									std::cout << "Created SYCL device is " << addedDevice->getName() << std::endl;
									std::cout << "maxComputeUnit= " << addedDevice->maxUnit << ", maxWorkGroupSize= " << addedDevice->maxGroup << ", globalMemSize= " << addedDevice->globalMem << std::endl;
									addedDevice->printDeviceInfo();
									std::cout << "\n";
								}
							}
						}
 #endif
						else
#endif
						{
							int numSlice = -1;
#if defined(SYCL_EXT_INTEL_CSLICE) || defined(SYCL_EXT_INTEL_QUEUE_INDEX)
 #ifdef SYCL_EXT_INTEL_CSLICE
							if (isSubSub && isPartitionableByCSlice(sub))
							{
								auto subsubs = sub.create_sub_devices<sycl::info::partition_property::ext_intel_partition_by_cslice>();
								if (subsubs.size() > 1)
									numSlice = subsubs.size();
							}
 #else
							if (isSubSub && numPartitionableByQueueIndex(sub) > 1)
								numSlice = numPartitionableByQueueIndex(sub);
 #endif
#endif
							selectedDevices.push_back(new devSYCL(ctx, sub, nDevice++, isInOrderQueue, isAsyncSubmission));
							auto addedDevice = dynamic_cast<devSYCL*>(selectedDevices.back());
							addedDevice->setCardID(nCard);
							addedDevice->setStackID(nStack);
							addedDevice->setNumStack(subs.size());
							addedDevice->setSliceID(-1);
							addedDevice->setNumSlice(numSlice);
							if (verbose)
							{
								std::cout << std::string(80, '*') << std::endl;
								std::cout << "Created SYCL device is " << addedDevice->getName() << std::endl;
								std::cout << "maxComputeUnit= " << addedDevice->maxUnit << ", maxWorkGroupSize= " << addedDevice->maxGroup << ", globalMemSize= " << addedDevice->globalMem << std::endl;
								addedDevice->printDeviceInfo();
								std::cout << "\n";
							}
						}
						nStack++;
					}
				}
				else
				{
#if defined(SYCL_EXT_INTEL_CSLICE) || defined(SYCL_EXT_INTEL_QUEUE_INDEX)
 #ifdef SYCL_EXT_INTEL_CSLICE
					if (isSubSub && isPartitionableByCSlice(device))
					{
						auto subsubs = device.create_sub_devices<sycl::info::partition_property::ext_intel_partition_by_cslice>();
						auto ctxctx = sycl::context(subsubs);
						int nSlice = 0;
						for (auto &subsub : subsubs)
						{
							selectedDevices.push_back(new devSYCL(ctxctx, subsub, nDevice++, isInOrderQueue, isAsyncSubmission));
							auto addedDevice = dynamic_cast<devSYCL*>(selectedDevices.back());
							addedDevice->setCardID(nCard);
							addedDevice->setStackID(0);
							addedDevice->setNumStack(1);
							addedDevice->setSliceID(nSlice++);
							addedDevice->setNumSlice(1);
							if (verbose)
							{
								std::cout << std::string(80, '*') << std::endl;
								std::cout << "Created SYCL device is " << addedDevice->getName() << std::endl;
								std::cout << "maxComputeUnit= " << addedDevice->maxUnit << ", maxWorkGroupSize= " << addedDevice->maxGroup << ", globalMemSize= " << addedDevice->globalMem << std::endl;
								addedDevice->printDeviceInfo();
								std::cout << "\n";
							}
						}
					}
 #else
					int maxQindex = numPartitionableByQueueIndex(device);
					if (isSubSub && maxQindex > 1)
					{
						sycl::context ctx;
						int nSlice = 0;
						for (int i = 0; i < maxQindex; i++)
						{
							if (i == 0)
								selectedDevices.push_back(new devSYCL(device, i, nDevice++, isInOrderQueue, isAsyncSubmission));
							else
								selectedDevices.push_back(new devSYCL(ctx, device, i, nDevice++, isInOrderQueue, isAsyncSubmission));
							auto addedDevice = dynamic_cast<devSYCL*>(selectedDevices.back());
							if (i == 0)
								ctx = addedDevice->getContext();
							addedDevice->setCardID(nCard);
							addedDevice->setStackID(0);
							addedDevice->setNumStack(1);
							addedDevice->setSliceID(nSlice++);
							addedDevice->setNumSlice(1);
							if (verbose)
							{
								std::cout << std::string(80, '*') << std::endl;
								std::cout << "Created SYCL device is " << addedDevice->getName() << std::endl;
								std::cout << "maxComputeUnit= " << addedDevice->maxUnit << ", maxWorkGroupSize= " << addedDevice->maxGroup << ", globalMemSize= " << addedDevice->globalMem << std::endl;
								addedDevice->printDeviceInfo();
								std::cout << "\n";
							}
						}
					}
 #endif
					else
#endif
					{
						int numSlice = -1;
#if defined(SYCL_EXT_INTEL_CSLICE) || defined(SYCL_EXT_INTEL_QUEUE_INDEX)
 #ifdef SYCL_EXT_INTEL_CSLICE
						if (isSubSub && isPartitionableByCSlice(device))
						{
							auto subsubs = device.create_sub_devices<sycl::info::partition_property::ext_intel_partition_by_cslice>();
							if (subsubs.size() > 1)
								numSlice = subsubs.size();
						}
 #else
						if (isSubSub && numPartitionableByQueueIndex(device) > 1)
							numSlice = numPartitionableByQueueIndex(device);
 #endif
#endif
						selectedDevices.push_back(new devSYCL(device, nDevice++, isInOrderQueue, isAsyncSubmission));
						auto addedDevice = dynamic_cast<devSYCL*>(selectedDevices.back());
						addedDevice->setCardID(nCard);
						addedDevice->setStackID(0);
						addedDevice->setNumStack(1);
						addedDevice->setSliceID(-1);
						addedDevice->setNumSlice(numSlice);
						if (verbose)
						{
							std::cout << std::string(80, '*') << std::endl;
							std::cout << "Created SYCL device is " << addedDevice->getName() << std::endl;
							std::cout << "maxComputeUnit= " << addedDevice->maxUnit << ", maxWorkGroupSize= " << addedDevice->maxGroup << ", globalMemSize= " << addedDevice->globalMem << std::endl;
							addedDevice->printDeviceInfo();
							std::cout << "\n";
						}
					}
				}
				nCard++;
			}
		}
#ifdef ACC_DOUBLE_PRECISION
		bool isFP64 = true;
		for (auto &select : selectedDevices)
		{
			if (! dynamic_cast<devSYCL*>(select)->canSupportFP64())
			{
				isFP64 = false;
				break;
			}
		}
		if (! isFP64)
		{
			for (auto select : selectedDevices)
				delete select;

			selectedDevices.clear();
			std::cerr << "Double-Precision for Accelerator is requested but there is device which cannot support FP64.\n";
		}
#endif
	}
	else if (select == syclDeviceType::cpu)
	{
		int nDevice = 0;
		auto devices = sycl::device::get_devices(sycl::info::device_type::cpu);
		for (auto &device : devices)
		{
/*
// TODO: For heterogeneous run, is this necessary?
			auto ndev = device.get_info<sycl::info::device::partition_max_sub_devices>();
			if (ndev > 1)
			{
				auto subs = device.create_sub_devices<sycl::info::partition_property::partition_by_affinity_domain>(sycl::info::partition_affinity_domain::numa);
				for (auto &sub : subs)
					selectedDevices.push_back(new devSYCL(sub, nDevice++, isInOrderQueue, isAsyncSubmission));
			}
			else
*/
				selectedDevices.push_back(new devSYCL(device, nDevice++, isInOrderQueue, isAsyncSubmission));
				auto addedDevice = dynamic_cast<devSYCL*>(selectedDevices.back());
				addedDevice->setCardID(-1);
				addedDevice->setStackID(-1);
				addedDevice->setNumStack(-1);
				addedDevice->setSliceID(-1);
				addedDevice->setNumSlice(-1);
				if (verbose)
				{
					std::cout << std::string(80, '*') << std::endl;
					std::cout << "Created SYCL device is " << addedDevice->getName() << std::endl;
					std::cout << "maxComputeUnit= " << addedDevice->maxUnit << ", maxWorkGroupSize= " << addedDevice->maxGroup << ", globalMemSize= " << addedDevice->globalMem << std::endl;
					addedDevice->printDeviceInfo();
					std::cout << "\n";
				}
		}
	}
	else
		std::cerr << "Only GPU and CPU devices are supported.\n";

	return selectedDevices;
}

void MlOptimiserSYCL::checkDevices()
{
	std::vector<sycl::device> devices = sycl::device::get_devices();
	if (devices.size() == 0)
		REPORT_ERROR("NO SYCL devices are found");

	for (auto &device : devices)
	{
		std::cout << "\nPlatform: " << device.get_platform().get_info<sycl::info::platform::name>() << std::endl;

		PRINT_SYCL_INFO(name);
/*
		if (! device.is_host()) PRINT_SYCL_INFO(max_clock_frequency);
		PRINT_SYCL_INFO(version);
		PRINT_SYCL_INFO(driver_version);
		PRINT_SYCL_INFO(partition_max_sub_devices);
		PRINT_SYCL_INFO(profiling_timer_resolution);
		PRINT_SYCL_INFO(queue_profiling);

		PRINT_SYCL_INFO(max_compute_units);
		PRINT_SYCL_INFO(max_work_group_size);
		PRINT_SYCL_INFO(max_work_item_dimensions);
		const auto isizes = device.get_info<sycl::info::device::max_work_item_sizes>();
		std::cout << "  max_work_item_sizes: " << isizes[0] << " x " << isizes[1] << " x " << isizes[2] << std::endl;

		PRINT_SYCL_INFO(preferred_vector_width_int);
		PRINT_SYCL_INFO(preferred_vector_width_long);
		PRINT_SYCL_INFO(preferred_vector_width_half);
		PRINT_SYCL_INFO(preferred_vector_width_float);
		PRINT_SYCL_INFO(preferred_vector_width_double);
		PRINT_SYCL_INFO(native_vector_width_int);
		PRINT_SYCL_INFO(native_vector_width_long);
		PRINT_SYCL_INFO(native_vector_width_half);
		PRINT_SYCL_INFO(native_vector_width_float);
		PRINT_SYCL_INFO(native_vector_width_double);

		PRINT_SYCL_INFO(max_mem_alloc_size);
		PRINT_SYCL_INFO(global_mem_cache_line_size);
		PRINT_SYCL_INFO(global_mem_cache_size);
		PRINT_SYCL_INFO(global_mem_size);
		PRINT_SYCL_INFO(local_mem_size);

		PRINT_SYCL_INFO(image_support);
		PRINT_SYCL_INFO(max_read_image_args);
		PRINT_SYCL_INFO(max_write_image_args);
		PRINT_SYCL_INFO(image2d_max_height);
		PRINT_SYCL_INFO(image2d_max_width);
		PRINT_SYCL_INFO(image3d_max_height);
		PRINT_SYCL_INFO(image3d_max_width);
		PRINT_SYCL_INFO(image3d_max_depth);
		PRINT_SYCL_INFO(image_max_buffer_size);
		PRINT_SYCL_INFO(image_max_array_size);
		PRINT_SYCL_INFO(max_samplers);
		PRINT_SYCL_INFO(max_parameter_size);

		const auto domains = device.get_info<sycl::info::device::partition_affinity_domains>();
		std::cout << "  partition_affinity_domain:";
		for (auto &domain : domains)
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

		const auto hconfigs = device.get_info<sycl::info::device::half_fp_config>();
		std::cout << "  half_fp_config:";
		for (auto &hconfig : hconfigs)
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

		const auto sconfigs = device.get_info<sycl::info::device::single_fp_config>();
		std::cout << "  single_fp_config:";
		for (auto &sconfig : sconfigs)
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

		const auto dconfigs = device.get_info<sycl::info::device::double_fp_config>();
		std::cout << "  double_fp_config:";
		for (auto &dconfig : dconfigs)
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
*/
		std::cout << std::endl;
	}
}

void MlOptimiserSYCL::setupDevice()
{
#ifdef USE_EXISTING_SYCL_DEVICE
	_devAcc = new devSYCL(dynamic_cast<devSYCL*>(bundle->getSyclDevice()));
#else
// This will create separate device queue while the above is using existing queue
	auto dev = dynamic_cast<devSYCL*>(bundle->getSyclDevice());
	auto qType = dev->getSyclQueueType();
	auto isAsync = dev->isAsyncQueue();
	auto computeIndex = dev->getComputeIndex();
	auto q = dev->getQueue();
	auto c = q->get_context();
	auto d = q->get_device();

 #ifdef SYCL_EXT_INTEL_QUEUE_INDEX
	if (computeIndex >= 0)
		_devAcc = new devSYCL(c, d, computeIndex, qType, dev->getDeviceID(), isAsync);
	else
 #endif
		_devAcc = new devSYCL(c, d, qType, dev->getDeviceID(), isAsync);

	_devAcc->setCardID(dev->getCardID());
	_devAcc->setStackID(dev->getStackID());
	_devAcc->setNumStack(dev->getNumStack());
	_devAcc->setSliceID(dev->getSliceID());
	_devAcc->setNumSlice(dev->getNumSlice());

	auto ndev = dynamic_cast<devSYCL*>(_devAcc);
	ndev->setSyclHostPool(dev->getSyclHostPool(), dev->getSyclHostBlockSize());
	ndev->setSyclDevicePool(dev->getSyclDevicePool(), dev->getSyclDeviceBlockSize());
#endif
}

void MlOptimiserSYCL::resetData()
{
	transformer1.clear();
	transformer2.clear();

	for (int i = 0; i < baseMLO->mymodel.nr_classes; i++)
	{
		if (_useStream)
		{
			auto dev = dynamic_cast<devSYCL*>(_devAcc);
			auto isAsync = dev->isAsyncQueue();
			auto computeIndex = dev->getComputeIndex();
			auto q = dev->getQueue();
			auto c = q->get_context();
			auto d = q->get_device();

 #ifdef SYCL_EXT_INTEL_QUEUE_INDEX
			if (computeIndex >= 0)
				classStreams.push_back(new devSYCL(c, d, computeIndex, i, true, isAsync));
 #endif
			else
				classStreams.push_back(new devSYCL(c, d, i, true, isAsync));
	
			classStreams.back()->setCardID(dynamic_cast<devSYCL*>(_devAcc)->getCardID());
			classStreams.back()->setStackID(dynamic_cast<devSYCL*>(_devAcc)->getStackID());
			classStreams.back()->setNumStack(dynamic_cast<devSYCL*>(_devAcc)->getNumStack());
			classStreams.back()->setSliceID(dynamic_cast<devSYCL*>(_devAcc)->getSliceID());
			classStreams.back()->setNumSlice(dynamic_cast<devSYCL*>(_devAcc)->getNumSlice());

			auto ndev = dynamic_cast<devSYCL*>(classStreams.back());
			ndev->setSyclHostPool(dev->getSyclHostPool(), dev->getSyclHostBlockSize());
			ndev->setSyclDevicePool(dev->getSyclDevicePool(), dev->getSyclDeviceBlockSize());
		}
		else
			classStreams.push_back(_devAcc);
	}
}

void MlOptimiserSYCL::expectationOneParticle(unsigned long my_part_id, const int thread_id)
{
	AccPtrFactory ptrFactory(AccType::accCPU);
	accDoExpectationOneParticle<MlOptimiserSYCL>(this, my_part_id, thread_id, ptrFactory);
}

void MlOptimiserSYCL::doThreadExpectationSomeParticles(const int thread_id)
{
	size_t first_ipart = 0, last_ipart = 0;
	while (baseMLO->exp_ipart_ThreadTaskDistributor->getTasks(first_ipart, last_ipart))
	{
		for (long unsigned ipart = first_ipart; ipart <= last_ipart; ipart++)
		{
			expectationOneParticle(baseMLO->exp_my_first_part_id + ipart, thread_id);
		}
	}
}
#endif // _SYCL_ENABLED

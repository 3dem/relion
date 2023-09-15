#ifdef _SYCL_ENABLED
// A large amount of this code is direct from cuda_ml_optimizer and so could
// be shared (but possibly with difficulty since it is enough different that
// we either need a lot of #ifdefs, or a lot of macros/some other mechanism to
// abstract the differences).  The biggest differences are the type of memory
// objects used (std::vector vs. CudaGlobalPtr and CudaCustomAllocator), the
// lack of transfers to/from the device, and on-device operations (which are
// replaced by loops/function calls).
//
// CudaFFT has been replaced with lib FFTW, if RELION is configured with mix
// precision, both single and double precision FFTW are linked into RELION.
// Install fftw-static.x86_64 and fftw-static.i686 to get the libraries without
// having to pull them at build time.  Over time we hope to replace FFTW with
// MKL.
//
// NOTE:  Since the GPU code was ported back to CPU there may be additional
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

#include "src/acc/sycl/sycl_dev.h"

#define PRINT_SYCL_INFO(x) std::cout << "  "#x": " << device.get_info<sycl::info::device::x>() << std::endl

template <typename RangeTy, typename ElemTy>
bool contains(RangeTy &&Range, const ElemTy &Elem) {
    return std::find(Range.begin(), Range.end(), Elem) != Range.end();
}

bool isPartitionableByCSlice(sycl::device &Dev) {
    return contains(Dev.get_info<sycl::info::device::partition_properties>(), sycl::info::partition_property::ext_intel_partition_by_cslice);
}

std::mutex	fft_mutex;

MlSyclDataBundle::MlSyclDataBundle(virtualSYCL *dev) : generateProjectionPlanOnTheFly {false}, _devAcc {dev}
{
}

MlSyclDataBundle::~MlSyclDataBundle()
{
	projectors.clear();
	backprojectors.clear();
	coarseProjectionPlans.clear();
//	delete _devAcc;
}

void MlSyclDataBundle::setup(MlOptimiser *baseMLO)
{
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

MlOptimiserSYCL::MlOptimiserSYCL(MlOptimiser *baseMLOptimiser, MlSyclDataBundle *b, const char *timing_fnm) :
	baseMLO {baseMLOptimiser},
	transformer1 {baseMLOptimiser->mymodel.data_dim},
	transformer2 {baseMLOptimiser->mymodel.data_dim},
	refIs3D {baseMLO->mymodel.ref_dim == 3},
	dataIs3D {baseMLO->mymodel.data_dim == 3},
	shiftsIs3D {baseMLO->mymodel.data_dim == 3 || baseMLO->mydata.is_tomo},
	generateProjectionPlanOnTheFly {bundle->generateProjectionPlanOnTheFly},
	bundle {b}
{
	setupDevice();
}

MlOptimiserSYCL::~MlOptimiserSYCL()
{
#ifdef USE_SYCL_STREAM
	for (auto q : classStreams)
		delete q;
#endif

	classStreams.clear();
#ifndef USE_EXISTING_SYCL_DEVICE
	delete _devAcc;
#endif
}

std::vector<virtualSYCL*> MlOptimiserSYCL::getDevices(const syclDeviceType select, const syclBackendType BE, const bool verbose)
{
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
						if (isPartitionableByCSlice(sub))
						{
							auto subsubs = sub.create_sub_devices<sycl::info::partition_property::ext_intel_partition_by_cslice>();
#ifdef USE_SUBSUB_DEVICE
							auto ctxctx = sycl::context(subsubs);
							int nSlice = 0;
							for (auto &subsub : subsubs)
							{
								selectedDevices.push_back(new devSYCL(ctxctx, subsub, nDevice++));
								auto addedDevice = dynamic_cast<devSYCL*>(selectedDevices.back());
								addedDevice->setCardID(nCard);
								addedDevice->setStackID(nStack);
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
#else
							selectedDevices.push_back(new devSYCL(ctx, sub, nDevice++));
							auto addedDevice = dynamic_cast<devSYCL*>(selectedDevices.back());
							addedDevice->setCardID(nCard);
							addedDevice->setStackID(nStack);
							addedDevice->setSliceID(-1);
							addedDevice->setNumSlice(subsubs.size());
							if (verbose)
							{
								std::cout << std::string(80, '*') << std::endl;
								std::cout << "Created SYCL device is " << addedDevice->getName() << std::endl;
								std::cout << "maxComputeUnit= " << addedDevice->maxUnit << ", maxWorkGroupSize= " << addedDevice->maxGroup << ", globalMemSize= " << addedDevice->globalMem << std::endl;
								addedDevice->printDeviceInfo();
								std::cout << "  Number of Slice: " << subsubs.size() << "\n\n";
							}
#endif
						}
						else
						{
							selectedDevices.push_back(new devSYCL(ctx, sub, nDevice++));
							auto addedDevice = dynamic_cast<devSYCL*>(selectedDevices.back());
							addedDevice->setCardID(nCard);
							addedDevice->setStackID(nStack);
							addedDevice->setSliceID(-1);
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
						nStack++;
					}
				}
				else
				{
					if (isPartitionableByCSlice(device))
					{
						auto subsubs = device.create_sub_devices<sycl::info::partition_property::ext_intel_partition_by_cslice>();
#ifdef USE_SUBSUB_DEVICE
						auto ctxctx = sycl::context(subsubs);
						int nSlice = 0;
						for (auto &subsub : subsubs)
						{
							selectedDevices.push_back(new devSYCL(ctxctx, subsub, nDevice++));
							auto addedDevice = dynamic_cast<devSYCL*>(selectedDevices.back());
							addedDevice->setCardID(nCard);
							addedDevice->setStackID(0);
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
#else
						selectedDevices.push_back(new devSYCL(device, nDevice++));
						auto addedDevice = dynamic_cast<devSYCL*>(selectedDevices.back());
						addedDevice->setCardID(nCard);
						addedDevice->setStackID(0);
						addedDevice->setSliceID(-1);
						addedDevice->setNumSlice(subsubs.size());
						if (verbose)
						{
							std::cout << std::string(80, '*') << std::endl;
							std::cout << "Created SYCL device is " << addedDevice->getName() << std::endl;
							std::cout << "maxComputeUnit= " << addedDevice->maxUnit << ", maxWorkGroupSize= " << addedDevice->maxGroup << ", globalMemSize= " << addedDevice->globalMem << std::endl;
							addedDevice->printDeviceInfo();
							std::cout << "  Number of Slice: " << subsubs.size() << "\n\n";
						}
#endif
					}
					else
					{
						selectedDevices.push_back(new devSYCL(device, nDevice++));
						auto addedDevice = dynamic_cast<devSYCL*>(selectedDevices.back());
						addedDevice->setCardID(nCard);
						addedDevice->setStackID(0);
						addedDevice->setSliceID(-1);
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
				delete dynamic_cast<devSYCL*>(select);

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
					selectedDevices.push_back(new devSYCL(sub, nDevice++));
			}
			else
*/
				selectedDevices.push_back(new devSYCL(device, nDevice++));
				auto addedDevice = dynamic_cast<devSYCL*>(selectedDevices.back());
				addedDevice->setCardID(-1);
				addedDevice->setStackID(-1);
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

	for (auto& device : devices)
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

		const auto hconfigs = device.get_info<sycl::info::device::half_fp_config>();
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

		const auto sconfigs = device.get_info<sycl::info::device::single_fp_config>();
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

		const auto dconfigs = device.get_info<sycl::info::device::double_fp_config>();
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
	auto q = dynamic_cast<devSYCL*>(bundle->getSyclDevice())->getQueue();
	auto c = q->get_context();
	auto d = q->get_device();
	_devAcc = new devSYCL(c, d, bundle->getSyclDevice()->getDeviceID());
	_devAcc->setCardID(bundle->getSyclDevice()->getCardID());
	_devAcc->setStackID(bundle->getSyclDevice()->getStackID());
	_devAcc->setSliceID(bundle->getSyclDevice()->getSliceID());
	_devAcc->setNumSlice(bundle->getSyclDevice()->getNumSlice());
//	_devAcc->printDeviceInfo();
#endif
}

void MlOptimiserSYCL::resetData()
{
	transformer1.clear();
	transformer2.clear();

	for (int i = 0; i < baseMLO->mymodel.nr_classes; i++)
	{
#ifdef USE_SYCL_STREAM
		auto d = dynamic_cast<devSYCL*>(_devAcc)->getDevice();
		auto c = dynamic_cast<devSYCL*>(_devAcc)->getContext();
		classStreams.push_back(new devSYCL(c, d, syclQueueType::inOrder, i));
		classStreams.back()->setCardID(dynamic_cast<devSYCL*>(_devAcc)->getCardID());
		classStreams.back()->setStackID(dynamic_cast<devSYCL*>(_devAcc)->getStackID());
		classStreams.back()->setSliceID(dynamic_cast<devSYCL*>(_devAcc)->getSliceID());
		classStreams.back()->setNumSlice(dynamic_cast<devSYCL*>(_devAcc)->getNumSlice());
#else
		classStreams.push_back(_devAcc);
#endif
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

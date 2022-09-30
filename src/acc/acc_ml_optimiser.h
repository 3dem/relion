#ifndef ACC_ML_OPTIMISER_H_
#define ACC_ML_OPTIMISER_H_

#include "src/acc/acc_ptr.h"

#ifdef ALTCPU
#include <tbb/spin_mutex.h>
#endif

/*
#ifdef ACC_DOUBLE_PRECISION
#define XFLOAT double
#else
#define XFLOAT float
#endif
*/

class SamplingParameters
{
public:
	unsigned long nr_dir,
	nr_psi,
	nr_trans,
	nr_oversampled_rot,
	nr_oversampled_trans,
	nr_images,
	current_oversampling,
	current_image_size,
	iclass_min, iclass_max,
	idir_min, idir_max,
	ipsi_min, ipsi_max,
	itrans_min, itrans_max;
	std::string current_img;

	SamplingParameters():
		nr_dir(0),
		nr_psi(0),
		nr_trans(0),
		nr_oversampled_rot(0),
		nr_oversampled_trans(0),
		nr_images(0),
		current_oversampling(0),
		current_image_size(0),
		iclass_min(0), iclass_max(0),
		idir_min(0), idir_max(0),
		ipsi_min(0), ipsi_max(0),
		itrans_min(0), itrans_max(0),
		current_img()
	{};
};

class Indices
{
public:
	size_t fineIdx,
	coarseIdx,
	iclass,
	idir,
	ipsi,
	itrans,
	ioverrot,
	iovertrans;

	Indices():
		fineIdx(0),
		coarseIdx(0),
		iclass(0),
		idir(0),
		ipsi(0),
		itrans(0),
		ioverrot(0),
		iovertrans(0)
	{};

	void fineIndexToFineIndices(SamplingParameters sp) // converts an "ihidden_over" (finely sampled) index to partial indices (and coarse index)
	{
		int oversamples = sp.nr_oversampled_rot*sp.nr_oversampled_trans;
		size_t t_idx = fineIdx;
		iclass = floor( t_idx / ( sp.nr_dir * sp.nr_psi * sp.nr_trans * oversamples ));
		t_idx   -= iclass     * ( sp.nr_dir * sp.nr_psi * sp.nr_trans * oversamples );
		idir   = floor( t_idx / ( sp.nr_psi * sp.nr_trans * oversamples ));
		t_idx   -= idir       * ( sp.nr_psi * sp.nr_trans * oversamples );
		ipsi   = floor( t_idx / ( sp.nr_trans * oversamples ));
		t_idx   -= ipsi       * ( sp.nr_trans * oversamples );
		itrans = floor( t_idx /  oversamples );
		t_idx   -= itrans     *  oversamples ;
		ioverrot = floor( t_idx / sp.nr_oversampled_trans );
		t_idx   -= ioverrot  *   sp.nr_oversampled_trans ;
		iovertrans = t_idx ;

		coarseIdx = sp.nr_trans * sp.nr_psi * idir   +   sp.nr_trans * ipsi   +   itrans;
	}

	void fineIndicesToFineIndex(SamplingParameters sp) // converts partial indices to an "ihidden_over" (finely sampled) index // FIXME Untested
	{
		int oversamples = sp.nr_oversampled_rot*sp.nr_oversampled_trans;
		size_t idx = 0;
		idx += iclass   * sp.nr_dir * sp.nr_psi * sp.nr_trans * oversamples;
		idx += idir     * sp.nr_psi * sp.nr_trans * oversamples;
		idx += ipsi     * sp.nr_trans * oversamples;
		idx += itrans   * oversamples;
		idx += ioverrot * sp.nr_oversampled_trans;
		idx += iovertrans;
		fineIdx = idx;
	}

	void coarseIndexToCoarseIndices(SamplingParameters sp) // converts an "ihidden" (coarsely sampled) index to coarse partial indices // FIXME Untested
	{
		size_t t_idx = coarseIdx;
		iclass = floor( t_idx / ( sp.nr_dir * sp.nr_psi * sp.nr_trans));
		t_idx   -= iclass     * ( sp.nr_dir * sp.nr_psi * sp.nr_trans);
		idir   = floor( t_idx / ( sp.nr_psi * sp.nr_trans ));
		t_idx   -= idir       * ( sp.nr_psi * sp.nr_trans  );
		ipsi   = floor( t_idx / ( sp.nr_trans ));
		t_idx   -= ipsi       * ( sp.nr_trans  );
		itrans = t_idx ;
		ioverrot   = 0;
		iovertrans = 0;
	}

	void coarseIndicesToCoarseIndex(SamplingParameters sp) // converts coarse partial indices to an "ihidden" (coarsely sampled) index // FIXME Untested
	{
		size_t idx = 0;
		idx += idir     * sp.nr_psi * sp.nr_trans;
		idx += ipsi     * sp.nr_trans;
		idx += itrans;
		coarseIdx = idx;
	}
};


class OptimisationParamters
{
public:
	unsigned metadata_offset;

	unsigned long part_id;

	std::vector<MultidimArray<Complex > > Fimg, Fimg_nomask, local_Fimgs_shifted, local_Fimgs_shifted_nomask;
	std::vector<MultidimArray<RFLOAT> > Fctf, local_Fctf, local_Minvsigma2, FstMulti;
	std::vector<int> pointer_dir_nonzeroprior, pointer_psi_nonzeroprior;
	std::vector<RFLOAT> directions_prior, psi_prior, local_sqrtXi2;
	std::vector<RFLOAT> highres_Xi2_img, min_diff2;
	MultidimArray<bool> Mcoarse_significant;
	// And from storeWeightedSums
	std::vector<RFLOAT> sum_weight, significant_weight, max_weight;
	std::vector< std::vector<RFLOAT> > sum_weight_class;
	std::vector<Matrix1D<RFLOAT> > old_offset, prior;
	std::vector<MultidimArray<RFLOAT> > power_img;
	MultidimArray<XFLOAT> Mweight;
	std::vector<Indices> max_index;

	OptimisationParamters (unsigned nr_images, unsigned long part_id):
		metadata_offset(0),
		part_id(part_id)
	{
		power_img.resize(nr_images);
		highres_Xi2_img.resize(nr_images);
		Fimg.resize(nr_images);
		Fimg_nomask.resize(nr_images);
		Fctf.resize(nr_images);
		old_offset.resize(nr_images);
		prior.resize(nr_images);
		max_index.resize(nr_images);
		sum_weight_class.resize(nr_images);
	};
};

class ProjectionParams
{

public:
	std::vector< size_t > orientation_num; 					// the number of significant orientation for each class
	size_t orientationNumAllClasses;							// sum of the above
	std::vector< RFLOAT > rots, tilts, psis;
	std::vector< size_t > iorientclasses, iover_rots;

	// These are arrays which detial the number of entries in each class, and where each class starts.
	// NOTE: There is no information about which class each class_idx refers to, there is only
	// a distinction between different classes.
	std::vector< size_t > class_entries, class_idx;
	inline
	ProjectionParams():

		rots(),
		tilts(),
		psis(),
		iorientclasses(),
		iover_rots(),

		class_entries(),
		class_idx(),
		orientation_num(),
		orientationNumAllClasses(0)

	{};

	inline
	ProjectionParams(size_t classes):

		rots(),
		tilts(),
		psis(),
		iorientclasses(),
		iover_rots(),

		class_entries(classes),
		class_idx(classes),
		orientation_num(classes),
		orientationNumAllClasses(0)
	{
		class_idx[0]=0;
		class_entries[0]=0;
	};


	// constructor that slices out a part of a parent ProjectionParams, assumed to contain a single (partial or entire) class
	inline
	ProjectionParams(ProjectionParams &parent, size_t start, size_t end):
		rots(				parent.rots.begin() 			+start,  	parent.rots.begin() 			+end),
		tilts(				parent.tilts.begin() 			+start, 	parent.tilts.begin() 			+end),
		psis(				parent.psis.begin() 			+start,  	parent.psis.begin() 			+end),
		iorientclasses( 	parent.iorientclasses.begin() 	+start,  	parent.iorientclasses.begin() 	+end),
		iover_rots(			parent.iover_rots.begin() 		+start,  	parent.iover_rots.begin() 		+end),
		orientation_num(1),
		orientationNumAllClasses(0),
		class_entries(1,end-start),
		class_idx(1,0) // NOTE: this is NOT the class, but rather where in these partial PrjParams to start, which is @ 0.
	{};

public:
	// Appends new values into the projection parameters for later use.
	// class_idx is used as such:
	// the n:th class (beginning with 0:th)
	// begins @ element class_idx[n]
	// ends   @ element class_idx[n]+class_entries[n]

	void pushBackAll(size_t iclass, RFLOAT NEWrot,RFLOAT NEWtilt ,RFLOAT NEWpsi, size_t NEWiorientclasses,size_t NEWiover_rots)
	{
		// incremement the counter for this class
		class_entries[iclass]++;
		// and push a new entry
		rots.push_back(NEWrot);
		tilts.push_back(NEWtilt);
		psis.push_back(NEWpsi);
		iorientclasses.push_back(NEWiorientclasses);
		iover_rots.push_back(NEWiover_rots);
	}
};

class IndexedDataArrayMask
{
public:
	// indexes of job partition
	//   every element in jobOrigin    is a reference to point to a position in a IndexedDataArray.weights array where that job starts RELATIVE to firstPos
	//   every element in jobExtent    specifies the number of weights for that job
	AccPtr<size_t> jobOrigin, jobExtent;

	size_t firstPos, lastPos; // positions in indexedDataArray data and index arrays to slice out
	size_t weightNum, jobNum; // number of weights and jobs this class
	

public:
		
	IndexedDataArrayMask(AccPtrFactory ptrFactory) :
		firstPos(), lastPos(), weightNum(), jobNum()
	{
		jobOrigin = ptrFactory.make<size_t>();
		jobExtent = ptrFactory.make<size_t>();
	}
	
	void setNumberOfJobs(size_t newSize)
	{
		jobNum=newSize;
		jobOrigin.setSize(newSize);
		jobExtent.setSize(newSize);
	}

	void setNumberOfWeights(size_t newSize)
	{
		weightNum=newSize;
	}

	inline
	 ~IndexedDataArrayMask()
	{
//		jobOrigin.free_host();
//		jobExtent.free_host();
	};
};

class IndexedDataArray
{
public:
	//actual data
	AccPtr<XFLOAT> weights;

	// indexes with same length as data
	// -- basic indices ---------------------------------
	//     rot_id  = id of rot     = which of all POSSIBLE orientations                               this weight signifies
	//     rot_idx = index of rot  = which in the sequence of the determined significant orientations this weight signifies
	//   trans_id  = id of trans   = which of all POSSIBLE translations                               this weight signifies
	// -- special indices ---------------------------------
	//   ihidden_overs  =  mapping to MWeight-based indexing for compatibility
	AccPtr<size_t> rot_id, rot_idx, trans_idx, ihidden_overs;

public:
	
	inline
	 IndexedDataArray(AccPtrFactory ptrFactory):
	 	weights(ptrFactory.make<XFLOAT>()),
		rot_id(ptrFactory.make<size_t>()),
		rot_idx(ptrFactory.make<size_t>()),
		trans_idx(ptrFactory.make<size_t>()),
		ihidden_overs(ptrFactory.make<size_t>())
	{};

	inline
	 IndexedDataArray(IndexedDataArray &parent, IndexedDataArrayMask &mask):
		weights(      parent.weights, mask.firstPos, mask.weightNum),
		rot_id(       parent.rot_id, mask.firstPos, mask.weightNum),
		rot_idx(      parent.rot_idx, mask.firstPos, mask.weightNum),
		trans_idx(    parent.trans_idx, mask.firstPos, mask.weightNum),
		ihidden_overs(parent.ihidden_overs, mask.firstPos, mask.weightNum)
	{};

public:

	void setDataSize(size_t newSize)
	{
		weights.setSize(newSize);
		rot_id.setSize(newSize);
		rot_idx.setSize(newSize);
		trans_idx.setSize(newSize);
		ihidden_overs.setSize(newSize);
	}

	void host_alloc_all()
	{
		weights.freeHostIfSet();
		weights.hostAlloc();
		rot_id.freeHostIfSet();
		rot_id.hostAlloc();
		rot_idx.freeHostIfSet();
		rot_idx.hostAlloc();
		trans_idx.freeHostIfSet();
		trans_idx.hostAlloc();
		ihidden_overs.freeHostIfSet();
		ihidden_overs.hostAlloc();
	}

	void device_alloc_all()
	{
		weights.freeDeviceIfSet();
		weights.deviceAlloc();
		rot_id.freeDeviceIfSet();
		rot_id.deviceAlloc();
		rot_idx.freeDeviceIfSet();
		rot_idx.deviceAlloc();
		trans_idx.freeDeviceIfSet();
		trans_idx.deviceAlloc();
		ihidden_overs.freeDeviceIfSet();
		ihidden_overs.deviceAlloc();
	}

	void dual_alloc_all()
	{
		host_alloc_all();
		device_alloc_all();
	}
	
	void dual_free_all()
	{
		weights.freeDeviceIfSet();
		rot_id.freeDeviceIfSet();
		rot_idx.freeDeviceIfSet();
		trans_idx.freeDeviceIfSet();
		ihidden_overs.freeDeviceIfSet();
		weights.freeHostIfSet();
		rot_id.freeHostIfSet();
		rot_idx.freeHostIfSet();
		trans_idx.freeHostIfSet();
		ihidden_overs.freeHostIfSet();
	}
	
	~IndexedDataArray()
	{
		dual_free_all();
	}


};

#endif

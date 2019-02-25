#ifndef RELION_H_
#define RELION_H_

//This is a main header - it includes everything else.


#include <external/Healpix_2.15a/arr.h>
#include <external/Healpix_2.15a/cxxutils.h>
#include <external/Healpix_2.15a/datatypes.h>
#include <external/Healpix_2.15a/geom_utils.h>
#include <external/Healpix_2.15a/healpix_base.h>
#include <external/Healpix_2.15a/lsconstants.h>
#include <external/Healpix_2.15a/message_error.h>
#include <external/Healpix_2.15a/openmp_support.h>
#include <external/Healpix_2.15a/pointing.h>
#include <external/Healpix_2.15a/vec3.h>
#include <src/args.h>
#include <src/backprojector.h>
#include <src/blobs.h>
#include <src/ctf.h>
#include <src/error.h>
#include <src/euler.h>
#include <src/exp_model.h>
#include <src/fftw.h>
#include <src/filename.h>
#include <src/funcs.h>
#include <src/gcc_version.h>
#include <src/healpix_map.h>
#include <src/healpix_sampling.h>
#include <src/image.h>
#include <src/macros.h>
#include <src/matrix1d.h>
#include <src/matrix2d.h>
#include <src/memory.h>
#include <src/metadata_container.h>
#include <src/metadata_label.h>
#include <src/metadata_table.h>
#include <src/ml_model.h>
#include <src/ml_optimiser.h>
#include <src/ml_optimiser_mpi.h>
#include <src/mpi.h>
#include <src/multidim_array.h>
#include <src/numerical_recipes.h>
#include <src/parallel.h>
#include <src/projector.h>
#include <src/rwIMAGIC.h>
#include <src/rwMRC.h>
#include <src/rwSPIDER.h>
#include <src/strings.h>
#include <src/symmetries.h>
#include <src/timer.h>
#include <src/transformations.h>
#include <src/util.h>

#endif // RELION_H_

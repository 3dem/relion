## LibRelion:

## Place generated object files (.o) into the same directory as their source
## files, in order to avoid collisions when non-recursive make is used.
AUTOMAKE_OPTIONS = subdir-objects

## Additional flags to pass to aclocal when it is invoked automatically at
## make time. The ${ACLOCAL_FLAGS} variable is picked up from the environment
## to provide a way for the user to supply additional arguments.
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

## Define an independent executable script for inclusion in the distribution
## archive.  However, it will not be installed on an end user's system due to
## the noinst_ prefix.
dist_noinst_SCRIPTS = autogen.sh

# Install my_script in $(bindir) and distribute it.
dist_bin_SCRIPTS = \
scripts/star_datablock_ctfdat \
scripts/star_datablock_stack \
scripts/star_datablock_singlefiles \
scripts/star_loopheader \
scripts/star_plottable \
scripts/star_printtable \
scripts/qsub.csh

### Include tar.gz of FFTW and FLTK in distribution, as well as all source files for the GUI
EXTRA_DIST = INSTALL.sh \
external/fftw-3.2.2.tar.gz \
external/fltk-1.3.0.tar.gz 

###############################################################################
#
# Core library
#
###############################################################################

## Define a libtool archive target "librelion-@RELION_API_VERSION@.la", with
## @RELION_API_VERSION@ substituted into the generated Makefile at configure
## time.
## The libtool archive file (.la) will be installed into the directory named
## by the predefined variable $(bindir), along with the actual shared library
## file (.so).
lib_LTLIBRARIES = librelion-@RELION_API_VERSION@.la

## Define the source file list for the "librelion-@RELION_API_VERSION@.la"
## target.  Note that @RELION_API_VERSION@ is not interpreted by Automake and
## will therefore be treated as if it were literally part of the target name,
## and the variable name derived from that.
## The file extension .cc is recognized by Automake, and makes it produce
## rules which invoke the C++ compiler to produce a libtool object file (.lo)
## from each source file.  Note that it is not necessary to list header files
## which are already listed elsewhere in a _HEADERS variable assignment.
librelion_@RELION_API_VERSION@_la_SOURCES = src/args.cpp \
src/assembly.cpp \
src/autopicker.cpp \
src/backprojector.cpp \
src/complex.cpp \
src/ctf.cpp \
src/ctffind_runner.cpp \
src/error.cpp \
src/euler.cpp \
src/exp_model.cpp \
src/fftw.cpp \
src/filename.cpp \
src/funcs.cpp \
src/healpix_sampling.cpp \
src/image.cpp \
src/mask.cpp \
src/matrix1d.cpp \
src/matrix2d.cpp \
src/memory.cpp \
src/metadata_container.cpp \
src/metadata_label.cpp \
src/metadata_table.cpp \
src/ml_model.cpp \
src/ml_optimiser.cpp \
src/multidim_array.cpp \
src/numerical_recipes.cpp \
src/parallel.cpp \
src/particle_polisher.cpp \
src/particle_sorter.cpp \
src/postprocessing.cpp \
src/preprocessing.cpp \
src/projector.cpp \
src/strings.cpp \
src/symmetries.cpp \
src/tabfuncs.cpp \
src/time.cpp \
src/transformations.cpp \
src/Healpix_2.15a/cxxutils.cc \
src/Healpix_2.15a/healpix_base.cc 

## Define the list of public header files and their install location.  The
## nobase_ prefix instructs Automake to not strip the directory part from each
## filename, in order to avoid the need to define separate file lists for each
## installation directory.  This only works if the directory hierarchy in the
## source tree matches the hierarchy at the install location, however.
relion_includedir = $(includedir)/relion-$(RELION_API_VERSION)
nobase_relion_include_HEADERS = relion.h \
src/args.h \
src/assembly.h \
src/autopicker.h \
src/backprojector.h \
src/complex.h \
src/ctf.h \
src/ctffind_runner.h \
src/displayer.h \
src/error.h \
src/euler.h \
src/exp_model.h \
src/fftw.h \
src/filename.h \
src/funcs.h \
src/gcc_version.h \
src/gui_entries.h \
src/gui_jobwindow.h \
src/gui_mainwindow.h \
src/healpix_sampling.h \
src/image.h \
src/macros.h \
src/mask.h \
src/matrix1d.h \
src/matrix2d.h \
src/memory.h \
src/metadata_container.h \
src/metadata_label.h \
src/metadata_table.h \
src/ml_model.h \
src/ml_optimiser.h \
src/multidim_array.h \
src/numerical_recipes.h \
src/parallel.h \
src/particle_polisher.h \
src/particle_sorter.h \
src/postprocessing.h \
src/preprocessing.h \
src/projector.h \
src/rwIMAGIC.h \
src/rwMRC.h \
src/rwSPIDER.h \
src/strings.h \
src/symmetries.h \
src/tabfuncs.h \
src/time.h \
src/transformations.h \
src/Healpix_2.15a/arr.h \
src/Healpix_2.15a/cxxutils.h \
src/Healpix_2.15a/datatypes.h \
src/Healpix_2.15a/geom_utils.h \
src/Healpix_2.15a/lsconstants.h \
src/Healpix_2.15a/message_error.h \
src/Healpix_2.15a/openmp_support.h \
src/Healpix_2.15a/pointing.h \
src/Healpix_2.15a/vec3.h \
src/Healpix_2.15a/healpix_base.h 

if HAVE_MPI
 librelion_@RELION_API_VERSION@_la_SOURCES += src/mpi.cpp \
 src/autopicker_mpi.cpp \
 src/ctffind_runner_mpi.cpp \
 src/particle_polisher_mpi.cpp \
 src/particle_sorter_mpi.cpp \
 src/preprocessing_mpi.cpp \
 src/ml_optimiser_mpi.cpp
 nobase_relion_include_HEADERS += src/mpi.h \
 src/autopicker_mpi.h \
 src/ctffind_runner_mpi.h \
 src/particle_polisher_mpi.h \
 src/particle_sorter_mpi.h \
 src/preprocessing_mpi.h \
 src/ml_optimiser_mpi.h
endif

if HAVE_FLTK
 librelion_@RELION_API_VERSION@_la_SOURCES += src/displayer.cpp \
 src/gui_entries.cpp \
 src/gui_jobwindow.cpp \
 src/gui_mainwindow.cpp \
 src/manualpicker.cpp 
 nobase_relion_include_HEADERS += src/displayer.h \
 src/gui_entries.h \
 src/gui_jobwindow.h \
 src/gui_mainwindow.h \
 src/manualpicker.h 
endif

## Instruct libtool to include ABI version information in the generated shared
## library file (.so).  The library ABI version is defined in configure.ac, so
## that all version information is kept in one place.
librelion_@RELION_API_VERSION@_la_LDFLAGS = -version-info $(RELION_SO_VERSION) $(FFT_LIBS)

## The generated configuration header is installed in its own subdirectory of
## $(libdir).  The reason for this is that the configuration information put
## into this header file describes the target platform the installed library
## has been built for.  Thus the file must not be installed into a location
## intended for architecture-independent files, as defined by the Filesystem
## Hierarchy Standard (FHS).
relion_libincludedir = $(libdir)/relion-$(RELION_API_VERSION)/include

## Install the generated pkg-config file (.pc) into the expected location for
## architecture-dependent package configuration information.  Occasionally,
## pkg-config files are also used for architecture-independent data packages,
## in which case the correct install location would be $(datadir)/pkgconfig.
pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = relion-$(RELION_API_VERSION).pc


###############################################################################
#
# Applications
#
###############################################################################

bin_PROGRAMS = project reconstruct refine preprocess postprocess autopick particle_polish particle_sort run_ctffind star_compare find_tiltpairs tiltpair_plot
if HAVE_MPI
 bin_PROGRAMS += refine_mpi preprocess_mpi autopick_mpi particle_polish_mpi particle_sort_mpi run_ctffind_mpi 
endif
if HAVE_FLTK
 bin_PROGRAMS += manualpick display maingui 	
endif

RelionLibs = librelion-$(RELION_API_VERSION).la 

refine_SOURCES = src/apps/refine.cpp
refine_LDADD = $(RelionLibs) 

refine_mpi_SOURCES = src/apps/refine_mpi.cpp
refine_mpi_LDADD = $(RelionLibs)

reconstruct_SOURCES = src/apps/reconstruct.cpp
reconstruct_LDADD = $(RelionLibs)

project_SOURCES = src/apps/project.cpp
project_LDADD = $(RelionLibs)

postprocess_SOURCES = src/apps/postprocess.cpp
postprocess_LDADD = $(RelionLibs)

preprocess_SOURCES = src/apps/preprocess.cpp
preprocess_LDADD = $(RelionLibs)

preprocess_mpi_SOURCES = src/apps/preprocess_mpi.cpp
preprocess_mpi_LDADD = $(RelionLibs)

autopick_SOURCES = src/apps/autopick.cpp
autopick_LDADD = $(RelionLibs)

autopick_mpi_SOURCES = src/apps/autopick_mpi.cpp
autopick_mpi_LDADD = $(RelionLibs)

manualpick_SOURCES = src/apps/manualpick.cpp
manualpick_LDADD = $(RelionLibs)

particle_polish_SOURCES = src/apps/particle_polish.cpp
particle_polish_LDADD = $(RelionLibs)

particle_polish_mpi_SOURCES = src/apps/particle_polish_mpi.cpp
particle_polish_mpi_LDADD = $(RelionLibs)

particle_sort_SOURCES = src/apps/particle_sort.cpp
particle_sort_LDADD = $(RelionLibs)

particle_sort_mpi_SOURCES = src/apps/particle_sort_mpi.cpp
particle_sort_mpi_LDADD = $(RelionLibs)

run_ctffind_SOURCES = src/apps/run_ctffind.cpp
run_ctffind_LDADD = $(RelionLibs)

run_ctffind_mpi_SOURCES = src/apps/run_ctffind_mpi.cpp
run_ctffind_mpi_LDADD = $(RelionLibs)

star_compare_SOURCES = src/apps/star_compare.cpp
star_compare_LDADD = $(RelionLibs)

display_SOURCES = src/apps/display.cpp
display_LDADD = $(RelionLibs)

maingui_SOURCES = src/apps/maingui.cpp
maingui_LDADD = $(RelionLibs)

bin_PROGRAMS += mask_create
mask_create_SOURCES = src/apps/mask_create.cpp
mask_create_LDADD = $(RelionLibs)

bin_PROGRAMS += stack_create
stack_create_SOURCES = src/apps/stack_create.cpp
stack_create_LDADD = $(RelionLibs)

bin_PROGRAMS += image_handler 
image_handler_SOURCES = src/apps/image_handler.cpp
image_handler_LDADD = $(RelionLibs)

bin_PROGRAMS += find_tiltpairs
find_tiltpairs_SOURCES = src/apps/find_tiltpairs.cpp
find_tiltpairs_LDADD = $(RelionLibs)

#bin_PROGRAMS += subtomo_proj3d_correct
#subtomo_proj3d_correct_SOURCES = src/apps/subtomo_proj3d_correct.cpp
#subtomo_proj3d_correct_LDADD = $(RelionLibs)

#bin_PROGRAMS += subtomo_angle_converter
#subtomo_angle_converter_SOURCES = src/apps/subtomo_angle_converter.cpp
#subtomo_angle_converter_LDADD = $(RelionLibs)

bin_PROGRAMS += tiltpair_plot
tiltpair_plot_SOURCES = src/apps/tiltpair_plot.cpp
tiltpair_plot_LDADD = $(RelionLibs)



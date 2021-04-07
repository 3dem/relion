---
name: Bug report
about: Report a problem
title: ''
labels: ''
assignees: ''

---

This is a template for reporting bugs. Please fill in as much information as you can.

**Describe your problem**

Please write a clear description of what the problem is.
Data processing questions should be posted to [the CCPEM mailing list](https://www.jiscmail.ac.uk/cgi-bin/webadmin?A0=CCPEM), not here.
**DO NOT** cross post a same question to multiple issues and/or many mailing lists (CCPEM, 3DEM, etc).

**Environment:**
 - OS: [e.g. Ubuntu 16.04 LTS]
 - MPI runtime: [e.g. OpenMPI 2.0.1]
 - RELION version [e.g. RELION-3.1-devel-commit-6ba935 (please see the title bar of the GUI)] 
 - Memory: [e.g. 128 GB]
 - GPU: [e.g. GTX 1080Ti]

**Dataset:**
 - Box size: [e.g. 256 px]
 - Pixel size: [e.g. 0.9 Å/px]
 - Number of particles: [e.g. 150,000]
 - Description: [e.g. A tetrameric protein of about 400 kDa in total]

**Job options:**
 - Type of job: [e.g. Refine3D]
 - Number of MPI processes: [e.g. 4]
 - Number of threads: [e.g. 6]
 - Full command (see `note.txt` in the job directory):
   ```
   `which relion_refine_mpi` --o Refine3D/job019/run --auto_refine --split_random_halves --i CtfRefine/job018/particles_ctf_refine.star --ref PostProcess/job001/postprocess.mrc --firstiter_cc --ini_high 12 --dont_combine_weights_via_disc --scratch_dir /ssd --pool 3 --pad 2  --ctf --ctf_corrected_ref --particle_diameter 142 --flatten_solvent --zero_mask --solvent_mask Result-by-Rado/run_class001_mask_th0.01_ns3_ngs7_box400.mrc --solvent_correct_fsc  --oversampling 1 --healpix_order 3 --auto_local_healpix_order 4 --offset_range 5 --offset_step 2 --sym O --low_resol_join_halves 40 --norm --scale  --j 8 --gpu "" --keep_scratch --pipeline_control Refine3D/job019/
   ```

**Error message:**

Please cite the *full* error message as the example below.

```
A line in the STAR file contains fewer columns than the number of labels. Expected = 3 Found = 2
Error in line: 0 0.0
in: /prog/relion-devel-lmb/src/metadata_table.cpp, line 966
=== Backtrace  ===
/prog/relion-devel-lmb/bin/relion_motion_refine_mpi(_ZN11RelionErrorC1ERKSsS1_l+0x41) [0x42e981]
/prog/relion-devel-lmb/bin/relion_motion_refine_mpi(_ZN13MetaDataTable12readStarLoopERSt14basic_ifstreamIcSt11char_traitsIcEEPSt6vectorI8EMDLabelSaIS6_EESsb+0xedd) [0x4361ad]
/prog/relion-devel-lmb/bin/relion_motion_refine_mpi(_ZN13MetaDataTable8readStarERSt14basic_ifstreamIcSt11char_traitsIcEERKSsPSt6vectorI8EMDLabelSaIS8_EESsb+0x580) [0x436f10]
/prog/relion-devel-lmb/bin/relion_motion_refine_mpi(_ZN10Micrograph4readE8FileNameb+0x5a3) [0x454bb3]
/prog/relion-devel-lmb/bin/relion_motion_refine_mpi(_ZN10MicrographC2E8FileNameS0_d+0x2e3) [0x4568b3]
/prog/relion-devel-lmb/bin/relion_motion_refine_mpi(_ZN17MicrographHandler14isMoviePresentERK13MetaDataTableb+0x180) [0x568280]
/prog/relion-devel-lmb/bin/relion_motion_refine_mpi(_ZN17MicrographHandler17cullMissingMoviesERKSt6vectorI13MetaDataTableSaIS1_EEi+0xe6) [0x568dc6]
/prog/relion-devel-lmb/bin/relion_motion_refine_mpi(_ZN13MotionRefiner4initEv+0x56f) [0x49e1ff]
/prog/relion-devel-lmb/bin/relion_motion_refine_mpi(main+0x31) [0x42a5e1]
/lib64/libc.so.6(__libc_start_main+0xf5) [0x2b7ac026e495]
/prog/relion-devel-lmb/bin/relion_motion_refine_mpi() [0x42b3cf]
==================
```

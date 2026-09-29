
# version 50001

data_scheme_general

_rlnSchemeName                       Schemes/amyproc/
_rlnSchemeCurrentNodeName            WAIT
 

# version 50001

data_scheme_floats

loop_ 
_rlnSchemeFloatVariableName #1 
_rlnSchemeFloatVariableValue #2 
_rlnSchemeFloatVariableResetValue #3 
 batchsize 100000.00000 100000.00000 
current_batch     1.000000     1.000000 
current_nr_mics     0.000000     0.000000 
maxtime_hr    48.000000    48.000000 
mybatchsize     0.000000     0.000000 
nr_batches     0.000000     0.000000 
       one     1.000000     1.000000 
prev_nr_mics     0.000000     0.000000 
  wait_sec   300.000000   300.000000 
 

# version 50001

data_scheme_bools

loop_ 
_rlnSchemeBooleanVariableName #1 
_rlnSchemeBooleanVariableValue #2 
_rlnSchemeBooleanVariableResetValue #3 
  has_mics            0            0 
has_mics_increased            0            0 
is_batch_full            1            1 
is_full_batch            0            0 
 

# version 50001

data_scheme_strings

loop_ 
_rlnSchemeStringVariableName #1 
_rlnSchemeStringVariableValue #2 
_rlnSchemeStringVariableResetValue #3 
batchnames Schemes/amyproc/split/particles_split*.star Schemes/amyproc/split/particles_split*.star 
input_mics Schemes/amyprep/select_mics/micrographs.star Schemes/amyprep/select_mics/micrographs.star 
micrographs micrographs micrographs 
   mybatch         ""         "" 
 mybatches         ""         "" 
 particles  particles  particles 
 

# version 50001

data_scheme_operators

loop_ 
_rlnSchemeOperatorName #1 
_rlnSchemeOperatorType #2 
_rlnSchemeOperatorOutput #3 
_rlnSchemeOperatorInput1 #4 
_rlnSchemeOperatorInput2 #5 
COUNT_mics float=count_images current_nr_mics input_mics micrographs 
COUNT_mybatchsize float=count_images mybatchsize    mybatch  particles 
EXIT_maxtime exit_maxtime  undefined maxtime_hr  undefined 
GLOB_batchnames string=glob  mybatches batchnames  undefined 
  HAS_mics bool=file_exists   has_mics input_mics  undefined 
HAS_mics_increased    bool=gt has_mics_increased current_nr_mics prev_nr_mics 
INCR_batch float=plus current_batch current_batch        one 
IS_batchfull    bool=eq is_full_batch mybatchsize  batchsize 
SET_mybatch string=nth_word    mybatch  mybatches current_batch 
SET_prev_nr_mics  float=set prev_nr_mics current_nr_mics  undefined 
      WAIT       wait  undefined   wait_sec  undefined 
 

# version 50001

data_scheme_jobs

loop_ 
_rlnSchemeJobNameOriginal #1 
_rlnSchemeJobName #2 
_rlnSchemeJobMode #3 
_rlnSchemeJobHasStarted #4 
autopick_trace autopick_trace   continue            0 
   class2d    class2d        new            0 
   extract    extract   continue            0 
     split      split   continue            0 
 

# version 50001

data_scheme_edges

loop_ 
_rlnSchemeEdgeInputNodeName #1 
_rlnSchemeEdgeOutputNodeName #2 
_rlnSchemeEdgeIsFork #3 
_rlnSchemeEdgeOutputNodeNameIfTrue #4 
_rlnSchemeEdgeBooleanVariable #5 
      WAIT EXIT_maxtime            0  undefined  undefined 
EXIT_maxtime   HAS_mics            0  undefined  undefined 
  HAS_mics       WAIT            1 COUNT_mics   has_mics 
COUNT_mics HAS_mics_increased            0  undefined  undefined 
HAS_mics_increased       WAIT            1 SET_prev_nr_mics has_mics_increased 
      WAIT EXIT_maxtime            0  undefined  undefined 
SET_prev_nr_mics autopick_trace            0  undefined  undefined 
autopick_trace    extract            0  undefined  undefined 
   extract      split            0  undefined  undefined 
     split GLOB_batchnames            0  undefined  undefined 
GLOB_batchnames SET_mybatch            0  undefined  undefined 
SET_mybatch COUNT_mybatchsize            0  undefined  undefined 
COUNT_mybatchsize IS_batchfull            0  undefined  undefined 
IS_batchfull       WAIT            1    class2d is_full_batch 
   class2d INCR_batch            0  undefined  undefined 
INCR_batch       WAIT            0  undefined  undefined 
 

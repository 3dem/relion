
# version 30001

data_schedule_general

_rlnScheduleName                       Schedules/preprocess/
_rlnScheduleCurrentNodeName            WAIT
 

# version 30001

data_schedule_floats

loop_ 
_rlnScheduleFloatVariableName #1 
_rlnScheduleFloatVariableValue #2 
_rlnScheduleFloatVariableResetValue #3 
do_at_most    50.000000    50.000000 
   maxtime    48.000000    48.000000 
  wait_sec   180.000000   180.000000 
 

# version 30001

data_schedule_bools

loop_ 
_rlnScheduleBooleanVariableName #1 
_rlnScheduleBooleanVariableValue #2 
_rlnScheduleBooleanVariableResetValue #3 
do_until_ctf            0            0 
has_topaz_model            0            0 
 

# version 30001

data_schedule_strings

loop_ 
_rlnScheduleStringVariableName #1 
_rlnScheduleStringVariableValue #2 
_rlnScheduleStringVariableResetValue #3 
topaz_model Schedules/refine/train_topaz/model_epoch10.sav Schedules/refine/train_topaz/model_epoch10.sav 
 

# version 30001

data_schedule_operators

loop_ 
_rlnScheduleOperatorName #1 
_rlnScheduleOperatorType #2 
_rlnScheduleOperatorOutput #3 
_rlnScheduleOperatorInput1 #4 
_rlnScheduleOperatorInput2 #5 
EXIT_maxtime exit_maxtime  undefined    maxtime  undefined 
HAS_topaz_model bool=file_exists has_topaz_model topaz_model  undefined 
      WAIT       wait  undefined   wait_sec  undefined 
 

# version 30001

data_schedule_jobs

loop_ 
_rlnScheduleJobNameOriginal #1 
_rlnScheduleJobName #2 
_rlnScheduleJobMode #3 
_rlnScheduleJobHasStarted #4 
   ctffind    ctffind   continue            0 
extract_logpick extract_logpick   continue            0 
extract_topazpick extract_topazpick   continue            0 
importmovies importmovies   continue            0 
 logpicker  logpicker   continue            0 
motioncorr motioncorr   continue            0 
split_logpick split_logpick   continue            0 
topazpicker topazpicker   continue            0 
 

# version 30001

data_schedule_edges

loop_ 
_rlnScheduleEdgeInputNodeName #1 
_rlnScheduleEdgeOutputNodeName #2 
_rlnScheduleEdgeIsFork #3 
_rlnScheduleEdgeOutputNodeNameIfTrue #4 
_rlnScheduleEdgeBooleanVariable #5 
      WAIT EXIT_maxtime            0  undefined  undefined 
EXIT_maxtime importmovies            0  undefined  undefined 
importmovies motioncorr            0  undefined  undefined 
motioncorr    ctffind            0  undefined  undefined 
   ctffind HAS_topaz_model            1       WAIT do_until_ctf 
HAS_topaz_model  logpicker            1 topazpicker has_topaz_model 
 logpicker extract_logpick            0  undefined  undefined 
extract_logpick split_logpick            0  undefined  undefined 
extract_logpick split_logpick            0  undefined  undefined 
split_logpick       WAIT            0  undefined  undefined 
topazpicker extract_topazpick            0  undefined  undefined 
extract_topazpick       WAIT            0  undefined  undefined 
 

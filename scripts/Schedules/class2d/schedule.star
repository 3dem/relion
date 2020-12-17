
# version 30001

data_schedule_general

_rlnScheduleName                       Schedules/class2d/
_rlnScheduleCurrentNodeName            WAIT
 

# version 30001

data_schedule_floats

loop_ 
_rlnScheduleFloatVariableName #1 
_rlnScheduleFloatVariableValue #2 
_rlnScheduleFloatVariableResetValue #3 
current_logbatch_size     0.000000     0.000000 
current_rest_size     0.000000     0.000000 
logbatch_size 10000.000000 10000.000000 
   maxtime    48.000000    48.000000 
prev_rest_size     0.000000     0.000000 
  wait_sec   180.000000   180.000000 
 

# version 30001

data_schedule_bools

loop_ 
_rlnScheduleBooleanVariableName #1 
_rlnScheduleBooleanVariableValue #2 
_rlnScheduleBooleanVariableResetValue #3 
has_larger_rest_size            0            0 
has_topaz_model            0            0 
logbatch_big_enough            0            0 
 

# version 30001

data_schedule_strings

loop_ 
_rlnScheduleStringVariableName #1 
_rlnScheduleStringVariableValue #2 
_rlnScheduleStringVariableResetValue #3 
  logbatch Schedules/preprocess/split_logpick/particles_split1.star Schedules/preprocess/split_logpick/particles_split1.star 
 particles  particles  particles 
rest_batch Schedules/preprocess/extract_topazpick/particles.star Schedules/preprocess/extract_topazpick/particles.star 
topaz_model Schedules/class2d/train_topaz/model_epoch10.sav Schedules/class2d/train_topaz/model_epoch10.sav 
 

# version 30001

data_schedule_operators

loop_ 
_rlnScheduleOperatorName #1 
_rlnScheduleOperatorType #2 
_rlnScheduleOperatorOutput #3 
_rlnScheduleOperatorInput1 #4 
_rlnScheduleOperatorInput2 #5 
CHECK_logbatch    bool=ge logbatch_big_enough current_logbatch_size logbatch_size 
COUNT_logbatch float=count_images current_logbatch_size   logbatch  particles 
COUNT_restbatch float=count_images current_rest_size rest_batch  particles 
EXIT_maxtime exit_maxtime  undefined    maxtime  undefined 
HAS_rest_increased    bool=gt has_larger_rest_size current_rest_size prev_rest_size 
HAS_topaz_model bool=file_exists has_topaz_model topaz_model  undefined 
SET_prev_rest_size  float=set prev_rest_size current_rest_size  undefined 
      WAIT       wait  undefined   wait_sec  undefined 
 

# version 30001

data_schedule_jobs

loop_ 
_rlnScheduleJobNameOriginal #1 
_rlnScheduleJobName #2 
_rlnScheduleJobMode #3 
_rlnScheduleJobHasStarted #4 
class2d_logbatch class2d_logbatch        new            0 
class2d_rest class2d_rest        new            0 
select_logbatch select_logbatch        new            0 
select_rest select_rest        new            0 
train_topaz train_topaz        new            0 
 

# version 30001

data_schedule_edges

loop_ 
_rlnScheduleEdgeInputNodeName #1 
_rlnScheduleEdgeOutputNodeName #2 
_rlnScheduleEdgeIsFork #3 
_rlnScheduleEdgeOutputNodeNameIfTrue #4 
_rlnScheduleEdgeBooleanVariable #5 
      WAIT EXIT_maxtime            0  undefined  undefined 
EXIT_maxtime HAS_topaz_model            0  undefined  undefined 
HAS_topaz_model COUNT_logbatch            1 COUNT_restbatch has_topaz_model 
COUNT_logbatch CHECK_logbatch            0  undefined  undefined 
CHECK_logbatch       WAIT            1 class2d_logbatch logbatch_big_enough 
class2d_logbatch select_logbatch            0  undefined  undefined 
select_logbatch train_topaz            0  undefined  undefined 
train_topaz       WAIT            0  undefined  undefined 
COUNT_restbatch HAS_rest_increased            0  undefined  undefined 
HAS_rest_increased       WAIT            1 class2d_rest has_larger_rest_size 
class2d_rest select_rest            0  undefined  undefined 
select_rest SET_prev_rest_size            0  undefined  undefined 
SET_prev_rest_size       WAIT            0  undefined  undefined 
 

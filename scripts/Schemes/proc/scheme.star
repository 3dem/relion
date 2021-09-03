# version 30001

data_scheme_general

_rlnSchemeName                       Schemes/proc/
_rlnSchemeCurrentNodeName            WAIT
 

# version 30001

data_scheme_floats

loop_ 
_rlnSchemeFloatVariableName #1 
_rlnSchemeFloatVariableValue #2 
_rlnSchemeFloatVariableResetValue #3 
current_nr_mics     0.000000     0.000000 
prev_nr_mics     0.000000     0.000000 
maxtime_hr    48.000000    48.000000 
wait_sec   180.000000   180.000000
current_nr_parts      0.000000     0.000000 
min_nr_parts_3d 5000.000   5000.000
 

# version 30001

data_scheme_bools

loop_ 
_rlnSchemeBooleanVariableName #1 
_rlnSchemeBooleanVariableValue #2 
_rlnSchemeBooleanVariableResetValue #3 
has_ctffind 0 0
do_3d 0 0
has_larger_nr_mics            0            0 
do_log            1    1
do_topaz          0   0
has_topaz_model            0            0 
has_iniref            0            0  
 

# version 30001

data_scheme_strings

loop_ 
_rlnSchemeStringVariableName #1 
_rlnSchemeStringVariableValue #2 
_rlnSchemeStringVariableResetValue #3 
ctffind_mics Schemes/prep/ctffind/micrographs_ctf.star Schemes/prep/ctffind/micrographs_ctf.star 
selected_mics Schemes/proc/select_mics/micrographs.star Schemes/proc/select_mics/micrographs.star
selected_parts Schemes/proc/select_parts/particles.star Schemes/proc/select_parts/particles.star
particles  particles  particles
micrographs  micrographs  micrographs 
topaz_model Schemes/proc/train_topaz/model_epoch10.sav Schemes/proc/train_topaz/model_epoch10.sav 
iniref None None
myref undefined undefined
inimodel_output Schemes/proc/inimodel3d/initial_model.mrc Schemes/proc/inimodel3d/initial_model.mrc

# version 30001

data_scheme_operators

loop_ 
_rlnSchemeOperatorName #1 
_rlnSchemeOperatorType #2 
_rlnSchemeOperatorOutput #3 
_rlnSchemeOperatorInput1 #4 
_rlnSchemeOperatorInput2 #5 
HAS_ctffind bool=file_exists has_ctffind ctffind_mics undefined
COUNT_mics float=count_images current_nr_mics   selected_mics  micrographs 
COUNT_parts float=count_images current_nr_parts   selected_parts particles
CHECK_do_3d bool=ge do_3d current_nr_parts min_nr_parts_3d
SET_prev_nr_mics  float=set prev_nr_mics current_nr_mics  undefined 
HAS_mics_increased    bool=gt has_larger_nr_mics current_nr_mics prev_nr_mics 
EXIT_maxtime exit_maxtime  undefined    maxtime_hr  undefined 
SET_do_topaz bool=not do_topaz do_log undefined
HAS_topaz_model bool=file_exists has_topaz_model topaz_model  undefined 
CHECK_iniref bool=file_exists has_iniref     iniref  undefined 
SET_myref_user string=set myref iniref undefined
SET_myref_inimodel string=set myref inimodel_output undefined
WAIT       wait  undefined   wait_sec  undefined 

# version 30001

data_scheme_jobs

loop_ 
_rlnSchemeJobNameOriginal #1 
_rlnSchemeJobName #2 
_rlnSchemeJobMode #3 
_rlnSchemeJobHasStarted #4 
select_mics select_mics continue    0
autopick autopick   continue            0 
extract extract   continue            0 
class2d class2d        new            0 
select_parts select_parts        new            0 
inimodel3d inimodel3d        new            0 
refine3d   refine3d        new            0 
 

# version 30001

data_scheme_edges

loop_ 
_rlnSchemeEdgeInputNodeName #1 
_rlnSchemeEdgeOutputNodeName #2 
_rlnSchemeEdgeIsFork #3 
_rlnSchemeEdgeOutputNodeNameIfTrue #4 
_rlnSchemeEdgeBooleanVariable #5 
WAIT EXIT_maxtime             0  undefined  undefined 
EXIT_maxtime HAS_ctffind      0  undefined  undefined 
HAS_ctffind WAIT             1 select_mics has_ctffind
select_mics COUNT_mics              0  undefined  undefined
COUNT_mics HAS_mics_increased              0  undefined  undefined 
HAS_mics_increased WAIT             1 SET_prev_nr_mics has_larger_nr_mics
SET_prev_nr_mics SET_do_topaz  0 undefined undefined
SET_do_topaz  autopick  0 undefined undefined
autopick extract           0  undefined  undefined 
extract class2d       0  undefined  undefined 
class2d select_parts            0  undefined  undefined 
select_parts COUNT_parts         0  undefined  undefined 
COUNT_parts CHECK_do_3d         0  undefined  undefined 
CHECK_do_3d   WAIT       1 CHECK_iniref  do_3d 
CHECK_iniref inimodel3d 1 SET_myref_user has_iniref
inimodel3d SET_myref_inimodel            0  undefined  undefined 
SET_myref_inimodel   refine3d            0  undefined  undefined 
SET_myref_user refine3d            0  undefined  undefined 
refine3d WAIT            0  undefined  undefined 

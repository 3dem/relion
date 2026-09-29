
# version 50001

data_scheme_general

_rlnSchemeName                       Schemes/amyprep/
_rlnSchemeCurrentNodeName            WAIT
 

# version 50001

data_scheme_floats

loop_ 
_rlnSchemeFloatVariableName #1 
_rlnSchemeFloatVariableValue #2 
_rlnSchemeFloatVariableResetValue #3 
current_nr_movies     0.000000     0.000000 
do_at_most 25000.000000 25000.000000 
maxtime_hr    48.000000    48.000000 
min_fom_skew     1.000000     1.000000 
prev_nr_movies     0.000000     0.000000 
  wait_sec   300.000000   300.000000 
 

# version 50001

data_scheme_bools

loop_ 
_rlnSchemeBooleanVariableName #1 
_rlnSchemeBooleanVariableValue #2 
_rlnSchemeBooleanVariableResetValue #3 
has_increased            0            0 
has_larger_nr_mics            1            1 
  has_mics            1            1 
has_movies            0            0 
 

# version 50001

data_scheme_strings

loop_ 
_rlnSchemeStringVariableName #1 
_rlnSchemeStringVariableValue #2 
_rlnSchemeStringVariableResetValue #3 
input_movies Schemes/amyprep/import/movies.star Schemes/amyprep/import/movies.star 
    movies     movies     movies 
 

# version 50001

data_scheme_operators

loop_ 
_rlnSchemeOperatorName #1 
_rlnSchemeOperatorType #2 
_rlnSchemeOperatorOutput #3 
_rlnSchemeOperatorInput1 #4 
_rlnSchemeOperatorInput2 #5 
COUNT_movies float=count_images current_nr_movies input_movies     movies 
EXIT_maxtime exit_maxtime  undefined maxtime_hr  undefined 
HAS_increased    bool=gt has_increased current_nr_movies prev_nr_movies 
SET_prev_nr_movies  float=set prev_nr_movies current_nr_movies  undefined 
      WAIT       wait  undefined   wait_sec  undefined 
 

# version 50001

data_scheme_jobs

loop_ 
_rlnSchemeJobNameOriginal #1 
_rlnSchemeJobName #2 
_rlnSchemeJobMode #3 
_rlnSchemeJobHasStarted #4 
autopick_fom autopick_fom   continue            0 
   ctffind    ctffind   continue            0 
    import     import   continue            0 
importmovies importmovies         ""            0 
motioncorr motioncorr   continue            0 
select_mics select_mics   continue            0 
 

# version 50001

data_scheme_edges

loop_ 
_rlnSchemeEdgeInputNodeName #1 
_rlnSchemeEdgeOutputNodeName #2 
_rlnSchemeEdgeIsFork #3 
_rlnSchemeEdgeOutputNodeNameIfTrue #4 
_rlnSchemeEdgeBooleanVariable #5 
      WAIT EXIT_maxtime            0  undefined  undefined 
EXIT_maxtime     import            0  undefined  undefined 
    import COUNT_movies            0  undefined  undefined 
COUNT_movies HAS_increased            0  undefined  undefined 
HAS_increased       WAIT            1 SET_prev_nr_movies has_larger_nr_mics 
SET_prev_nr_movies motioncorr            0  undefined  undefined 
motioncorr    ctffind            0  undefined  undefined 
   ctffind autopick_fom            0  undefined  undefined 
autopick_fom select_mics            0  undefined  undefined 
select_mics       WAIT            0  undefined  undefined 
 

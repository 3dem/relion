# The MIT License (MIT)
#
# Copyright (c) 2023 Dari Kimanius
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

# This script finds Conda environment and its Python executable
# and puts it into PYTHON_EXE_PATH.
# It requires that CONDA_ENV_NAME is set.
# It also find the Torch home and puts it into TORCH_HOME_PATH.


set(PYTHON_EXE_PATH "" CACHE PATH "Path to Python executable")

if(NOT PYTHON_EXE_PATH STREQUAL "" AND NOT EXISTS "${PYTHON_EXE_PATH}")
  message(WARNING "Provided path to the Python executable does not exist: ${PYTHON_EXE_PATH}")
  set(PYTHON_EXE_PATH "")
endif ()

if (PYTHON_EXE_PATH STREQUAL "")
  message(STATUS "Will try to find Conda environment (name: ${CONDA_ENV_NAME}) and its Python executable...")

  find_program(CONDA_EXECUTABLE conda)

  if(CONDA_EXECUTABLE)
    message(STATUS "Found Conda executable: ${CONDA_EXECUTABLE}")

    # Run the conda command and capture its output
    execute_process(COMMAND ${CONDA_EXECUTABLE} config --show envs_dirs
            OUTPUT_VARIABLE CONDA_ENVS_OUTPUT
            OUTPUT_STRIP_TRAILING_WHITESPACE
            RESULT_VARIABLE CONDA_RESULT)

    # Check if the conda command was successful
    if(CONDA_RESULT EQUAL 0)
      # Split the output into lines
      string(REPLACE "\n" ";" CONDA_ENVS_LIST ${CONDA_ENVS_OUTPUT})

      # Loop through the list of paths
      foreach(path ${CONDA_ENVS_LIST})
        string(STRIP "${path}" path)
        string(REPLACE " " ";" path "${path}")
        list(GET path -1 path)
        set(path "${path}/${CONDA_ENV_NAME}/bin/python")
        if(EXISTS "${path}")
          set(PYTHON_EXE_PATH "${path}")
          message(STATUS "Found Conda environment (name: ${CONDA_ENV_NAME}) and its Python executable.")
          break()  # Exit the loop once a valid path is found
        endif()
      endforeach()
    endif()
    if (PYTHON_EXE_PATH STREQUAL "")
      message(
              WARNING
              "Could NOT find path to Conda environment (name: ${CONDA_ENV_NAME}) and its Python executable.\n"
              "You can specify it directly with -DPYTHON_EXE_PATH=<path>\n"
      )
    endif()
  else (CONDA_EXECUTABLE)
    message(WARNING "Could NOT find Conda executable...")
  endif()
endif ()

set(TORCH_HOME_PATH "" CACHE PATH "Path to Torch home directory for storage of trained models.")

if(NOT TORCH_HOME_PATH STREQUAL "" AND NOT EXISTS "${TORCH_HOME_PATH}")
  message(WARNING "Provided path to the Torch home directory does not exist: ${TORCH_HOME_PATH}")
  set(TORCH_HOME_PATH "")
endif ()

message(STATUS "Using Python executable: ${PYTHON_EXE_PATH}")

if (NOT PYTHON_EXE_PATH STREQUAL "")
  message(STATUS "Will try to find Torch home directory...")

  if (TORCH_HOME_PATH STREQUAL "")
    execute_process(
            COMMAND ${PYTHON_EXE_PATH} -c "import torch; print(torch.hub._get_torch_home())"
            OUTPUT_VARIABLE TORCH_HOME_PATH_OUTPUT
            OUTPUT_STRIP_TRAILING_WHITESPACE
            RESULT_VARIABLE PYTHON_RESULT
    )
    if(EXISTS "${TORCH_HOME_PATH_OUTPUT}")
      set(TORCH_HOME_PATH ${TORCH_HOME_PATH_OUTPUT} CACHE PATH "Path to Torch home directory for storage of trained models." FORCE)
    else ()
      message(
              WARNING
              "Could NOT find Torch home directory for Conda environment (name: ${CONDA_ENV_NAME}).\n"
              "You can specify it directly with -DTORCH_HOME_PATH=<path>\n"
      )
    endif ()
  endif ()
endif ()

message(STATUS "Using Torch home: ${TORCH_HOME_PATH}")

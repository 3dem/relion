/***************************************************************************
 *
 * Author: "Sjors H.W. Scheres"
 * MRC Laboratory of Molecular Biology
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#ifndef PIPELINE_CONTROL_H_
#define PIPELINE_CONTROL_H_

#include <string>
#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include "src/macros.h"

// pipeliner
extern std::string pipeline_control_outputname;
#define RELION_JOB_EXIT_SUCCESS "RELION_JOB_EXIT_SUCCESS"
#define RELION_JOB_EXIT_FAILURE "RELION_JOB_EXIT_FAILURE"
#define RELION_JOB_EXIT_ABORTED "RELION_JOB_EXIT_ABORTED"
#define RELION_JOB_ABORT_NOW    "RELION_JOB_ABORT_NOW"

#define RELION_EXIT_SUCCESS pipeline_control_relion_exit(0)
#define RELION_EXIT_FAILURE pipeline_control_relion_exit(1)
#define RELION_EXIT_ABORTED pipeline_control_relion_exit(2)

int pipeline_control_relion_exit(int mode);

bool is_under_pipeline_control();

bool pipeline_control_check_abort_job();

void pipeline_control_delete_exit_files();

#endif /* PIPELINE_CONTROL_H_ */

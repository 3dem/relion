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
/***************************************************************************
 *
 * Authors:     Roberto Marabini (roberto@cnb.csic.es)
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.csic.es'
 ***************************************************************************/

#include "batchrunner_mpi.h"

//#define MAX_LINE 2048
//char szline[MAX_LINE+1];

void BatchRunnerMpi::read(int argc, char **argv)
{
	// Define a new MpiNode
	node = new MpiNode(argc, argv);

	if (node->isLeader())
		PRINT_VERSION_INFO();

	// First read in non-parallelisation-dependent variables
	BatchRunner::read(argc, argv);

	// Don't put any output to screen for mpi followers
	if (!node->isLeader())
		verb = 0;

	// Print out MPI info
	printMpiNodesMachineNames(*node);
}


void BatchRunnerMpi::run()
{
	MPI_Status status;
	if (node->rank == 0)
    {
    	int barstep, done=0;
    	if (verb > 0)
    	{
    		std::cout << " Executing " << number_of_commands << " commands ..." << std::endl;
    		init_progress_bar(number_of_commands);
    		barstep = XMIPP_MAX(1, number_of_commands / 60);
    	}

		FileName line;
        int number_of_node_waiting = 0; // max is nprocs -1
        for (int iline = 0; iline < mylines.size(); iline++)
        {

            //wait until a server is free
            int signal_back;
            MPI_Recv(&signal_back, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
            if (signal_back >= 0)
            {
            	// Write continuation file
            	done_lines[signal_back] = true;
            	writeContinuationFile();
            }

            number_of_node_waiting++;

            FileName line=mylines[iline];
            //strcpy(szline, line.c_str());

            std::string::size_type loc = line.find("MPI_BARRIER", 0);
            if (loc != std::string::npos)
            {
                while (number_of_node_waiting < (node->size - 1))
                {
                    MPI_Recv(0, 0, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
                    number_of_node_waiting++;
                }
                while (number_of_node_waiting > 0)
                {
                    MPI_Send(&iline, 1, MPI_INT, number_of_node_waiting, MPITAG_WAIT, MPI_COMM_WORLD);
                    number_of_node_waiting--;
                }
                continue;
            }

            //send work: which line do I need to process?
            MPI_Send(&iline, 1, MPI_INT, status.MPI_SOURCE, MPITAG_JOB_REQUEST, MPI_COMM_WORLD);
            number_of_node_waiting--;

            done++;
    		if (verb > 0 && done % barstep == 0) progress_bar(done);

        }

    	if (verb > 0) progress_bar(number_of_commands);

        for (size_t i = 1; i < node->size; i++)
            MPI_Send(0, 0, MPI_INT, i, MPITAG_STOP, MPI_COMM_WORLD);
        if (verb > 0) std::cout << " done! " << std::endl;
    }
    else
    {
        int signal_back = -1;
    	while (1)
        {
            //I am free
            MPI_Send(&signal_back, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);

            //get your next task
            MPI_Probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

            if (status.MPI_TAG == MPITAG_STOP) //I am free
            {
            	break;
            }
            else if (status.MPI_TAG == MPITAG_WAIT) //wait
            {
                int tmp;
                MPI_Recv(&tmp, 1, MPI_INT, 0, MPITAG_WAIT, MPI_COMM_WORLD, &status);
                continue;
            }
            else if (status.MPI_TAG == MPITAG_JOB_REQUEST) //work to do
            {
                MPI_Recv(&signal_back, 1, MPI_INT, 0, MPITAG_JOB_REQUEST, MPI_COMM_WORLD, &status);
                // do the job
                executeCommand(mylines[signal_back], node->rank);
            }
            else
            {
            	std::cerr << "WRONG TAG RECEIVED" << std::endl;
            }

        }
    }

}

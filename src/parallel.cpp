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
 * Authors:     J.M. De la Rosa Trevin (jmdelarosa@cnb.csic.es)
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
#include "src/parallel.h"

// ================= MUTEX ==========================
Mutex::Mutex()
{
    omp_init_lock(&mutex);
}

void Mutex::lock()
{
    omp_set_lock(&mutex);
}

Mutex::~Mutex()
{
    omp_destroy_lock(&mutex);
}

void Mutex::unlock()
{
    omp_unset_lock(&mutex);
}

// =================== TASK_DISTRIBUTOR ============================

ParallelTaskDistributor::ParallelTaskDistributor(size_t nTasks, size_t bSize)
{

	resize(nTasks, bSize);
}

void ParallelTaskDistributor::resize(size_t nTasks, size_t bSize)
{
	if (!(nTasks && bSize && bSize <= nTasks))
	{
		std::cerr << " nTasks= " << nTasks << " bSize= " << bSize << std::endl;
		REPORT_ERROR("nTasks and bSize should be > 0, also bSize <= nTasks");
	}

	numberOfTasks = nTasks;
	blockSize = bSize;
	assignedTasks = 0;

}

void ParallelTaskDistributor::reset()
{
    lock();
    assignedTasks = 0;
    unlock();
}

void ParallelTaskDistributor::setBlockSize(size_t bSize)
{
    lock();
    blockSize = bSize;
    unlock();
}

int ParallelTaskDistributor::getBlockSize() const
{
    return blockSize;
}

bool ParallelTaskDistributor::getTasks(size_t &first, size_t &last)
{
    lock();
    bool result = distribute(first, last);
    unlock();
    return result;
}

bool ParallelTaskDistributor::setAssignedTasks(size_t tasks)
{
    if (tasks < 0 || tasks >= numberOfTasks)
        return false;
    lock();
    assignedTasks = tasks;
    unlock();
    return true;
}

void ThreadTaskDistributor::lock()
{
    mutex.lock();
}

void ThreadTaskDistributor::unlock()
{
    mutex.unlock();
}

bool ThreadTaskDistributor::distribute(size_t &first, size_t &last)
{
    bool result = true;
    first = last = 0;
    if (assignedTasks >= numberOfTasks)
    {
        result = false;
    }
    else
    {
        first = assignedTasks;
        assignedTasks
        = (assignedTasks + blockSize < numberOfTasks) ? (assignedTasks
                + blockSize) : numberOfTasks;
        last = assignedTasks - 1;
    }
    return result;
}

/** Divides a number into most equally groups */
long int divide_equally(long int N, int size, int rank, long int &first, long int &last)
{
    long int jobs_per_worker = N / size;
    long int jobs_resting = N % size;
    if (rank < jobs_resting)
    {
        first = rank * (jobs_per_worker + 1);
        last = first + jobs_per_worker;
    }
    else
    {
        first = rank * jobs_per_worker + jobs_resting;
        last = first + jobs_per_worker - 1;
    }
    return last - first + 1;
}

/** In which group from divide_equally is myself? */
int divide_equally_which_group(long int N, int size, long int myself)
{
    long int first, last;
    for (int rank = 0; rank < size; rank++)
    {
        divide_equally(N, size, rank, first, last);
        if (myself >= first && myself <= last)
            return rank;
    }
    return -1;
}


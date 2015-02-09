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
    pthread_mutex_init(&mutex, NULL);
}

void Mutex::lock()
{
    pthread_mutex_lock(&mutex);
}

Mutex::~Mutex()
{
    pthread_mutex_destroy(&mutex);
}

void Mutex::unlock()
{
    pthread_mutex_unlock(&mutex);
}

// ================= BARRIER ==========================
Barrier::Barrier(int numberOfThreads)
{
    needed = numberOfThreads + 1;
    called = 0;
    pthread_mutex_init(&mutex, NULL);
    pthread_cond_init(&cond, NULL);
}

Barrier::~Barrier()
{
    pthread_mutex_destroy(&mutex);
    pthread_cond_destroy(&cond);
}

void Barrier::wait()
{
    pthread_mutex_lock(&mutex);
    ++called;
    if (called == needed)
    {
        called = 0;
        pthread_cond_broadcast(&cond);
    }
    else
    {
        pthread_cond_wait(&cond, &mutex);
    }
    pthread_mutex_unlock(&mutex);
}

// ================= THREAD MANAGER =======================

ThreadArgument::ThreadArgument()
{
    thread_id = -1;
    manager = NULL;
    data = NULL;
}

ThreadArgument::ThreadArgument(int id, ThreadManager * manager, void * data)
{
    this->thread_id = id;
    this->manager = manager;
    this->data = data;
}

void * _threadMain(void * data)
{
    ThreadArgument * thArg = (ThreadArgument*) data;
    ThreadManager * thMgr = thArg->manager;

    while (true)
    {
        //Wait for start working or leave
        thMgr->wait();
        //After awaked check what to do
        if (thMgr->workFunction != NULL)
        {
        	try
        	{
        		thMgr->workFunction(*thArg);
        		thMgr->wait(); //wait for finish together
        	}
        	catch (RelionError XE)
        	{
        		std::cerr << XE << std::endl
        		          << "In thread " << thArg->thread_id << std::endl;
        		pthread_exit(NULL);
        	}
        }
        else //exit thread
        {
            pthread_exit(NULL);
        }
    }
}

ThreadManager::ThreadManager(int numberOfThreads, void * workClass)
{
    threads = numberOfThreads;
    barrier = new Barrier(threads);
    workFunction = NULL;
    ids = new pthread_t[threads];
    arguments = new ThreadArgument[threads];
    started = false;
    this->workClass = workClass;

}

void ThreadManager::startThreads()
{
    //Create threads
    int result;

    for (int i = 0; i < threads; ++i)
    {
        arguments[i].thread_id = i;
        arguments[i].manager = this;
        arguments[i].data = NULL;
        arguments[i].workClass = workClass;

#ifdef DEBUG_THREADS
        std::cerr << " startThreads: i= " << i << std::endl;
#endif
        result = pthread_create(ids + i, NULL, _threadMain, (void*) (arguments + i));

        if (result != 0)
        {
            std::cerr << "ThreadManager: can't create threads." << std::endl;
            exit(1);
        }
    }
    started = true;
}

ThreadManager::~ThreadManager()
{
    //Destroy the threads
    workFunction = NULL;
    if (started)
    {
        wait();
        for (int i = 0; i < threads; ++i)
        	pthread_join(ids[i], NULL);
    }

    delete barrier;
    delete[] ids;
    delete[] arguments;
}

void ThreadManager::run(ThreadFunction function)
{
	runAsync(function);
    //Wait on barrier to wait for threads finish
    wait();
}

void ThreadManager::runAsync(ThreadFunction function)
{
    workFunction = function;
    if (!started)
        startThreads();
    //Wait on barrier to threads starts working
    wait();
}

void ThreadManager::wait()
{
    barrier->wait();
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


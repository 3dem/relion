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
#ifndef __PARALLEL_H
#define __PARALLEL_H

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "src/error.h"

// This code was copied from a developmental version of Xmipp-3.0
// which is developed at the Biocomputing Unit of the National Center for Biotechnology - CSIC
// in Madrid , Spain

/** Class wrapping around the OpenMPI mutex.
 * This class will provide a more object oriented implementation
 * of a mutex, to ensure the unique access to critical regions of
 * code and other syncronization problems.
 */
class Mutex
{
private:
    omp_lock_t mutex; //our OpenMPI mutex

public:
    /** Default constructor.
     * This constructor just initialize the OpenMPI mutex structure
     * with its defaults values, just like static initialization with PTHREAD_MUTEX_INITIALIZER
     */
    Mutex();

    /** Destructor. */
    virtual ~Mutex();

    /** Function to get the access to the mutex.
     * If some thread has the mutex and an other one
     * asks to lock it, it will be waiting until the first one
     * releases the mutex
     */
    virtual void lock();

    /** Function to release the mutex.
     * This allows access to the mutex to other
     * threads that are waiting for it.
     */
    virtual void unlock();
};//end of class Mutex

/** This class distributes dynamically N tasks between parallel workers.
 * @ingroup ParallelLibrary
 * This class is a generalization of a common task in a parallel
 * environment of dynamically distribute N tasks between workers(threads or mpi proccess).
 * Each worker will ask for a group of tasks, proccess it and ask for more tasks
 * until there is not more task to process.
 *
 * This class is abstract and only serves as base for
 * concrete implementations, which will provides the
 * specific lock mechanisms and the way of distribution.
 */
class ParallelTaskDistributor
{
//protected:
public:
    //How many tasks give in each request
    size_t blockSize;
    //The number of tasks that have been assigned
    size_t assignedTasks;

public:
    //The total number of tasks to be distributed
    size_t numberOfTasks;
    /** Constructor.
     * The number of jobs and block size should be provided.
     */
    ParallelTaskDistributor(size_t nTasks, size_t bSize);

    /* Resets the size of the task distributor
     */
    virtual void resize(size_t nTasks, size_t bSize);

    /** Destructor.
     */
    virtual ~ParallelTaskDistributor()
    {}
    ;

    /** Restart the number of assigned tasks and distribution again. (Resets assignedTasks = 0)
     * This method should only be called in the main thread
     * before start distributing the tasks between the workers
     * threads.
     */
    void reset();

    /** Set the number of tasks assigned in each request */
    void setBlockSize(size_t bSize);

    /** Return the number of tasks assigned in each request */
    int getBlockSize() const;

    /** Gets parallel tasks.
     *  @ingroup ParallelJobHandler
     *  This function will be called by workers for asking tasks
     *  until there are not more tasks to process.
     *  Example:
     *  @code
     *  //...
     *  ParallelTaskDistributor * td = new ThreadTaskDistributor(1000, 100);
     *  //...
     *  //function to perform some operation
     *  //to N images executed in parellel
     *  void processSeveralImages()
     *  {
     *      size_t firstImage, lastImage;
     *      while (td->getTasks(firstImage, lastImage))
     *          for (size_t image = firstImage; image <= lastImage; ++image)
     *          {
     *              //...
     *              processOneImage(image);
     *              //...
     *          }
     *  }
     *  @endcode
     */
    bool getTasks(size_t &first, size_t &last); // False = no more jobs, true = more jobs
    /* This function set the number of completed tasks.
     * Usually this not need to be called. Its more useful
     * for restarting work, when usually the leader detects
     * the number of tasks already done.
     */
    bool setAssignedTasks(size_t tasks);

protected:
    //Virtual functions that should be implemented in
    //subclasses, providing a mechanism of lock and
    //the specific way of distribute tasks.
    virtual void lock() = 0;
    virtual void unlock() = 0;
    virtual bool distribute(size_t &first, size_t &last) = 0;

};//class ParallelTaskDistributor

/** This class is a concrete implementation of ParallelTaskDistributor for POSIX threads.
 * It uses mutex as the locking mechanism
 * and distributes tasks from 0 to numberOfTasks.
 */
class ThreadTaskDistributor: public ParallelTaskDistributor
{
public:
    ThreadTaskDistributor(size_t nTasks, size_t bSize):ParallelTaskDistributor(nTasks, bSize) {}
    virtual ~ThreadTaskDistributor(){};

protected:
    Mutex mutex; ///< Mutex to syncronize access to critical region
    virtual void lock();
    virtual void unlock();
    virtual bool distribute(size_t &first, size_t &last);
};//end of class ThreadTaskDistributor

/// @name Miscellaneous functions
//@{
/** Divides a number into most equally groups
 *
 * For example you want to distribute N jobs between M workers
 * so each worker will have N/M jobs and some of them(N % M first)
 * will have N/M + 1 jobs
 * So for the worker 'rank' will be computed the first and last job to do
 * Return the number of jobs assigned, that could be N/M + 1 or N/M
 *
 */
long int divide_equally(long int N, int size, int rank, long int &first, long int &last);

/** In which group (of divide_equally) is myself situated?
 */
int divide_equally_which_group(long int N, int size, long int myself);

// Class which use local scope locks on higher-scope mutexes,
// resulting in zero risk of leaving locks on. Effectively,
// a mutex is treated like an container capable of holding a
// single Lock. We create a Lock and put it inside the mutex,
// and as soon as the Lock goes out out scope, the mutex is
// emptied (unlocked).
class Lock
{
public:
    explicit Lock(omp_lock_t *pm)
    : mutexPtr(pm)
    { omp_set_lock(mutexPtr); }

    ~Lock() { omp_unset_lock(mutexPtr); }

private:
    omp_lock_t *mutexPtr;
};

#endif

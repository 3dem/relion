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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * This complete copyright notice must be included in any revised version of the
 * source code. Additional authorship citations may be added, but existing
 * author citations must be preserved.
 ***************************************************************************/

#ifndef TIME_H_
#define TIME_H_

#include <sys/time.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <climits>
#include <vector>
#include <typeinfo>

// For timing functions
// Uncomment next line timing functions are giving problems in your system
//#define _NO_TIME
#ifndef _NO_TIME
#include <unistd.h>
#include <sys/times.h>
#ifdef _IRIX65
#include <sys/types.h>
#include <time.h>
#endif
#endif

/** @name Time managing
 *
 * These functions are used to make time measures of the algorithms. If you know
 * the total amount of work to do then some estimation can be done about how
 * much time is left before finishing. The time functions are very machine
 * dependent, we've tried to accomodate the compilation for several machines,
 * but if still programs do not work, you may configure Xmipp to avoid these
 * time measurements, functions are then substituted by null functions doing
 * nothing.
 *
 * @code
 * // Variable declaration
 * TimeStamp t0;
 *
 * // Beginning of the program
 * time_config();
 * ...
 *
 * annotate_time(&t0);
 * // Part to be measured
 * ...
 *
 * // End of part to be measured
 * print_elapsed_time(t0);
 * @endcode
 *
 * While for an estimation of time to go you can make it in two ways:  one
 * analytical, and another graphical.
 *
 * Analytical:
 *
 * @code
 * // Variable declaration
 * TimeStamp t0;
 * float to_go;
 *
 * // Beginning of the program
 * time_config();
 * ...
 *
 * annotate_time(&t0);
 * // Part to be measured
 * for (int i=0; i<60; i++)
 * {
 *	   ...
 *	   // Compute the time to go with the fraction of work already done
 *	   to_go = time_to_go(t0, (float) (i + 1) / 60);
 *	   std::cout << "I think you will be here " << to_go << "seconds more\n";
 * }
 * @endcode
 *
 * Graphical:
 * @code
 * // Beginning of the program
 * time_config();
 * ...
 *
 * // Init the progress bar with the total amount of work to do
 * // It is very important that there is no print out to stdout but
 * // the progress bar
 * init_progress_bar(60);
 *
 * // Part to be measured
 * for (int i=0; i<60; i++)
 * {
 *	   ...
 *	   progress_bar(i+1);
 * }
 *
 * // In this case the following call is useless since it has been
 * // already done in the loop, but there are cases where a final call
 * // with the total amount of work is not performed and although the
 * // whole task has been finished it seems that it hasn't as the
 * // progress bar hasn't been called with the final work but with
 * // a quantity a little smaller.
 * progress_bar(60);
 * @endcode
 *
 */
//@{
typedef struct tms TimeStamp; // Renaming of the time structure
/** Read the system clock frequency
 *
 * This operation is needed only once in a program always we want to have a time
 * measure, or an estimation of remaining time.
 *
 * @code
 * time_config();
 * @endcode
 *
 */
void time_config();
/** Annotate actual time
 *
 * This annotation is used later to compute the elapsed time.
 *
 * @code
 * TimeStamp t0;
 * annotate_time(&t0);
 * @endcode
 *
 */
void annotate_time(TimeStamp* time);
/** Acumulate time
 *
 * Initially dest_time should be set to orig time. Then you acumulate succesive
 * times calling this function (Destination time=destination_time + (now -
 * original time)) and finally the elapsed time is the dest time minus the first
 * one (the one which initiliazed the dest time.
 *
 */
void acum_time(TimeStamp* orig, TimeStamp* dest);
/** Compute elapsed time since a given annotation
 *
 * Given an annotation of time, this function computes the time elapsed since
 * then in seconds. The annotation is not modified. Usually the time is shown in
 * seconds, but you might specify to show it in clock ticks setting the variable
 * _IN_SECS to FALSE.
 *
 * @code
 * TimeStamp t0;
 * annotate_time(&t0);
 * ...;
 * float elapsed = elapsed_time(t0);
 *
 * TimeStamp t0;
 * annotate_time(&t0);
 * ...;
 * float elapsed = elapsed_time(t0, FALSE);
 * @endcode
 *
 */
float elapsed_time(TimeStamp& time, bool _IN_SECS = true);
/** Show on screen the elapsed time since a given annotation
 *
 * The format of the printing is "Elapsed time: User(13) System(1)" that means
 * that the user has used 13 seconds and the system 1, a total of 14 seconds
 * since the last annotation in this TimeStamp variable.
 *
 * @code
 * TimeStamp t0;
 * annotate_time(&t0);
 * ...;
 * print_elapsed_time(t0);
 * @endcode
 *
 * Usually the time is shown in seconds, but you might specify to show it in
 * clock ticks setting the variable _IN_SECS to FALSE.
 *
 */
void print_elapsed_time(TimeStamp& time, bool _IN_SECS = true);
/** Returns the estimated time left to finish
 *
 * To make this estimation the starting time must have been annotated before and
 * the fraction of the total amount of work must be estimated by the programmer.
 * See Time managing for an example.
 *
 */
float time_to_go(TimeStamp& time, float fraction_done);
/** Initialise the progress bar
 *
 * The progress bar is initialised to count for a total amount of work. For
 * instance, if we are to do something 60 times, the progress bar should be
 * initialised to that value. At the same time the bar is printed with the
 * initial guess of time left (ie, nothing "0000/????"). The number before the
 * slash is the elapsed time since initialisation of the progress bar, while the
 * second number is the estimation of total time that this task will take. See
 * Time managing for a more detailed example.
 *
 * @code
 * init_progress_bar(60);
 * @endcode
 */
void init_progress_bar(long total);
/** Update progress bar
 *
 * With this function you can change the already done amount of work, if
 * something is to be done 60 times and now we have already done 13 then we
 * could tell this to the progress bar with
 *
 * @code
 * progress_bar(13);
 * @endcode
 *
 * The information that this thing was to be done 60 times was given at the
 * initialisation of the progress bar. It is very important that during the use
 * of the progress bar, nobody prints anything to stdout as it is being used by
 * the progress bar. At the end you could make a call to progress_bar with the
 * total amount of work just to make sure that the printout is pretty enough.
 */
void progress_bar(long act_time);


/* Class to do some profiling
 *
 */
class Timer
{
public:
	///Start times for all individual timers
	std::vector<timeval> start_times;

	// General end time
	timeval end_time;

	// How many times has each tic/toc been called.
	std::vector<int> counts;

	// Total number of microseconds
	std::vector< long int> times;

	// Labels
	std::vector<std::string> tags;

	Timer()
	{
		clear();
	}

	~Timer()
	{
		clear();
	}

	void clear();

	void initZero();

	int setNew(const std::string tag);

	void tic(int timer);

	void toc(int timer);

	void printTimes(bool doClear);
};


#endif /* TIME_H_ */

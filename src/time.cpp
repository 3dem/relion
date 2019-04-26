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

#include "src/time.h"

/* Time managing ----------------------------------------------------------- */
// A global ................................................................
int XmippTICKS;

// Time configuration ......................................................
// The clock frequency for each machine must be known
void time_config()
{
	XmippTICKS = sysconf(_SC_CLK_TCK);
}

// Annotate actual time ....................................................
void annotate_time(TimeStamp *time)
{
	times(time);
}

// Acumulative time
void acum_time(TimeStamp *orig, TimeStamp *dest)
{
	TimeStamp now;
	times(&now);
	(*dest).tms_utime += (*dest).tms_utime + (now.tms_utime - (*orig).tms_utime);
	(*dest).tms_stime += (*dest).tms_stime + (now.tms_utime - (*orig).tms_utime);
}

// Show elapsed time since last annotation .................................
void print_elapsed_time(TimeStamp &time, bool _IN_SECS)
{
	TimeStamp now;
	times(&now);
	float userTime = now.tms_utime - time.tms_utime;
	float sysTime = now.tms_stime - time.tms_stime;
	if (_IN_SECS)
	{
		userTime /= XmippTICKS;
		sysTime /= XmippTICKS;
	}
	std::cout << "Elapsed time: User(" << userTime << ") System(" << sysTime << ")\n";
}

// Calculate elapsed time since last annotation .............................
float elapsed_time(TimeStamp &time, bool _IN_SECS)
{
	TimeStamp now;
	times(&now);
	float userTime = now.tms_utime - time.tms_utime;
	float sysTime = now.tms_stime - time.tms_stime;
	if (_IN_SECS)
	{
		userTime /= XmippTICKS;
		sysTime /= XmippTICKS;
	}
	return userTime + sysTime;
}

// Compute the predicted time left .........................................
float time_to_go(TimeStamp &time, float fraction_done)
{
	TimeStamp now;
	times(&now);
	float totalTime = (now.tms_utime - time.tms_utime +
	                   now.tms_stime - time.tms_stime) / XmippTICKS;
	return totalTime*(1 - fraction_done) / fraction_done;
}

// Show a message with the time it is produced .............................
void TimeMessage(const std::string & message)
{
	struct tm *T;
	time_t	   seconds;

	if (time(&seconds) < 0)
		seconds = 0;
	T = localtime(&seconds);

	printf("%2d:%2d:%2d (day=%2d) =>%s ", T->tm_hour,
	       T->tm_min, T->tm_sec, T->tm_mday, message.c_str());
}

// Init progress bar
void init_progress_bar(long total)
{
	progress_bar(-(total));
}

// Show a bar with the progress in time ....................................
// When the input is negative then we are setting the progress bar, this
// will be the total of elements to process. Afterwards the call to this
// routine must be in ascending order, ie, 0, 1, 2, ... No. elements
void progress_bar(long rlen)
{
	static time_t startt, prevt;
	time_t currt;
	static long totlen;
	long t1, t2;
	int min, i, hour;
	float h1, h2, m1, m2;

	if (rlen == 0)
		return;
	currt = time(NULL);

	if (rlen < 0)
	{
		totlen = -rlen;
		prevt = startt = currt;
		fprintf(stdout, "000/??? sec ");
		fprintf(stdout, "~~(,_,\">");
		for (i = 1; i < 10; i++)
			fprintf(stdout, "      ");
		fprintf(stdout, "    [oo]");
		fflush(stdout);
	}
	else if (totlen > 0)
	{
		t1 = currt - startt; // Elapsed time
		t2 = (long)(t1 * (float)totlen / rlen); // Total time

		hour = 0;
		min = 0;
		if (t2 > 60)
		{
			m1 = (float)t1 / 60.0;
			m2 = (float)t2 / 60.0;
			min = 1;
			if (m2 > 60)
			{
				h1 = (float)m1 / 60.0;
				h2 = (float)m2 / 60.0;
				hour = 1;
				min = 0;
			}
			else
				hour = 0;
		}
		else
			min = 0;

		if (hour)
			fprintf(stdout, "\r%3.2f/%3.2f %s ", h1, h2, "hrs");
		else if (min)
			fprintf(stdout, "\r%3.2f/%3.2f %s ", m1, m2, "min");
		else
			fprintf(stdout, "\r%4u/%4u %s ", (int)t1, (int)t2, "sec");

		i = (int)(60 * (1 - (float)(totlen - rlen) / totlen));
		while (i--)
			fprintf(stdout, ".");
		fprintf(stdout, "~~(,_,\">");
		if (rlen == totlen)
		{
			fprintf(stdout, "\n");
			totlen = 0;
		}
		fflush(stdout);
		prevt = currt;
	}
}



void Timer::clear()
{
	start_times.clear();
	counts.clear();
	times.clear();
	tags.clear();
}

void Timer::initZero()
{
	for (int i = 0; i < counts.size(); i++)
	{
		counts[i] = 0;
		times[i] = 0;
	}
}

int Timer::setNew(const std::string tag)
{
	//std::cerr << " tag = " << tag << std::endl;
	start_times.push_back(end_time);
	counts.push_back(0);
	times.push_back(0);
	tags.push_back(tag);
	return start_times.size() - 1;
}

void Timer::tic(int timer)
{
	gettimeofday(&(start_times[timer]), NULL);
	counts[timer]++;
}

void Timer::toc(int timer)
{
	gettimeofday(&end_time, NULL);
	times[timer] += (end_time.tv_sec - start_times[timer].tv_sec) * 1000000 +
				   (end_time.tv_usec - start_times[timer].tv_usec);
}

void Timer::printTimes(bool doClear)
{
	for (int i = 0; i < tags.size(); i++)
	{
		if (counts[i] > 0)
		{
			std::cout.width(35);
			std::cout << std::left << tags[i] << ": " << (times[i]/1000)/1000.0 << " sec (" << times[i] / counts[i] << " microsec/operation)"<<std::endl;
		}
	}
	if (doClear)
		Timer::clear();
	else
		Timer::initZero();
}

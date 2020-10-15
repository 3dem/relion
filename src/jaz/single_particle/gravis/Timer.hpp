#ifndef __LIBGRAVIS_TIMER_HPP__
#define __LIBGRAVIS_TIMER_HPP__
#ifndef WIN32
#include <sys/time.h>
#else
#include <time.h>
#endif
#include <iostream>
#include <iomanip>

namespace gravis
{
  class Timer
  {
      friend std::ostream& operator<<(std::ostream& os, Timer& t);

    private:
      clock_t start_clock;
#ifndef WIN32
      timeval start_time;
#endif

    public:
      Timer() : start_clock(), start_time()
      {
        restart();
      }
      // Copy and assignment are fine

      double wall_time() const
      {
#ifndef WIN32
        timeval current;
        gettimeofday(&current, 0);
        return (current.tv_sec - start_time.tv_sec) + (current.tv_usec - start_time.tv_usec)*1e-6;
#else
        return -1;
#endif
      }

      inline double cpu_time() const
      {
        return ticks_to_seconds(ticks());
      }

      static inline double ticks_to_seconds(const clock_t& ticks)
      {
        return double(ticks) / double(CLOCKS_PER_SEC);
      }

      inline clock_t ticks() const
      {
        return (clock() - start_clock);
      }

      inline void restart()
      {
#ifndef WIN32
        gettimeofday(&start_time, 0);
#endif
        start_clock = clock();
      }

      ~Timer() { };

  };
}


//===========================================================================
// Allow timers to be printed to ostreams using the syntax 'os << t'
// for an ostream 'os' and a timer 't'.  For example, "cout << t" will
// print out the total amount of time 't' has been "running".

inline std::ostream& operator<<(std::ostream& os, gravis::Timer& t)
{
  double wall_time = double(t.wall_time());
  double cpu_time = t.cpu_time();

  double min_time=wall_time;
  if (cpu_time<wall_time) min_time=cpu_time;

  if (min_time<1)
  {
    os <<
       "[" << std::setw(3) << std::setprecision(0) << std::setiosflags(std::ios::fixed) << wall_time*1000.0 <<
       "/" << std::setw(4) << std::setprecision(0) << std::setiosflags(std::ios::fixed) << cpu_time*1000.0 <<
       "ms]" << std::setprecision(5);
  }
  else if (min_time<10)
  {
    os <<
       "[" << std::setw(4) << std::setprecision(2) << std::setiosflags(std::ios::fixed) << wall_time <<
       "/" << std::setw(4) << std::setprecision(2) << std::setiosflags(std::ios::fixed) << cpu_time <<
       "s]" << std::setprecision(5);
  }
  else if (min_time < 100)
  {
    os <<
       "[" << std::setw(4) << std::setprecision(1) << std::setiosflags(std::ios::fixed) << wall_time <<
       "/" << std::setw(4) << std::setprecision(1) << std::setiosflags(std::ios::fixed) << cpu_time <<
       "s]" << std::setprecision(5);
  }
  else if (min_time < 60*100)
  {
    os <<
       "[" << std::setw(4) << std::setprecision(1) << std::setiosflags(std::ios::fixed) << wall_time/60 <<
       "/" << std::setw(4) << std::setprecision(1) << std::setiosflags(std::ios::fixed) << cpu_time/60 <<
       "m]" << std::setprecision(5);
  }
  else
  {
    os <<
       "[" << std::setw(4) << std::setprecision(1) << std::setiosflags(std::ios::fixed) << wall_time/60/60 <<
       "/" << std::setw(4) << std::setprecision(1) << std::setiosflags(std::ios::fixed) << cpu_time/60/60 <<
       "h]" << std::setprecision(5);
  };

  return os;
}


//===========================================================================

#endif

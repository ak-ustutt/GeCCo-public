/*---------------------------------------------------------------------------*
 *                       timer routines                                      *
 *  atim_cs_ returns cpu and system time and is reasonably fast              *
 *  atim_csw_ return in addition the real time (wall-clock time)             *
 *  the call to gettimeofday is expensive, however, so do not call           *
 *  this routine to often                                                    *
 *---------------------------------------------------------------------------*/

#include <unistd.h>
#include <sys/times.h>
#include <sys/time.h>
#include <timing.h>

/*---------------------------------------------------------------------------*/

void atim_cs_(double * cpu, double * syst)
{
  struct  tms  buf;

  times(&buf);

  *cpu  = buf.tms_utime / clocks_per_second;
  *syst = buf.tms_stime / clocks_per_second;

}

void atim_csw_(double * cpu, double * syst, double * realt)
{
  struct  tms  buf;
  struct  timeval  timv;

  times(&buf);

  gettimeofday(&timv,NULL);

  *realt  = timv.tv_usec ;
  *realt /= 1000000.0;
  *realt += (timv.tv_sec-time_of_day_initial);

  *cpu  = buf.tms_utime / clocks_per_second;
  *syst = buf.tms_stime / clocks_per_second;

}

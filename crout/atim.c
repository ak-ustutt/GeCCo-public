/*---------------------------------------------------------------------------*
 *                  timer routine                                            *
 *---------------------------------------------------------------------------*/

#include <unistd.h>
#include <sys/times.h>
#include <sys/time.h>
#include <timing.h>

/*---------------------------------------------------------------------------*/

void    atim_(double * cpu, double * syst, double * realt)
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


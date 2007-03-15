/*-----------------------------------------------*
 *    set global parameters for timer            *
 *-----------------------------------------------*/
#include<timing.h>
#include<unistd.h>
#include<sys/time.h>

void init_time_(void)
{
  struct timeval tmv;
 
  /* set clocks per second (for cpu and system time) */
  clocks_per_second = sysconf(_SC_CLK_TCK);
  /* get inital offset for wall time measuring */
  gettimeofday(&tmv,NULL);
  time_of_day_initial = tmv.tv_sec; 
} 



#include <stdio.h>

int flush_(iunit)
   int *iunit;
{
   if (*iunit == 6) fflush(stdout);
}

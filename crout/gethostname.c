#include <unistd.h>
void gethostname_(hostname,lhost)
char *hostname;
int lhost;
{
   gethostname(hostname,lhost);
}

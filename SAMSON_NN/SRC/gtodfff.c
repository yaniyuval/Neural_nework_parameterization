/*
 Interface to getttimeofday for Fortran
 taken from http://www.ncsa.uiuc.edu/UserInfo/Resources/Hardware/IA32LinuxCluster/Doc/timing.html
 Re-adapted from timing_fgettod at http://www.hpc.uio.no/hpc/SoftwareDev/opttune/gettimeofday.f
 Adapted __cplusplus ifdefs from Marat Khairoutdinov's systemf function */

#include <stddef.h> 
#include <sys/time.h>

#ifdef __cplusplus 
extern "C" {
#endif

void gtodfff(int times[2])
{
    struct timeval tv;
    int rtn;
    rtn=gettimeofday(&tv, NULL);
    times[0] = tv.tv_sec;
    times[1] = tv.tv_usec;
}

void gtodfff_(int times[2])
{
    struct timeval tv;
    int rtn;
    rtn=gettimeofday(&tv, NULL);
    times[0] = tv.tv_sec;
    times[1] = tv.tv_usec;
}

#ifdef __cplusplus
}
#endif

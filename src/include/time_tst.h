#ifndef TIME_TST_H_
#define TIME_TST_H_

#include <sys/time.h>
#include <ctime>

uint64_t GetTimeMs64()
{
  /* Linux */
  struct timeval tv;

  gettimeofday(&tv, NULL);

  uint64_t ret = tv.tv_usec;
  /* Convert from micro seconds (10^-6) to milliseconds (10^-3) */
  ret /= 1000;

  /* Adds the seconds (10^0) after converting them to milliseconds (10^-3) */
  ret += (tv.tv_sec * 1000);
  return ret;
}

#endif // TIME_TST_H_
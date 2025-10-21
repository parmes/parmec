/*
   The MIT License (MIT)

   Copyright (c) 2015 Tomasz Koziara

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in all
   copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
 */

#include <sys/time.h>

#ifndef __timer__
#define __timer__

struct timing
{
  struct timeval time;
  double sec, total;
};

static inline void timerstart (struct timing *t)
{
  gettimeofday (&t->time, NULL);
}

static inline double timerend (struct timing *t)
{
  struct timeval newtime;
  gettimeofday (&newtime, NULL);
  t->sec = (((double)newtime.tv_sec - (double)t->time.tv_sec) * 1000000. +
      ((double)newtime.tv_usec - (double)t->time.tv_usec)) / 1000000.;
  t->total += t->sec;
  return t->sec;
}

#endif

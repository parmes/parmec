/*
The MIT License (MIT)

Copyright (c) 2017 Tomasz Koziara

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

#include <stdio.h>
#include "macros.h"

#ifndef __tms__
#define __tms__

typedef struct time_series TMS;

struct time_series
{
  REAL value; /* constant value => used if size == 0 */
  REAL (*points) [2]; /* vector of (time, value) pairs */
  int marker; /* index of the last read interval */
  int size; /* total number of pairs */
};

/* create a copy */
TMS* TMS_Copy (TMS *ts);

/* create time series */
TMS* TMS_Create (int size, REAL *times, REAL *values);

/* create time series from a text file */
TMS* TMS_File (char *path);

/* wrapper for a constant value */
TMS* TMS_Constant (REAL value);

/* create from another series through integration */
TMS* TMS_Integral (TMS *ts);

/* as above, but calculate derivative */
TMS* TMS_Derivative (TMS *ts);

/* obtain value at a specific time */
REAL TMS_Value (TMS *ts, REAL time);

/* free time series memory */
void TMS_Destroy (TMS *ts);

#endif

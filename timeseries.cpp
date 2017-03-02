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

#include <stdlib.h>
#include <string.h>
#include "timeseries.h"
#include "mem.h"

#define CHUNK 512

static int findmarker (REAL (*begin)[2], REAL (*end)[2], REAL time)
{
  REAL (*low)[2] = begin, (*high)[2] = end, (*mid)[2];

  while (low <= high)
  {
    mid = low + (high-low) / 2;

    if (time >= (*mid)[0] &&
        time < (*(mid+1))[0])
    {
      return (mid - begin);
    }
    else if (time < (*mid)[0])
    {
      high = mid - 1;
    }
    else
    {
      low = mid + 1;
    }
  }

  return 0;
}

static REAL linterp (REAL (*point) [2], REAL time)
{
  REAL dt = time - point[0][0];
  return point[0][1] + (point[1][1]-point[0][1]) * (dt / (point[1][0] - point[0][0]));
}

TMS* TMS_Copy (TMS *ts)
{
  TMS *out;

  ERRMEM (out = (TMS*)malloc (sizeof (TMS)));
  out->value = ts->value;
  ERRMEM (out->points = (REAL(*)[2])malloc (sizeof (REAL [2]) * ts->size));
  out->marker = ts->marker;
  out->size = ts->size;
  memcpy (out->points, ts->points, sizeof (REAL [2]) * ts->size);

  return out;
}

TMS* TMS_Create (int size, REAL *times, REAL *values)
{
  TMS *ts;
  int i;

  ASSERT (size > 1, "Time series size must be greater than one");
  ERRMEM (ts = (TMS*)malloc (sizeof (TMS)));
  ERRMEM (ts->points = (REAL(*)[2])malloc (sizeof (REAL [2]) * size));
  ts->marker = 0;
  ts->size = size;

  for (i = 0; i < size; i ++)
  {
    ts->points [i][0] = times [i];
    ts->points [i][1] = values [i];

    if (i)
    {
      ASSERT (times [i] > times [i-1], "Time must increase with index");
    }
  }

  return ts;
}

TMS* TMS_File (char *path)
{
  char buf [4096];
  FILE *fp;
  TMS *ts;
  int np;

  if (!(fp = fopen (path, "r"))) return NULL;
  ERRMEM (ts = (TMS*)malloc (sizeof (TMS)));
  ERRMEM (ts->points = (REAL(*)[2])MEM_CALLOC (sizeof (REAL [2]) * CHUNK));

  np = CHUNK;
  ts->size = 0;
  while (fgets(buf, sizeof(buf), fp) != NULL) /* read line */
  {
    if (buf[0] == '#') continue; /* skip comments */

    if (sscanf (buf, "%lf%lf", &ts->points [ts->size][0], &ts->points [ts->size][1]) == EOF) break;

    if (ts->size)
    {
      ASSERT (ts->points [ts->size][0] > ts->points [ts->size-1][0], "Time must increase with index");
    }

    if (++ ts->size >= np)
    {
      np += CHUNK;
      ERRMEM (ts->points = (REAL(*)[2])realloc (ts->points, sizeof (REAL [2]) * np));
    }
  }

  if (ts->size == 0)
  {
    free (ts->points);
    free (ts);
    ASSERT (0, "Empty time series file");
  }
  else
  {
    ERRMEM (ts->points = (REAL(*)[2])realloc (ts->points, sizeof (REAL [2]) * ts->size)); /* trim */
  }

  ts->marker = 0;

  return ts;
}

TMS* TMS_Constant (REAL value)
{
  TMS *ts;

  ERRMEM (ts = (TMS*)malloc (sizeof (TMS)));
  ts->value = value;
  ts->size = 0; /* indicates constant value */

  return ts;
}

TMS* TMS_Integral (TMS *ts)
{
  REAL (*pin) [2], (*pout) [2];
  TMS *out;
  int n;

  ASSERT (ts->size > 0, "Cannot integrate constant time series");
  ERRMEM (out = (TMS*)malloc (sizeof (TMS)));
  out->size = ts->size;

  if (out->size == 0)
  {
    free (out);
    return TMS_Constant ((ts->points[1][0] - ts->points[0][0]) * 0.5 *  (ts->points[0][1] + ts->points[1][1]));
  }

  ERRMEM (out->points = (REAL(*)[2])MEM_CALLOC (sizeof (REAL [2]) * out->size));

  out->marker = 0;
  pin = ts->points;
  pout = out->points;

  pout [0][0] = pin[0][0];
  pout [0][1] = pin[0][1];

  for (n = 1; n < out->size; n ++)
  {
    pout [n][0] = pin [n][0];
    pout [n][1] = pout [n-1][1] + (pin[n][0] - pin[n-1][0]) * 0.5 * (pin[n-1][1] + pin[n][1]);
  }

  return out;
}

TMS* TMS_Derivative (TMS *ts)
{
  REAL (*pin) [2], (*pout) [2];
  TMS *out;
  int n;

  if (ts->size == 0) return TMS_Constant (0.0);

  if (ts->size == 1)
  {
    return TMS_Constant ((ts->points[1][1] - ts->points[0][1])/(ts->points[1][0] - ts->points[0][0]));
  }

  ERRMEM (out = (TMS*)malloc (sizeof (TMS)));
  out->size = 2*(ts->size-1); /* constant value per interval with linear trend: two ends of the interval */

  ERRMEM (out->points = (REAL(*)[2])MEM_CALLOC (sizeof (REAL [2]) * out->size));

  out->marker = 0;
  pin = ts->points;
  pout = out->points;

  for (n = 0; n < ts->size-1; n ++)
  {
    REAL dt = pin[n+1][0] - pin[n][0], eps = 1E-10 * dt;

    pout [2*n][0] = pin[n][0];
    pout [2*n+1][0] = pin[n+1][0] - eps;
    pout [2*n][1] =
    pout [2*n+1][1] = (pin[n+1][1] - pin[n][1]) / dt;
  }

  return out;
}

REAL TMS_Value (TMS *ts, REAL time)
{
  REAL lo, hi;

  if (ts->size == 0) return ts->value;

  if (time < ts->points[0][0]) return ts->points[0][1];
  else if (time > ts->points[ts->size-1][0]) return ts->points[ts->size-1][1];

  lo = ts->points[ts->marker > 0 ? ts->marker - 1 : ts->marker][0];
  hi = ts->points[ts->marker < ts->size - 1 ? ts->marker + 1 : ts->marker][0];

  if (time < lo || time > hi)
  {
    ts->marker = findmarker (ts->points, ts->points + ts->size - 1, time);
  }
  else if (time >= lo && ts->marker &&
	   time < ts->points[ts->marker][0]) ts->marker --;
  else if (time >= ts->points[ts->marker+1][0] &&
	   time < hi && ts->marker < ts->size - 1) ts->marker ++;

  return linterp (&ts->points[ts->marker], time);
}

void TMS_Destroy (TMS *ts)
{
  if (ts->size == 0) free (ts);
  else
  {
    free (ts->points);
    free (ts);
  }
}

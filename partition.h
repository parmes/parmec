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

#ifndef __partition__
#define __partition__


#define CUTOFF 64 /* TODO: tune radix tree cutoff size */

#define LSIZE 96 /* TODO: tune leaf buffer size */

/* leaf ellipsoids */
struct leaf_data
{
  uniform int size;
  uniform int color[LSIZE];
  uniform int part[LSIZE];
  uniform int ell[LSIZE];
  uniform REAL center[3][LSIZE];
  uniform REAL radii[3][LSIZE];
  uniform REAL orient[9][LSIZE];
};

/* partitioning tree */
struct partitioning
{
  uniform REAL coord;
  uniform int dimension;
  uniform int left;
  uniform int right;

  uniform leaf_data * uniform data;
};

#endif

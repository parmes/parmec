/*
The MIT License (MIT)

Copyright (c) 2016 EDF Energy

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

#ifndef __spring__
#define __spring__

#define SPRBUF 8 /* TODO: tune spring buffer size */

/* master springs; global array is used */
struct master_spring
{
  uniform int sprid[SPRBUF]; /* spring number returned to user */
  uniform int master[SPRBUF]; /* first particle index */
  uniform int slave[SPRBUF]; /* second particle index */
  uniform REAL mpnt[6][SPRBUF]; /* spatial and referential master points */
  uniform REAL spnt[6][SPRBUF]; /* spatial and referential slave points */
  uniform REAL * uniform spring[2][2]; /* spring force lookup tables: [0][][] loading, [1][][] unloading */
  uniform int spridx[2][SPRBUF+1]; /* spring force lookup start index: [0][] loading, [1][] unloading */
  uniform int8 sprtype[SPRBUF]; /* spring force type: SPRFRC_ in constants.h */
  uniform int spring_size[2]; /* size of spring force lookup tables: [0] loading, [1] unloading */
  uniform REAL * uniform dashpot[2]; /* dashpot force lookup tables */
  uniform int dashidx[SPRBUF+1]; /* dashpot force lookup start index */
  uniform int dashpot_size; /* size of dashpot force lookup tables */
  uniform REAL dir[3][SPRBUF]; /* spring direction */
  uniform int8 dirup[SPRBUF]; /* spring direction update flag: SPRDIR_ in constants.h */
  uniform REAL stroke0[SPRBUF]; /* initial spring stroke */
  uniform REAL stroke[SPRBUF]; /* current stroke */
  uniform REAL orient[3][SPRBUF]; /* current spring orientation */
  uniform REAL frc[2][SPRBUF]; /* total and spring force magnitude */
  uniform int size; /* SPRBUF used size */

  uniform master_spring * uniform next; /* local list */
  uniform int lock; /* list update lock */
};

/* slave springs; they are created by
 * symmetrically coppying master springs */
struct slave_spring
{
  uniform int master[SPRBUF];
  uniform REAL point[3][SPRBUF];
  uniform REAL force[3][SPRBUF];
  uniform int size;

  uniform slave_spring * uniform next; /* local list */
  uniform int lock; /* list update lock */
};

#endif

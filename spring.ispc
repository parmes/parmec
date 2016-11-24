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

/* Contributors: Tomasz Koziara */

#include "macros.h"
#include "spring.h"

/* allocate global array of master springs */
export uniform master_spring * uniform master_spring_alloc (uniform master_spring * uniform old, uniform int nold, uniform int size)
{
  uniform master_spring * uniform spr = uniform new uniform master_spring [size];

  if (nold)
  {
    memcpy (spr, old, nold * sizeof (uniform master_spring));

    delete old;
  }

  for (uniform int i = nold; i < size; i ++)
  {
    spr[i].spring[0][0] = NULL;
    spr[i].spring[0][1] = NULL;
    spr[i].spring[1][0] = NULL;
    spr[i].spring[1][1] = NULL;
    spr[i].spridx[0][0] = 0;
    spr[i].spridx[1][0] = 0;
    spr[i].spring_size[0] = 0;
    spr[i].spring_size[1] = 0;
    spr[i].dashpot[0] = NULL;
    spr[i].dashpot[1] = NULL;
    spr[i].dashidx[0] = 0;
    spr[i].dashpot_size = 0;

    spr[i].size = 0;
    spr[i].next = NULL;
    spr[i].lock = -1;
  }
 
  return spr;
}

/* allocate new local master spring that can be written to */
export uniform master_spring * uniform newspring (uniform master_spring * uniform master, uniform int taskindex, uniform int * uniform k)
{
  if (taskindex >= 0)
  {
    while (master->lock != taskindex)
      atomic_compare_exchange_global (&master->lock, -1, taskindex); /* lock access */
  }

  uniform master_spring * uniform spr = master;

  while (spr->size == SPRBUF && spr->next != NULL) spr = spr->next; /* find available item or rewind to end */

  if (spr->size < SPRBUF)
  {
    *k = spr->size ++;
  }
  else
  {
    uniform master_spring * uniform ptr = uniform new uniform master_spring;

    ptr->spring[0][0] = NULL;
    ptr->spring[0][1] = NULL;
    ptr->spring[1][0] = NULL;
    ptr->spring[1][1] = NULL;
    ptr->spridx[0][0] = 0;
    ptr->spridx[1][0] = 0;
    ptr->spring_size[0] = 0;
    ptr->spring_size[1] = 0;
    ptr->dashpot[0] = NULL;
    ptr->dashpot[1] = NULL;
    ptr->dashidx[0] = 0;
    ptr->dashpot_size = 0;

    ptr->size = 0;
    ptr->next = NULL;
    spr->next = ptr; /* append new item at the end */
    spr = ptr; /* return new item */
    *k = 0;
  }

  if (taskindex >= 0) master->lock = -1; /* unlock */

  return spr;
}

/* allocate local data of master springs */
export void master_spring_alloc_data (uniform master_spring * uniform ptr, uniform int spring_inc[2], uniform int dashpot_inc)
{
  if (ptr->spring_size[0] < ptr->spridx[0][ptr->size-1] + spring_inc[0])
  {
    uniform int n = 2*(ptr->spridx[0][ptr->size-1] + spring_inc[0]);

    uniform REAL * uniform p0 = uniform new uniform REAL[n];
    uniform REAL * uniform p1 = uniform new uniform REAL[n];

    memcpy (p0, ptr->spring[0][0], ptr->spridx[0][ptr->size-1] * sizeof (uniform REAL));
    memcpy (p1, ptr->spring[0][1], ptr->spridx[0][ptr->size-1] * sizeof (uniform REAL));

    delete ptr->spring[0][0];
    delete ptr->spring[0][1];

    ptr->spring[0][0] = p0;
    ptr->spring[0][1] = p1;

    ptr->spring_size[0] = n;
  }

  if (ptr->spring_size[1] < ptr->spridx[1][ptr->size-1] + spring_inc[1])
  {
    uniform int n = 2*(ptr->spridx[1][ptr->size-1] + spring_inc[1]);

    uniform REAL * uniform p0 = uniform new uniform REAL[n];
    uniform REAL * uniform p1 = uniform new uniform REAL[n];

    memcpy (p0, ptr->spring[1][0], ptr->spridx[1][ptr->size-1] * sizeof (uniform REAL));
    memcpy (p1, ptr->spring[1][1], ptr->spridx[1][ptr->size-1] * sizeof (uniform REAL));

    delete ptr->spring[1][0];
    delete ptr->spring[1][1];

    ptr->spring[1][0] = p0;
    ptr->spring[1][1] = p1;

    ptr->spring_size[1] = n;
  }

  if (ptr->dashpot_size < ptr->dashidx[ptr->size-1] + dashpot_inc)
  {
    uniform int n = 2*(ptr->dashidx[ptr->size-1] + dashpot_inc);

    uniform REAL * uniform p0 = uniform new uniform REAL[n];
    uniform REAL * uniform p1 = uniform new uniform REAL[n];

    memcpy (p0, ptr->dashpot[0], ptr->dashidx[ptr->size-1] * sizeof (uniform REAL));
    memcpy (p1, ptr->dashpot[1], ptr->dashidx[ptr->size-1] * sizeof (uniform REAL));

    delete ptr->dashpot[0];
    delete ptr->dashpot[1];

    ptr->dashpot[0] = p0;
    ptr->dashpot[1] = p1;

    ptr->dashpot_size = n;
  }
}

/* free local data of master springs */
static void master_spring_free_data (uniform master_spring * uniform ptr)
{
  if (ptr->spring_size[0])
  {
    delete ptr->spring[0][0];
    delete ptr->spring[0][1];
  }

  if (ptr->spring_size[1])
  {
    delete ptr->spring[1][0];
    delete ptr->spring[1][1];
  }

  if (ptr->dashpot_size)
  {
    delete ptr->dashpot[0];
    delete ptr->dashpot[1];
  }
}

/* free global array of master springs */
export void master_spring_free (uniform master_spring * uniform spr, uniform int size)
{
  for (uniform int i = 0; i < size; i ++)
  {
    master_spring_free_data (&spr[i]);
    uniform master_spring * uniform ptr = spr[i].next;
    while (ptr)
    {
      master_spring_free_data (ptr);
      uniform master_spring * uniform next = ptr->next;
      delete ptr;
      ptr = next;
    }
  }

  delete spr;
}

/* allocate global array of slave springs */
export uniform slave_spring * uniform slave_spring_alloc (uniform slave_spring * uniform old, uniform int nold, uniform int size)
{
  uniform slave_spring * uniform spr = uniform new uniform slave_spring [size];

  if (nold)
  {
    memcpy (spr, old, nold * sizeof (uniform slave_spring));

    delete old;
  }

  for (uniform int i = 0; i < size; i ++)
  {
    spr[i].size = 0;
    spr[i].next = NULL;
    spr[i].lock = -1;
  }
 
  return spr;
}

/* free global array of slave springs */
export void slave_spring_free (uniform slave_spring * uniform spr, uniform int size)
{
  for (uniform int i = 0; i < size; i ++)
  {
    uniform slave_spring * uniform ptr = spr[i].next;
    while (ptr)
    {
      uniform slave_spring * uniform next = ptr->next;
      delete ptr;
      ptr = next;
    }
  }

  delete spr;
}

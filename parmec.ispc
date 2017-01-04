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

#include "constants.h"
#include "macros.h"

/* optimal number of hardware threads */
export uniform int ispc_num_cores ()
{
  return num_cores();
};

/* aligned real allocator */
export uniform REAL * uniform  aligned_real_alloc (uniform int n)
{
  uniform REAL * uniform ptr = uniform new uniform REAL [n];
  memset (ptr, 0, n * sizeof (uniform REAL));
  return ptr;
}

export void aligned_real_free (uniform REAL * uniform ptr)
{
  delete ptr;
}

/* aligned int allocator */
export uniform int * uniform  aligned_int_alloc (uniform int n)
{
  uniform int * uniform ptr = uniform new uniform int [n];
  memset (ptr, 0, n * sizeof(uniform int));
  return ptr;
}

export void aligned_int_free (uniform int * uniform ptr)
{
  delete ptr;
}

/* invert inertia properties */
task void invert (uniform int span, uniform int size,
                  uniform REAL * uniform inertia[9], uniform REAL * uniform inverse[9],
                  uniform REAL mass[], uniform REAL invm[])
{
  uniform int start = taskIndex*span;
  uniform int end = taskIndex == taskCount-1 ? size: start+span;

  foreach (i = start ... end)
  {
    REAL a[9], x[9], det;

    cfor (int j = 0; j < 9; j ++)
    {
      a[j] = inertia[j][i];
    }
    
    INVERT (a, x, det);

    cfor (int j = 0; j < 9; j ++)
    {
      inverse[j][i] = x[j];
    }

    invm[i] = 1.0 / mass[i];
  }
}

/* invert inertia properties */
export void invert_inertia (uniform int threads, uniform int size, uniform REAL * uniform inertia[9],
  uniform REAL * uniform inverse[9], uniform REAL mass[], uniform REAL invm[])
{
  uniform int num = threads;
  uniform int span = size / num;

  launch[num] invert (span, size, inertia, inverse, mass, invm);
}

/* estimate critical time step */
export uniform REAL critical (uniform int parnum, uniform REAL mass[], uniform int pairnum, uniform REAL * uniform iparam[NIPARAM],
  uniform int sprnum, uniform int * uniform sprpart[2], uniform REAL * uniform spring[2], uniform int spridx[],
  uniform REAL * uniform dashpot[2], uniform int dashidx[])
{
  uniform REAL kmax, emax, mmin, omax, step;
  uniform int j;
  REAL k, e, m;

  k = 0.0;
  e = 0.0;

  foreach (i = 0 ... pairnum) /* contact springs */
  {
    k = max(iparam[SPRING][i],k);
    e = max(iparam[DAMPER][i],e);
  }

  kmax = reduce_max(k);
  emax = reduce_max(e);

  k = 0.0;
  e = 0.0;

  for (j = 0; j < sprnum; j ++) /* spring constraints */
  {
    uniform REAL kspr = 0.0;

    for (uniform int i = spridx[j]; i < spridx[j+1]-1; i ++)
    {
      uniform REAL k = (spring[1][i+1]-spring[1][i]) / (spring[0][i+1]-spring[0][i]);

      kspr = max(abs(k), kspr);
    }

    kmax = max(kspr, kmax);

    uniform REAL espr = 0.0;

    for (uniform int i = dashidx[j]; i < dashidx[j+1]-1; i ++)
    {
      uniform REAL e = (dashpot[1][i+1]-dashpot[1][i]) / (dashpot[0][i+1]-dashpot[0][i]);

      espr = max(abs(e), espr);
    }

    uniform REAL mspr = 0.0;
    
    if (sprpart[1][j] == -1)
    {
      mspr = mass[sprpart[0][j]];
    }
    else
    {
      uniform REAL im0 = 1.0/mass[sprpart[0][j]], im1 = 1.0/mass[sprpart[1][j]];

      mspr = 1.0/(im0+im1);
    }

    uniform REAL omega = sqrt(kspr/mspr);

    uniform REAL ratio = espr / (2*mspr*omega);

    emax = max(ratio, emax);
  }

  m = 1E99;

  foreach (i = 0 ... parnum)
  {
    m = min(mass[i],m);
  }

  mmin = reduce_min(m);

  omax = sqrt (kmax/mmin);

  /* http://www.dynasupport.com/tutorial/ls-dyna-users-guide/time-integration */
  step  = (2.0/omax)*(sqrt(1.0+emax*emax) - emax);

  return step;
}

/* Euler task */
task void euler_task (uniform int span, uniform int parnum, uniform REAL * uniform angular[6],
  uniform REAL * uniform linear[3], uniform REAL * uniform rotation[9], uniform REAL * uniform position[3], uniform REAL step)
{
  uniform int start = taskIndex*span;
  uniform int end = taskIndex == taskCount-1 ? parnum : start+span;

  foreach (i = start ... end)
  {
    REAL O[3], L1[9], DL[9], L2[9], o[3];

    O[0] = step*angular[0][i];
    O[1] = step*angular[1][i];
    O[2] = step*angular[2][i];

    expmap (O[0], O[1], O[2], DL[0], DL[1], DL[2], DL[3], DL[4], DL[5], DL[6], DL[7], DL[8]);

    L1[0] = rotation[0][i];
    L1[1] = rotation[1][i];
    L1[2] = rotation[2][i];
    L1[3] = rotation[3][i];
    L1[4] = rotation[4][i];
    L1[5] = rotation[5][i];
    L1[6] = rotation[6][i];
    L1[7] = rotation[7][i];
    L1[8] = rotation[8][i];

    NNMUL (L1, DL, L2);

    rotation[0][i] = L2[0];
    rotation[1][i] = L2[1];
    rotation[2][i] = L2[2];
    rotation[3][i] = L2[3];
    rotation[4][i] = L2[4];
    rotation[5][i] = L2[5];
    rotation[6][i] = L2[6];
    rotation[7][i] = L2[7];
    rotation[8][i] = L2[8];

    NVMUL (L2, O, o);

    angular[3][i] = o[0];
    angular[4][i] = o[1];
    angular[5][i] = o[2];

    position[0][i] += step * linear[0][i];
    position[1][i] += step * linear[0][i];
    position[2][i] += step * linear[0][i];
  }
}

/* explicit Euler step */
export void euler (uniform int threads, uniform int parnum, uniform REAL * uniform angular[6],
  uniform REAL * uniform linear[3], uniform REAL * uniform rotation[9], uniform REAL * uniform position[3], uniform REAL step)
{
  launch [threads] euler_task (parnum/threads, parnum, angular, linear, rotation, position, step);
}

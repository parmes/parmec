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

/* optimal number of hardware ntasks */
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
export void invert_inertia (uniform int ntasks, uniform int size, uniform REAL * uniform inertia[9],
  uniform REAL * uniform inverse[9], uniform REAL mass[], uniform REAL invm[])
{
  uniform int num = ntasks;
  uniform int span = size / num;

  launch[num] invert (span, size, inertia, inverse, mass, invm);

  sync;
}

/* estimate critical time step */
export uniform REAL critical (uniform int parnum, uniform REAL mass[], uniform int pairnum, uniform REAL * uniform iparam[NIPARAM],
  uniform int sprnum, uniform int * uniform sprpart[2], uniform REAL * uniform spring[2], uniform int spridx[],
  uniform REAL * uniform dashpot[2], uniform int dashidx[])
{
  uniform REAL kmax, emax, mmin, omax, step;
  uniform int j, sppmax;
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

  uniform int * uniform spp = uniform new uniform int [parnum];
  foreach (i = 0 ... parnum) spp[i] = 0; /* springs per particle */
  for (j = 0; j < sprnum; j ++)
  {
    if (sprpart[1][j] == -1)
    {
      spp[sprpart[0][j]] ++;
    }
    else
    {
      spp[sprpart[0][j]] ++;
      spp[sprpart[1][j]] ++;
    }
  }
  for (sppmax = j = 0; j < parnum; j ++)
  {
    sppmax = max(spp[j], sppmax);
  }
  delete spp;

  /* TODO: include spring inertia and points and calculate projected
           point masses and more accurate critical steps per spring */

  for (j = 0; j < sprnum; j ++) /* spring constraints */
  {
    uniform REAL kspr = 0.0;

    for (uniform int i = spridx[j]; i < spridx[j+1]-1; i ++)
    {
      uniform REAL k = (spring[1][i+1]-spring[1][i]) / (spring[0][i+1]-spring[0][i]);

      kspr = max(abs(k), kspr);
    }

    kspr *= sppmax;

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

/* time step determination task */
task void time_step_task (uniform int span, uniform int parnum, uniform REAL mass[], uniform REAL * uniform inertia[9],
     uniform REAL * uniform kact[3], uniform REAL kmax[], uniform REAL emax[], uniform REAL * uniform krot[6], uniform REAL step[])
{
  uniform int start = taskIndex*span;
  uniform int end = taskIndex == taskCount-1 ? parnum : start+span;

  REAL step0 = 1E99;

  foreach (i = start ... end)
  {
    /* minimum step based on
     * translational free vibration;
     * eigenvalue problem is scalar: m * x'' + n * k * x = 0;
     * calculated eigenvalue is: omega = sqrt (n * k / m);
     * time step is: h = 2 / omega */

    REAL nx = kact[0][i];
    REAL ny = kact[1][i];
    REAL nz = kact[2][i];
    REAL nxy = max(nx, ny);
    REAL nxyz = max (nxy, nz);

    if (nxyz > 0.0)
    {
      REAL k = kmax[i];
      REAL mxyz = mass[i];
      REAL omegaxyz = sqrt(nxyz*k/mxyz);

      /* minimum step based on 
       * rotational free vibration;
       * genralized symmetric eigenvalue problem is: J * O'' + K * O = 0;
       * which can be rewritten as: K * x - omega^2 * J * x = 0;
       * or using lambda = omega^2: K * x - lambda * J * x = 0;
       * this is further transformed using Cholesky decomposition of J = L*L^T;
       * transformation:     K * x - lambda * L*L^T * x = 0;
       *                 L^-1 * K * L^-T * x - lambda * I * x = 0;
       *     same as:                  A * x - lambda * I * x = 0;
       *      with:                       A = L^-1 * K * L^-T;
       *
       * From http://rosettacode.org/wiki/Cholesky_decomposition
       * Cholesky factorization of a 3x3 matrix reads:
       *    [l11   0    0]
       * L =[l21 l22    0]
       *    [l31 l32  l33] , where:
       *
       * l11 = sqrt(J11)
       * l22 = sqrt(J22 - l21^2)
       * l33 = sqrt(J33 - (l31^2 + l32^2))
       * l21 = (1/l11) * J21;
       * l31 = (1/l11) * J31;
       * l32 = (1/l22) * (J32 - l31*l21);
       *
       * L can be inverted directly; Having calculate A we can use a direct algorithm from 
       * Wikipedia, https://en.wikipedia.org/wiki/Eigenvalue_algorithm, in order to calculate
       * all eigenvalues of a symmetric 3x3 problem; we then select the maximum one; */

      REAL K[9], L[9], inv, invL[9], temp[9], A[9], J[9];

      K[0] = k*krot[0][i];
      K[1] = k*krot[3][i];
      K[2] = k*krot[4][i];
      K[3] = K[1];
      K[4] = k*krot[1][i];
      K[5] = k*krot[5][i];
      K[6] = K[2];
      K[7] = K[5];
      K[8] = k*krot[2][i];

      J[0] = inertia[0][i];
      J[1] = inertia[1][i];
      J[2] = inertia[2][i];
      J[3] = inertia[3][i];
      J[4] = inertia[4][i];
      J[5] = inertia[5][i];
      J[6] = inertia[6][i];
      J[7] = inertia[7][i];
      J[8] = inertia[8][i];

      REAL mrot = extremum_eigenvalue(J, 1); /* minimum */

      L[0] = sqrt(J[0]);
      inv = 1.0/L[0];
      L[1] = inv * J[1];
      L[2] = inv * J[2];
      L[3] = 0.0;
      L[4] = sqrt(J[4] - L[1]*L[1]);
      inv = 1.0/L[4];
      L[5] = inv * (J[5] - L[2]*L[1]);
      L[6] = 0.0;
      L[7] = 0.0;
      L[8] = sqrt(J[8] - (L[2]*L[2] + L[5]*L[5]));

      INVERT (L, invL, inv);
      NTMUL (K, invL, temp);
      NNMUL (invL, temp, A);

      REAL lambda = extremum_eigenvalue(A, 0); /* maximum */
      REAL omegarot = sqrt(lambda);

      /* http://www.dynasupport.com/tutorial/ls-dyna-users-guide/time-integration */
      REAL damping = emax[i];
      REAL ratioxyz = damping / (2.0*mxyz*omegaxyz);
      REAL ratiorot = damping / (2.0*mrot*omegarot);
      REAL hxyz = (2.0/omegaxyz)*(sqrt(1.0+ratioxyz*ratioxyz) - ratioxyz);
      REAL hrot = (2.0/omegarot)*(sqrt(1.0+ratiorot*ratiorot) - ratiorot);
      REAL hcri = min (hxyz, hrot);
      step0 = min(hcri, step0);
    }
  }

  step[taskIndex] = reduce_min (step0);
}

/* explicit Euler step */
export void euler (uniform int ntasks, uniform int parnum, uniform REAL * uniform angular[6],
  uniform REAL * uniform linear[3], uniform REAL * uniform rotation[9], uniform REAL * uniform position[3], uniform REAL step)
{
  launch [ntasks] euler_task (parnum/ntasks, parnum, angular, linear, rotation, position, step);

  sync;
}

/* determine stable time step */
export uniform REAL determine_time_step (uniform int ntasks, uniform int parnum, uniform REAL mass[],
  uniform REAL * uniform inertia[9], uniform REAL * uniform kact[3], uniform REAL kmax[], uniform REAL emax[],
  uniform REAL * uniform krot[6], uniform REAL step_current, uniform REAL adaptive)
{
  if (adaptive > 0.0 && adaptive <= 1.0)
  {
    uniform REAL step_min = 1E99;

    uniform REAL * uniform step = uniform new uniform REAL [ntasks];

    launch [ntasks] time_step_task (parnum/ntasks, parnum, mass, inertia, kact, kmax, emax, krot, step);

    sync;

    for (uniform int i = 0; i < ntasks; i ++)
    {
      if (step[i] < step_min) step_min = step[i];
    }

    delete step;

    if (step_min < 1E99)
    {
      uniform REAL step_new = adaptive*step_min;

      if (step_new < 0.1*step_current) return 0.1*step_current; /* XXX */
      else if (step_new > 10*step_current) return 10*step_current; /* XXX */
      else return step_new;
    }
    else return step_current;
  }
  else return step_current;
}

#if 0
/* estimate critical time step */
export uniform REAL critical (uniform int ntasks, uniform int parnum, uniform REAL * uniform position[6],
  uniform REAL mass[], uniform REAL * uniform inertia[9], uniform REAL * uniform kact[3], uniform REAL kmax[],
  uniform REAL emax[], uniform REAL * uniform krot[6], uniform int pairnum,  uniform REAL * uniform iparam[NIPARAM],
  uniform int sprnum, uniform int * uniform sprpart[2], uniform REAL * uniform spring[2], uniform int spridx[],
  uniform REAL * uniform dashpot[2], uniform int dashidx[], uniform REAL * uniform sprpnt[2][6])
{
}
#endif

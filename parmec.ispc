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
#include "condet.h"

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

/* superseeded by critical_perparticle; does not include torsion springs; can this be made relevant again? */
#if 0
/* estimate critical time step */
export uniform REAL critical (uniform int parnum, uniform REAL mass[], uniform int pairnum, uniform REAL * uniform iparam[NIPARAM],
    uniform int sprnum, uniform int * uniform sprpart[2], uniform REAL * uniform spring[2], uniform int spridx[],
    uniform REAL * uniform dashpot[2], uniform int dashidx[])
{
  uniform REAL kmax, emax, mmin, omax;
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

    if (dashidx[j]+1 == dashidx[j+1]) /* damping ratio given directly */
    {
      uniform REAL ratio = dashpot[0][dashidx[j]];

      emax = max(ratio, emax);
    }
    else /* damping ratio from damping curve */
    {
      uniform REAL espr = 0.0;

      for (uniform int i = dashidx[j]; i < dashidx[j+1]-1; i ++)
      {
        uniform REAL e = (dashpot[1][i+1]-dashpot[1][i]) / (dashpot[0][i+1]-dashpot[0][i]);

        espr = max(abs(e), espr);
      }

      uniform REAL omega = sqrt(kspr/mspr);

      if (omega > 0.0)
      {
        uniform REAL ratio = espr / (2*mspr*omega);

        emax = max(ratio, emax);
      }
    }
  }

  m = 1E99;

  foreach (i = 0 ... parnum)
  {
    m = min(mass[i],m);
  }

  mmin = reduce_min(m);

  omax = sqrt (kmax/mmin);

  if (omax > 0.0)
  {
    /* http://www.dynasupport.com/tutorial/ls-dyna-users-guide/time-integration */
    return (2.0/omax)*(sqrt(1.0+emax*emax) - emax);
  }
  else return 1.0;
}
#endif

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
    position[1][i] += step * linear[1][i];
    position[2][i] += step * linear[2][i];
  }
}

/* calculate per-particle cirtical time step estimates hcri[], circular frequency ocri[], and damping ratio rcri[]  */
task void time_step_task (uniform int span, uniform int parnum, uniform REAL mass[], uniform REAL * uniform inertia[9],
    uniform REAL * uniform kact[3], uniform REAL kmax[], uniform REAL emax[], uniform REAL * uniform krot[6],
    uniform REAL hcri[], uniform REAL ocri[], uniform REAL rcri[])
{
  uniform int start = taskIndex*span;
  uniform int end = taskIndex == taskCount-1 ? parnum : start+span;

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
      hcri[i] = min (hxyz, hrot);
      cif (ocri) ocri[i] = max (omegaxyz, omegarot);
      cif (rcri) rcri[i] = max (ratioxyz, ratiorot);
    }
    else hcri[i] = REAL_MAX;
  }
}

/* explicit Euler step */
export void euler (uniform int ntasks, uniform int parnum, uniform REAL * uniform angular[6],
    uniform REAL * uniform linear[3], uniform REAL * uniform rotation[9], uniform REAL * uniform position[3], uniform REAL step)
{
  launch [ntasks] euler_task (parnum/ntasks, parnum, angular, linear, rotation, position, step);

  sync;
}

/* calculate adaptive time step */
export uniform REAL adaptive_timestep (uniform int ntasks, uniform int parnum, uniform REAL mass[],
    uniform REAL * uniform inertia[9], uniform REAL * uniform kact[3], uniform REAL kmax[],
    uniform REAL emax[], uniform REAL * uniform krot[6], uniform REAL curstep, uniform REAL adaptive)
{
  uniform REAL * uniform hcri = uniform new uniform REAL [parnum];

  launch [ntasks] time_step_task (parnum/ntasks, parnum, mass, inertia, kact, kmax, emax, krot, hcri, NULL, NULL);

  sync;

  REAL hmin = REAL_MAX;

  foreach (i = 0 ... parnum)
  {
    hmin = min(hcri[i], hmin);
  }

  uniform REAL step = reduce_min (hmin);

  delete hcri;

  if (step < REAL_MAX) return adaptive*step; /* TODO --> a more robust strategy ? */
  else return curstep;
}

/* calculate equivalent relative point mass */
export uniform REAL eqm (uniform int part1, uniform REAL * uniform point1, uniform int part2, uniform REAL * uniform point2,
    uniform REAL * uniform direction, uniform REAL * uniform inverse[9], uniform REAL invm[], uniform REAL * uniform position[3])
{
  uniform REAL Jiv[9], im, A[3], Ask[9], tmp[9], Miv[9], mass;

  /* XXX: include rotation multiplication like in joints.cpp */

  Jiv[0] = inverse[0][part1];
  Jiv[1] = inverse[1][part1];
  Jiv[2] = inverse[2][part1];
  Jiv[3] = inverse[3][part1];
  Jiv[4] = inverse[4][part1];
  Jiv[5] = inverse[5][part1];
  Jiv[6] = inverse[6][part1];
  Jiv[7] = inverse[7][part1];
  Jiv[8] = inverse[8][part1];
  im = invm[part1];
  A[0] = position[3][part1] - point1[0];
  A[1] = position[4][part1] - point1[1];
  A[2] = position[5][part1] - point1[2];
  VECSKEW (A, Ask);
  TNMUL(Ask, Jiv, tmp);
  NNMUL (tmp, Ask, Miv);
  Miv[0] += im;
  Miv[4] += im;
  Miv[8] += im;

  if (part2 >= 0)
  {
    Jiv[0] = inverse[0][part2];
    Jiv[1] = inverse[1][part2];
    Jiv[2] = inverse[2][part2];
    Jiv[3] = inverse[3][part2];
    Jiv[4] = inverse[4][part2];
    Jiv[5] = inverse[5][part2];
    Jiv[6] = inverse[6][part2];
    Jiv[7] = inverse[7][part2];
    Jiv[8] = inverse[8][part2];
    im = invm[part2];
    A[0] = position[3][part2] - point2[0];
    A[1] = position[4][part2] - point2[1];
    A[2] = position[5][part2] - point2[2];
    VECSKEW (A, Ask);
    TNMUL(Ask, Jiv, tmp);
    NNMUL (tmp, Ask, Jiv);
    NNADD (Miv, Jiv, Miv);
    Miv[0] += im;
    Miv[4] += im;
    Miv[8] += im;
  }

  if (direction && DOT(direction,direction) > 0.0)
  {
    NVMUL (Miv, direction, tmp);
    mass = 1.0/DOT(tmp, direction);
  }
  else mass = 1.0/extremum_eigenvalue (Miv, 1);

  return mass;
}

/* calculate per-spring cirtical time step estimates hcri[], circular frequency ocri[], and damping ratio rcri[]  */
export void critical_perspring (uniform int sprnum, uniform int * uniform sprpart[2], uniform REAL * uniform sprpnt[2][6],
    uniform REAL * uniform sprdir[6], uniform REAL * uniform spring[2], uniform int spridx[], uniform REAL * uniform dashpot[2],
    uniform int dashidx[], uniform REAL * uniform inverse[9], uniform REAL invm[], uniform REAL * uniform position[3],
    uniform REAL hcri[], uniform REAL ocri[], uniform REAL rcri[])
{
  for (uniform int j = 0; j < sprnum; j ++)
  {
    uniform int i = sprpart[0][j], k = sprpart[1][j];
    uniform REAL pi[3], pk[3], d[3];

    pi[0] = sprpnt[0][3][j];
    pi[1] = sprpnt[0][4][j];
    pi[2] = sprpnt[0][5][j];

    pk[0] = sprpnt[1][3][j];
    pk[1] = sprpnt[1][4][j];
    pk[2] = sprpnt[1][5][j];

    d[0] = sprdir[0][j];
    d[1] = sprdir[1][j];
    d[2] = sprdir[2][j];

    uniform REAL kspr = 0.0;

    for (uniform int i = spridx[j]; i < spridx[j+1]-1; i ++)
    {
      uniform REAL k = (spring[1][i+1]-spring[1][i]) / (spring[0][i+1]-spring[0][i]);

      kspr = max(abs(k), kspr);
    }

    uniform REAL mspr = eqm (i, pi, k, pk, d, inverse, invm, position);

    uniform REAL omega = sqrt(kspr/mspr), ksi = 0.0;

    if (dashidx[j]+1 == dashidx[j+1]) /* damping ratio given directly */
    {
      ksi = dashpot[0][dashidx[j]];
    }
    else /* damping ratio from damping curve */
    {
      uniform REAL espr = 0.0;

      for (uniform int i = dashidx[j]; i < dashidx[j+1]-1; i ++)
      {
        uniform REAL e = (dashpot[1][i+1]-dashpot[1][i]) / (dashpot[0][i+1]-dashpot[0][i]);

        espr = max(abs(e), espr);
      }

      if (omega > 0.0) ksi = espr / (2*mspr*omega);
    }

    /* http://www.dynasupport.com/tutorial/ls-dyna-users-guide/time-integration */
    if (omega > 0.0) hcri[j] = (2.0/omega)*(sqrt(1.0+ksi*ksi) - ksi);
    else hcri[j] = 1.0;
    ocri[j] = omega;
    rcri[j] = ksi;
  }
}

/* ciritical per-particle step contacts task */
task void critical_perparticle_contacts (uniform int span, uniform master_conpnt master[], uniform slave_conpnt slave[],
    uniform int parnum, uniform REAL * uniform rotation[9], uniform REAL * uniform position[6],
    uniform REAL * uniform kact[3], uniform REAL kmax[], uniform REAL emax[], uniform REAL * uniform krot[6])
{
  uniform int start = taskIndex*span;
  uniform int end = taskIndex == taskCount-1 ? parnum: start+span;

  for (uniform int i = start; i < end; i ++)
  {
    REAL A[3], f[3], a[3], dot, krot0[6], kact0[3], kmax0, kcur, emax0, ecur;
    uniform REAL L[9], po[3];

    L[0] = rotation[0][i];
    L[1] = rotation[1][i];
    L[2] = rotation[2][i];
    L[3] = rotation[3][i];
    L[4] = rotation[4][i];
    L[5] = rotation[5][i];
    L[6] = rotation[6][i];
    L[7] = rotation[7][i];
    L[8] = rotation[8][i];

    po[0] = position[0][i];
    po[1] = position[1][i];
    po[2] = position[2][i];

    SET6 (krot0, 0.0);
    SET(kact0, 0.0);
    kmax0 = 0.0;
    emax0 = 0.0;

    for (uniform master_conpnt * uniform m = &master[i]; m; m = m->next)
    {
      foreach (j = 0 ... m->size)
      {
        a[0] = m->point[0][j]-po[0];
        a[1] = m->point[1][j]-po[1];
        a[2] = m->point[2][j]-po[2];

        kact0[0] += abs(m->normal[0][j]);
        kact0[1] += abs(m->normal[1][j]);
        kact0[2] += abs(m->normal[2][j]);
        kcur = m->kcur[j];
        kmax0 = max(kcur, kmax0);
        ecur = m->ecur[j];
        emax0 = max(ecur, emax0);

        NVMUL (L, a, A);
        dot = DOT(A, A);
        krot0[0] += dot - A[0]*A[0];
        krot0[1] += dot - A[1]*A[1];
        krot0[2] += dot - A[2]*A[2];
        krot0[3] += -A[0]*A[1];
        krot0[4] += -A[0]*A[2];
        krot0[5] += -A[1]*A[2];
      }
    }

    for (uniform slave_conpnt * uniform s = &slave[i]; s; s = s->next)
    {
      foreach (j = 0 ... s->size)
      {
        a[0] = s->point[0][j]-po[0];
        a[1] = s->point[1][j]-po[1];
        a[2] = s->point[2][j]-po[2];

        f[0] = s->force[0][j];
        f[1] = s->force[1][j];
        f[2] = s->force[2][j];
        REAL len = LEN(f);
        cif (len != 0.0)
        {
          len = 1.0/len;
          SCALE (f,len);
        }
        else SET (f, 0.0);
        kact0[0] += abs(f[0]);
        kact0[1] += abs(f[1]);
        kact0[2] += abs(f[2]);
        kcur = s->kcur[j];
        kmax0 = max(kcur, kmax0);
        ecur = s->ecur[j];
        emax0 = max(ecur, emax0);

        NVMUL (L, a, A);
        dot = DOT(A, A);
        krot0[0] += dot - A[0]*A[0];
        krot0[1] += dot - A[1]*A[1];
        krot0[2] += dot - A[2]*A[2];
        krot0[3] += -A[0]*A[1];
        krot0[4] += -A[0]*A[2];
        krot0[5] += -A[1]*A[2];
      }
    }

    krot[0][i] = reduce_add (krot0[0]);
    krot[1][i] = reduce_add (krot0[1]);
    krot[2][i] = reduce_add (krot0[2]);
    krot[3][i] = reduce_add (krot0[3]);
    krot[4][i] = reduce_add (krot0[4]);
    krot[5][i] = reduce_add (krot0[5]);

    kact[0][i] = reduce_add (kact0[0]);
    kact[1][i] = reduce_add (kact0[1]);
    kact[2][i] = reduce_add (kact0[2]);

    kmax[i] = reduce_max (kmax0);
    emax[i] = reduce_max (emax0);
  }
}

/* ciritical per-particle step springs task */
task void critical_perparticle_springs (uniform int span, uniform REAL * uniform rotation[9], uniform REAL * uniform position[6],
    uniform REAL * uniform inverse[9], uniform REAL invm[], uniform int sprnum, uniform int sprflg[], uniform int * uniform sprpart[2],
    uniform REAL * uniform sprpnt[2][6], uniform REAL * uniform sprdir[6], uniform REAL * uniform spring[2], uniform int spridx[],
    uniform REAL * uniform dashpot[2], uniform int dashidx[], uniform REAL * uniform kact[3], uniform REAL kmax[], uniform REAL emax[],
    uniform REAL * uniform krot[6], uniform int lock[])
{
  uniform int start = taskIndex*span;
  uniform int end = taskIndex == taskCount-1 ? sprnum: start+span;

  for (uniform int j = start; j < end; j ++)
  {
    uniform int i = sprpart[0][j], k = sprpart[1][j];
    uniform REAL pi[3], pk[3], d[3], a[3], A[3], dot;
    uniform REAL Li[9], xi[3], Lk[9], xk[3];

    Li[0] = rotation[0][i];
    Li[1] = rotation[1][i];
    Li[2] = rotation[2][i];
    Li[3] = rotation[3][i];
    Li[4] = rotation[4][i];
    Li[5] = rotation[5][i];
    Li[6] = rotation[6][i];
    Li[7] = rotation[7][i];
    Li[8] = rotation[8][i];

    xi[0] = position[0][i];
    xi[1] = position[1][i];
    xi[2] = position[2][i];

    if (k >= 0)
    {
      Lk[0] = rotation[0][k];
      Lk[1] = rotation[1][k];
      Lk[2] = rotation[2][k];
      Lk[3] = rotation[3][k];
      Lk[4] = rotation[4][k];
      Lk[5] = rotation[5][k];
      Lk[6] = rotation[6][k];
      Lk[7] = rotation[7][k];
      Lk[8] = rotation[8][k];

      xk[0] = position[0][k];
      xk[1] = position[1][k];
      xk[2] = position[2][k];
    }

    pi[0] = sprpnt[0][0][j];
    pi[1] = sprpnt[0][1][j];
    pi[2] = sprpnt[0][2][j];

    pk[0] = sprpnt[1][0][j];
    pk[1] = sprpnt[1][1][j];
    pk[2] = sprpnt[1][2][j];

    d[0] = sprdir[0][j];
    d[1] = sprdir[1][j];
    d[2] = sprdir[2][j];

    uniform REAL kspr = 0.0;

    for (uniform int i = spridx[j]; i < spridx[j+1]-1; i ++)
    {
      uniform REAL k = (spring[1][i+1]-spring[1][i]) / (spring[0][i+1]-spring[0][i]);

      kspr = max(abs(k), kspr);
    }

    uniform REAL espr = 0.0;

    if (dashidx[i]+1 == dashidx[i+1]) /* critical damping ratio */
    {
      uniform REAL Pi[3], Pk[3];

      Pi[0] = sprpnt[0][3][j];
      Pi[1] = sprpnt[0][4][j];
      Pi[2] = sprpnt[0][5][j];

      Pk[0] = sprpnt[1][3][j];
      Pk[1] = sprpnt[1][4][j];
      Pk[2] = sprpnt[1][5][j];

      uniform REAL mass = eqm (i, Pi, k, Pk, d, inverse, invm, position);

      uniform REAL ratio = dashpot[0][dashidx[i]];

      espr = ratio * 2.0 * sqrt (kspr * mass);
    }
    else for (uniform int i = dashidx[j]; i < dashidx[j+1]-1; i ++) /* lookup table */
    {
      uniform REAL e = (dashpot[1][i+1]-dashpot[1][i]) / (dashpot[0][i+1]-dashpot[0][i]);

      espr = max(abs(e), espr);
    }

    /* update particle i data */
    while (lock[i] != taskIndex)
      atomic_compare_exchange_global (&lock[i], -1, taskIndex); /* lock access i data */

    if (DOT (d,d) > 0.0)
    {
      if (sprflg[i] & SPRDIR_PLANAR)
      {
        uniform REAL t1[3], t2[3];
        COPY (d, t1);
        t1[0] += 0.123;
        t1[1] += 0.346;
        t1[2] += 0.789;
        PRODUCT (d, t1, t2);
        uniform REAL l2 = LEN(t2);
        uniform REAL i2 = 1.0/l2;
        SCALE (t2, i2);
        PRODUCT (t2, d, t1);
        kact[0][i] += abs(t1[0]) + abs(t2[0]);
        kact[1][i] += abs(t1[1]) + abs(t2[1]);
        kact[2][i] += abs(t1[2]) + abs(t2[2]);
      }
      else
      {
        kact[0][i] += abs(d[0]);
        kact[1][i] += abs(d[1]);
        kact[2][i] += abs(d[2]);
      }
    } 
    else
    {
      kact[0][i] += 1.0;
      kact[1][i] += 1.0;
      kact[2][i] += 1.0;
    }

    SUB (pi, xi, a);
    TVMUL (Li, a, A);
    dot = DOT(A,A);
    krot[0][i] += dot-A[0]*A[0]; /* K(1,1) */
    krot[1][i] += dot-A[1]*A[1]; /* K(2,2) */
    krot[2][i] += dot-A[2]*A[2]; /* K(3,3) */
    krot[3][i] += -A[0]*A[1]; /* K(2,1) */
    krot[4][i] += -A[0]*A[2]; /* K(3,1) */
    krot[5][i] += -A[1]*A[2]; /* K(3,3) --> see LS-DYNA's Theory Manual Sec 28.1 (rev. >= 7709) */

    kmax[i] = max(kmax[i], kspr);
    emax[i] = max(emax[i], espr);

    lock[i] = -1; /* unlock */

    if (k >= 0)
    {
      while (lock[k] != taskIndex)
        atomic_compare_exchange_global (&lock[k], -1, taskIndex); /* lock access k data */

      if (DOT (d,d) > 0.0)
      {
        if (sprflg[k] & SPRDIR_PLANAR)
        {
          uniform REAL t1[3], t2[3];
          COPY (d, t1);
          t1[0] += 0.123;
          t1[1] += 0.346;
          t1[2] += 0.789;
          PRODUCT (d, t1, t2);
          uniform REAL l2 = LEN(t2);
          uniform REAL i2 = 1.0/l2;
          SCALE (t2, i2);
          PRODUCT (t2, d, t1);
          kact[0][k] += abs(t1[0]) + abs(t2[0]);
          kact[1][k] += abs(t1[1]) + abs(t2[1]);
          kact[2][k] += abs(t1[2]) + abs(t2[2]);
        }
        else
        {
          kact[0][k] += abs(d[0]);
          kact[1][k] += abs(d[1]);
          kact[2][k] += abs(d[2]);
        }
      } 
      else
      {
        kact[0][k] += 1.0;
        kact[1][k] += 1.0;
        kact[2][k] += 1.0;
      }

      SUB (pk, xk, a);
      TVMUL (Lk, a, A);
      dot = DOT(A,A);
      krot[0][k] += dot-A[0]*A[0]; /* K(1,1) */
      krot[1][k] += dot-A[1]*A[1]; /* K(2,2) */
      krot[2][k] += dot-A[2]*A[2]; /* K(3,3) */
      krot[3][k] += -A[0]*A[1]; /* K(2,1) */
      krot[4][k] += -A[0]*A[2]; /* K(3,1) */
      krot[5][k] += -A[1]*A[2]; /* K(3,3) --> see LS-DYNA's Theory Manual Sec 28.1 (rev. >= 7709) */

      kmax[k] = max(kmax[k], kspr);
      emax[k] = max(emax[k], espr);

      lock[k] = -1; /* unlock */
    }
  }
}

/* calculate per-particle cirtical time step estimates hcri[], circular frequency ocri[], and damping ratio rcri[]  */
export uniform REAL critical_perparticle (uniform int ntasks,
    uniform master_conpnt master[], uniform slave_conpnt slave[],
    uniform int parnum, uniform REAL mass[], uniform REAL * uniform inertia[9],
    uniform REAL invm[], uniform REAL * uniform inverse[9], uniform REAL * uniform rotation[9],
    uniform REAL * uniform position[6], uniform int sprnum, uniform int sprflg[], uniform int * uniform sprpart[2],
    uniform REAL * uniform sprpnt[2][6], uniform REAL * uniform sprdir[6], uniform REAL * uniform spring[2],
    uniform int spridx[], uniform REAL * uniform dashpot[2], uniform int dashidx[], uniform REAL * uniform kact[3],
    uniform REAL kmax[], uniform REAL emax[], uniform REAL * uniform krot[6], uniform REAL hcri[],
    uniform REAL ocri[], uniform REAL rcri[])
{
  /* process contacts and initialise krot, kact, kmax, emax */
  launch [ntasks] critical_perparticle_contacts (parnum/ntasks, master, slave, parnum,
      rotation, position, kact, kmax, emax, krot);
  sync;

  /* process springs and update krot, kact, kmax, emax */
  uniform int * uniform lock = uniform new uniform int [parnum];
  foreach (i = 0 ... parnum) lock[i] = -1;
  launch [ntasks] critical_perparticle_springs (sprnum/ntasks, rotation, position, inverse,
      invm, sprnum, sprflg, sprpart, sprpnt, sprdir, spring, spridx,
      dashpot, dashidx, kact, kmax, emax, krot, lock);
  sync;
  delete lock;

  /* per-particle critical time step task */
  launch [ntasks] time_step_task (parnum/ntasks, parnum, mass, inertia, kact, kmax, emax, krot, hcri, ocri, rcri);

  sync;
}

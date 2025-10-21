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

/* dynamics task */
task void dynamics_task (uniform int span, uniform master_conpnt master[], uniform slave_conpnt slave[],
    uniform int parnum, uniform REAL * uniform angular[6], uniform REAL * uniform linear[3],
    uniform REAL * uniform rotation[9], uniform REAL * uniform position[6],
    uniform REAL * uniform inertia[9], uniform REAL * uniform inverse[9],
    uniform REAL mass[], uniform REAL invm[], uniform REAL damping[6],
    uniform REAL * uniform force[3], uniform REAL * uniform torque[3],
    uniform int flags[], uniform REAL step0, uniform REAL step1)
{
  uniform int start = taskIndex*span;
  uniform int end = taskIndex == taskCount-1 ? parnum: start+span;
  uniform REAL half = 0.5*step0;
  uniform REAL step = 0.5*(step0+step1);

  foreach (i = start ... end) /* time integration */
  {
    REAL O[3], o[3], v[3], L1[9], J[9], I[9], ma, im, f[3], t[3], T[3], DL[9], L2[9], A[3], B[3];

    if (flags[i] & SKIP) continue;

    O[0] = angular[0][i];
    O[1] = angular[1][i];
    O[2] = angular[2][i];

    v[0] = linear[0][i];
    v[1] = linear[1][i];
    v[2] = linear[2][i];

    L1[0] = rotation[0][i];
    L1[1] = rotation[1][i];
    L1[2] = rotation[2][i];
    L1[3] = rotation[3][i];
    L1[4] = rotation[4][i];
    L1[5] = rotation[5][i];
    L1[6] = rotation[6][i];
    L1[7] = rotation[7][i];
    L1[8] = rotation[8][i];

    J[0] = inertia[0][i];
    J[1] = inertia[1][i];
    J[2] = inertia[2][i];
    J[3] = inertia[3][i];
    J[4] = inertia[4][i];
    J[5] = inertia[5][i];
    J[6] = inertia[6][i];
    J[7] = inertia[7][i];
    J[8] = inertia[8][i];

    I[0] = inverse[0][i];
    I[1] = inverse[1][i];
    I[2] = inverse[2][i];
    I[3] = inverse[3][i];
    I[4] = inverse[4][i];
    I[5] = inverse[5][i];
    I[6] = inverse[6][i];
    I[7] = inverse[7][i];
    I[8] = inverse[8][i];

    im = invm[i];

    cif (damping[0] != 0.0 || damping[1] != 0.0 || damping[2] != 0.0)
    {
      ma = mass[i];

      force[0][i] -= ma * damping[0] * v[0];
      force[1][i] -= ma * damping[1] * v[1];
      force[2][i] -= ma * damping[2] * v[2];
    }
    cif (damping[3] != 0.0 || damping[4] != 0.0 || damping[5] != 0.0)
    {
      o[0] = damping[3]*angular[3][i];
      o[1] = damping[4]*angular[4][i];
      o[2] = damping[5]*angular[5][i];

      TVMUL (L1, o, t);
      NVMUL (J, t, T);
      NVMUL (L1, T, t);

      torque[0][i] -= t[0];
      torque[1][i] -= t[1];
      torque[2][i] -= t[2];
    }

    f[0] = force[0][i];
    f[1] = force[1][i];
    f[2] = force[2][i];

    t[0] = torque[0][i];
    t[1] = torque[1][i];
    t[2] = torque[2][i];

    TVMUL (L1, t, T);

    expmap (-half*O[0], -half*O[1], -half*O[2], DL[0], DL[1], DL[2], DL[3], DL[4], DL[5], DL[6], DL[7], DL[8]);

    NVMUL (J, O, A);
    NVMUL (DL, A, B);
    ADDMUL (B, half, T, B);
    NVMUL (I, B, A); /* O(t+h/2) */

    NVMUL (J, A, B);
    PRODUCTSUB (A, B, T); /* T - O(t+h/2) x J O(t+h/2) */

    SCALE (T, step0);
    NVADDMUL (O, I, T, O); /* O(t+h) */

    im *= step0;
    ADDMUL (v, im, f, v); /* v(t+h) */

    expmap (step*O[0], step*O[1], step*O[2], DL[0], DL[1], DL[2], DL[3], DL[4], DL[5], DL[6], DL[7], DL[8]);

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

    position[0][i] += step*v[0];
    position[1][i] += step*v[1];
    position[2][i] += step*v[2];

    angular[0][i] = O[0];
    angular[1][i] = O[1];
    angular[2][i] = O[2];
    angular[3][i] = o[0];
    angular[4][i] = o[1];
    angular[5][i] = o[2];

    linear[0][i] = v[0];
    linear[1][i] = v[1];
    linear[2][i] = v[2];
  }
}

/* update dynamics */
export void dynamics (uniform int ntasks, uniform master_conpnt master[], uniform slave_conpnt slave[],
    uniform int parnum, uniform REAL * uniform angular[6], uniform REAL * uniform linear[3],
    uniform REAL * uniform rotation[9], uniform REAL * uniform position[6],
    uniform REAL * uniform inertia[9], uniform REAL * uniform inverse[9],
    uniform REAL mass[], uniform REAL invm[], uniform REAL damping[6],
    uniform REAL * uniform force[3], uniform REAL * uniform torque[3],
    uniform int flags[], uniform REAL step0, uniform REAL step1)
{
  launch [ntasks] dynamics_task (parnum/ntasks, master, slave, parnum, angular, linear, rotation,
      position, inertia, inverse, mass, invm, damping, force, torque, flags, step0, step1);

  sync;
}

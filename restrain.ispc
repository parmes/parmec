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

/* restrain particle motion via velocities task */
task void restrain_velocities_task (uniform int span, uniform int rstnum, uniform int rstpart[],
    uniform REAL * uniform rstlin[9], uniform REAL * uniform rstang[9],
    uniform REAL * uniform linear[3], uniform REAL * uniform angular[6], uniform REAL * uniform rotation[9])
{
  uniform int start = taskIndex*span;
  uniform int end = taskIndex == taskCount-1 ? rstnum: start+span;

  foreach (i = start ... end)
  {
    int j = rstpart[i];
    REAL v[3] = {linear[0][j], linear[1][j], linear[2][j]};
    REAL O[3] = {angular[0][j], angular[1][j], angular[2][j]};
    REAL o[3] = {angular[3][j], angular[4][j], angular[5][j]};
    REAL L[9] = {rotation[0][j], rotation[1][j], rotation[2][j],
      rotation[3][j], rotation[4][j], rotation[5][j],
      rotation[6][j], rotation[7][j], rotation[8][j]};
    REAL l0[3] = {rstlin[0][i], rstlin[1][i], rstlin[2][i]};
    REAL l1[3] = {rstlin[3][i], rstlin[4][i], rstlin[5][i]};
    REAL l2[3] = {rstlin[6][i], rstlin[7][i], rstlin[8][i]};
    REAL a0[3] = {rstang[0][i], rstang[1][i], rstang[2][i]};
    REAL a1[3] = {rstang[3][i], rstang[4][i], rstang[5][i]};
    REAL a2[3] = {rstang[6][i], rstang[7][i], rstang[8][i]};
    REAL A0[3], A1[3], A2[3];
    REAL dot;

    dot = DOT(v,l0);
    SUBMUL (v, dot, l0, v);
    dot = DOT(v,l1);
    SUBMUL (v, dot, l1, v);
    dot = DOT(v,l2);
    SUBMUL (v, dot, l2, v);

    dot = DOT(o,a0);
    SUBMUL (o, dot, a0, o);
    dot = DOT(o,a1);
    SUBMUL (o, dot, a1, o);
    dot = DOT(o,a2);
    SUBMUL (o, dot, a2, o);

    TVMUL (L, a0, A0);
    TVMUL (L, a1, A1);
    TVMUL (L, a2, A2);

    dot = DOT(O,A0);
    SUBMUL (O, dot, A0, O);
    dot = DOT(O,A1);
    SUBMUL (O, dot, A1, O);
    dot = DOT(O,A2);
    SUBMUL (O, dot, A2, O);

    linear[0][j] = v[0];
    linear[1][j] = v[1];
    linear[2][j] = v[2];
    angular[0][j] = O[0];
    angular[1][j] = O[1];
    angular[2][j] = O[2];
    angular[3][j] = o[0];
    angular[4][j] = o[1];
    angular[5][j] = o[2];
  }
}

/* restrain particle motion via forces task */
task void restrain_forces_task (uniform int span, uniform int rstnum, uniform int rstpart[],
    uniform REAL * uniform rstlin[9], uniform REAL * uniform rstang[9],
    uniform REAL * uniform force[3], uniform REAL * uniform torque[3])
{
  uniform int start = taskIndex*span;
  uniform int end = taskIndex == taskCount-1 ? rstnum: start+span;

  foreach (i = start ... end)
  {
    int j = rstpart[i];
    REAL f[3] = {force[0][j], force[1][j], force[2][j]};
    REAL t[3] = {torque[0][j], torque[1][j], torque[2][j]};
    REAL l0[3] = {rstlin[0][i], rstlin[1][i], rstlin[2][i]};
    REAL l1[3] = {rstlin[3][i], rstlin[4][i], rstlin[5][i]};
    REAL l2[3] = {rstlin[6][i], rstlin[7][i], rstlin[8][i]};
    REAL a0[3] = {rstang[0][i], rstang[1][i], rstang[2][i]};
    REAL a1[3] = {rstang[3][i], rstang[4][i], rstang[5][i]};
    REAL a2[3] = {rstang[6][i], rstang[7][i], rstang[8][i]};
    REAL dot;

    dot = DOT(f,l0);
    SUBMUL (f, dot, l0, f);
    dot = DOT(f,l1);
    SUBMUL (f, dot, l1, f);
    dot = DOT(f,l2);
    SUBMUL (f, dot, l2, f);

    dot = DOT(t,a0);
    SUBMUL (t, dot, a0, t);
    dot = DOT(t,a1);
    SUBMUL (t, dot, a1, t);
    dot = DOT(t,a2);
    SUBMUL (t, dot, a2, t);

    force[0][j] = f[0];
    force[1][j] = f[1];
    force[2][j] = f[2];
    torque[0][j] = t[0];
    torque[1][j] = t[1];
    torque[2][j] = t[2];
  }
}

/* restrain particle motion by adjusting velocities */
export void restrain_velocities (uniform int ntasks, uniform int rstnum, uniform int rstpart[],
    uniform REAL * uniform rstlin[9], uniform REAL * uniform rstang[9],
    uniform REAL * uniform linear[3], uniform REAL * uniform angular[6], uniform REAL * uniform rotation[9])
{
  launch [ntasks] restrain_velocities_task (rstnum/ntasks, rstnum, rstpart, rstlin, rstang, linear, angular, rotation);

  sync;
}

/* restrain particle motion by adjusting forces */
export void restrain_forces (uniform int ntasks, uniform int rstnum, uniform int rstpart[],
    uniform REAL * uniform rstlin[9], uniform REAL * uniform rstang[9],
    uniform REAL * uniform force[3], uniform REAL * uniform torque[3])
{
  launch [ntasks] restrain_forces_task (rstnum/ntasks, rstnum, rstpart, rstlin, rstang, force, torque);

  sync;
}

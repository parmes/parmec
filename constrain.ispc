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

/* constrain particle motion task */
task void constrain_task (uniform int span, uniform int cnsnum, uniform int cnspart[],
  uniform REAL * uniform cnslin[9], uniform REAL * uniform cnsang[9],
  uniform REAL * uniform force[3], uniform REAL * uniform torque[3])
{
  uniform int start = taskIndex*span;
  uniform int end = taskIndex == taskCount-1 ? cnsnum: start+span;

  foreach (i = start ... end)
  {
    int j = cnspart[i];
    REAL f[3] = {force[0][j], force[1][j], force[2][j]};
    REAL t[3] = {torque[0][j], torque[1][j], torque[2][j]};
    REAL l0[3] = {cnslin[0][i], cnslin[1][i], cnslin[2][i]};
    REAL l1[3] = {cnslin[3][i], cnslin[4][i], cnslin[5][i]};
    REAL l2[3] = {cnslin[6][i], cnslin[7][i], cnslin[8][i]};
    REAL a0[3] = {cnsang[0][i], cnsang[1][i], cnsang[2][i]};
    REAL a1[3] = {cnsang[3][i], cnsang[4][i], cnsang[5][i]};
    REAL a2[3] = {cnsang[6][i], cnsang[7][i], cnsang[8][i]};
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

/* constrain particle motion */
export void constrain (uniform int threads, uniform int cnsnum, uniform int cnspart[],
  uniform REAL * uniform cnslin[9], uniform REAL * uniform cnsang[9],
  uniform REAL * uniform force[3], uniform REAL * uniform torque[3])
{
  launch [threads] constrain_task (cnsnum/threads, cnsnum, cnspart, cnslin, cnsang, force, torque);
}

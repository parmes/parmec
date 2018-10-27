/*
The MIT License (MIT)

Copyright (c) 2018 Tomasz Koziara

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

#ifndef __joints__
#define __joints__

#ifdef __cplusplus
namespace parmec { /* namespace */
#endif

/* reset joints symbolic information after change */
void reset_joints_symbolic (int jnum, int *jpart[2]);

/* solve joints and update forces */
void solve_joints (int jnum, int *jpart[2], REAL *jpoint[3], REAL *jreac[3], int parnum,
  REAL *position[6], REAL *rotation[9], REAL *inertia[9], REAL *inverse[9], REAL mass[], REAL invm[],
  REAL damping[6], REAL *linear[3], REAL *angular[6], REAL *force[6], REAL *torque[6], REAL step0, REAL step1);

#ifdef __cplusplus
} /* namespace */
#endif
#endif

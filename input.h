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

#ifndef __input__
#define __input__

namespace parmec
{
/* update obstacles time histories from callbacks */
void obstaclev (int obsnum, REAL *obsang, REAL *obslin, pointer_t anghis[], pointer_t linhis[], REAL time);

/* prescribe particle velocity */
void prescribe_velocity (int prsnum, pointer_t tms[], int prspart[], pointer_t prslin[], int *tmslin[3], int linkind[],
  pointer_t prsang[], int *tmsang[3], int angkind[], REAL time, REAL *rotation[9], REAL *linear[3], REAL *angular[6]);

/* prescribe particle acceleration */
void prescribe_acceleration (int prsnum, pointer_t tms[], int prspart[], pointer_t prslin[], int *tmslin[3], int linkind[],
  pointer_t prsang[], int *tmsang[3], int angkind[], REAL time, REAL mass[], REAL *inertia[9], REAL *force[3], REAL *torque[3]);

/* read gravity and global damping */
void read_gravity_and_damping (REAL time, pointer_t *tms, pointer_t gravfunc[3],
  int gravtms[3], REAL gravity[3], pointer_t lindamp, pointer_t angdamp, REAL damping[6]);

/* call interval callback */
REAL current_interval (pointer_t func, REAL time);
}

#endif

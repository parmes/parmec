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

#ifndef __constants__
#define __constants__

#ifdef __cplusplus
namespace parmec
{
#endif

enum {DENSITY = 0, YOUNG, POISSON, NMAT}; /* bulk material constants */

enum {GRANULAR_FORCE, USER_FORCE}; /* interaction force kind */

enum {SPRING = 0, DAMPER, FRISTAT, FRIDYN, FRIROL, FRIDRIL, KSKN, NIPARAM}; /* surface material constants */

enum {ANALYTICAL = 1, OUTREST = 2}; /* particle flags */

enum {HIS_LIST = 1, HIS_SPHERE = 2, HIS_BOX = 4, HIS_POINT = 8}; /* history kind flags */

enum {HIS_PX, HIS_PY, HIS_PZ, HIS_PL, HIS_DX, HIS_DY, HIS_DZ, HIS_DL,
      HIS_VX, HIS_VY, HIS_VZ, HIS_VL, HIS_OX, HIS_OY, HIS_OZ, HIS_OL,
      HIS_FX, HIS_FY, HIS_FZ, HIS_FL, HIS_TX, HIS_TY, HIS_TZ, HIS_TL, 
      HIS_STROKE, HIS_STF, HIS_SF, HIS_TIME}; /* history entities */

enum {OUT_FORMAT_VTK = 1, OUT_FORMAT_XDMF = 2}; /* output format */

enum {OUT_MODE_SPH = 1, OUT_MODE_MESH = 2, OUT_MODE_RB = 4, OUT_MODE_CD = 8, OUT_MODE_SD = 16}; /* output modes */

enum {OUT_NUMBER = 1, OUT_COLOR = 2, OUT_DISPL = 4, OUT_ORIENT = 8, OUT_LINVEL = 16, OUT_ANGVEL = 32, OUT_FORCE = 64,
      OUT_TORQUE = 128, OUT_F = 256, OUT_FN = 512, OUT_FT = 1024, OUT_SF = 2048, OUT_AREA = 4096, OUT_PAIR = 8192}; /* output entities */

enum {SPRING_NONLINEAR_ELASTIC, SPRING_GENERAL_NONLINEAR}; /* spring type */

enum {SPRING_YIELDED = 1, SPRING_UNLOADING = 2, SPRDIR_CONSTANT = 4, SPRDIR_PLANAR = 8, SPRDIR_FOLLOWER = 16}; /* spring flags */

#ifdef __cplusplus
}
#endif

#endif

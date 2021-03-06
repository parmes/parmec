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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <algorithm>
#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <omp.h>
#include "timeseries.h"
#include "macros.h"
#include "parmec.h"
#include "timer.h"
#include "mesh.h"
#include "input.h"
#include "output.h"
#include "joints.h"
#include "constants.h"
#include "parmec_ispc.h"
#include "partition_ispc.h"
#include "forces_ispc.h"
#include "dynamics_ispc.h"
#include "shapes_ispc.h"
#include "obstacles_ispc.h"
#include "restrain_ispc.h"

using namespace parmec;
using namespace ispc; /* ISPC calls are used below */

#ifdef __cplusplus
namespace parmec { /* namespace */
#endif

char **argv; /* input arguments */
int argc; /* input arguments count */

int ntasks; /* number of tasks */

char *output_path; /* output path */
int output_frame; /* output files frame */

int stepnum; /* current step number */
REAL curtime; /* current time */
REAL curstep; /* current step */
REAL curtime_output; /* current output time */
REAL curtime_history; /* current history time */

int matnum; /* number of materials */
REAL *mparam[NMAT]; /* material parameters */
int material_buffer_size; /* size of the buffer */

int pairnum; /* number of pairings */
int *pairs; /* color pairs */
int *ikind; /* interaction kind */
REAL *iparam[NIPARAM]; /* interaction parameters */
pointer_t *uforce; /* user force callbacks */
int pair_buffer_size; /* size of the buffer */

int ellnum; /* number of ellipsoids */
int ellcon; /* index of the first ellipsoid used in contact detection */
int *ellcol; /* ellipsoid color */
int *part; /* ellipsoid to particle map */
REAL *center[6]; /* ellipsoid current and reference centers */
REAL *radii[3]; /* ellipsoid radii; radii[1][] < 0 indicates sphere */
REAL *orient[18]; /* ellipsoid current and reference orientations */
int ellipsoid_buffer_size; /* size of the buffer */

int parnum; /* particle count */
int *parmat; /* particle material */
REAL *angular[6]; /* angular velocities (referential, spatial) */
REAL *linear[3]; /* linear velocities */
REAL *rotation[9]; /* rotation operators */
REAL *position[6]; /* mass center current and reference positions */
REAL *inertia[9]; /* inertia tensors */
REAL *inverse[9]; /* inverse inertia tensors */
REAL *mass; /* scalar mass */
REAL *invm; /* inverse scalar mass */
REAL *force[3]; /* total spatial force */
REAL *torque[3]; /* total spatial torque */
REAL *kact[3]; /* time step control --> number of active constraints per particle per dimension */
REAL *kmax; /* time step control --> maximum stiffness coefficient per particle */
REAL *emax; /* time step control --> maximum damper coefficient per particle */
REAL *krot[6]; /* time step control --> symmetric rotational unit stiffness matrix per particle */
int *flags; /* particle flags */
ispc::master_conpnt *master; /* master contact points */
ispc::slave_conpnt *slave; /* slave contact points */
int particle_buffer_size; /* size of the buffer */

int trinum; /* number of triangles */
int tricon; /* index of the first triangle used in contact detection */
int *tricol; /* triangle color */
int *triobs; /* triangle obstacle --> <0 - moving obstalce (-index-2), -1 - static obstacle, >= 0 - triangulated particle */
REAL *tri[3][3]; /* triangle vertices */
int triangle_buffer_size; /* size of the buffer */

int nodnum; /* number of nodes */
REAL *nodes[6]; /* current and reference nodes */
int *nodpart; /* node particle index */
int elenum; /* number of elements */
int *eletype; /* element type (4, 5, 6, 8) */
int *elenod; /* element nodes */
int *eleidx; /* element nodes start index */
int *elepart; /* element particle index */
int *elemat; /* element material index */
int facnum; /* number of faces (triangulated) */
int faccon; /* index of the first face used in contact detection */
int *facnod[3]; /* face nodes */
int *facpart; /* face particle index */
int *factri; /* face to triangle mapping */
int node_buffer_size; /* size of the nodes buffer */
int element_node_buffer_size; /* size of the element nodes buffer */
int element_buffer_size; /* size of the element buffers */
int face_buffer_size; /* size of the face buffer */

int obsnum; /* number of obstacles */
int *trirng; /* triangles range */
REAL *obspnt; /* obstacle spatial points */
REAL *obsang; /* obstacle angular velocities at t and t+h */
REAL *obslin; /* obstacle linear velocities at t and t+h */
pointer_t *anghis; /* angular velocity history */
pointer_t *linhis; /* linear velocity history */
int obstacle_buffer_size; /* size of the buffer */

int sprnum; /* number of spring constraints */
int *sprid; /* spring id --> number returned to user */
int *sprmap; /* map of spring ids to spring indices */
int *sprtype; /* spring type */
int *unspring; /* unspring action; -3 unused, -2 inactive, -1 zero force, >= 0 use lcurve */
int *sprpart[2]; /* spring constraint particle numbers */
REAL *sprpnt[2][6]; /* spring constraint current and reference points */
REAL *spring[2]; /* spring force lookup tables */
int *spridx; /* spring force lookup start index */
REAL *dashpot[2]; /* dashpot force lookup tables */
int *dashidx; /* dashpot force lookup start index */
REAL *unload[2]; /* spring unloading lookup tables */
int *unidx; /* spring unloading lookup start index */
REAL *yield[2]; /* spring yield limits: 0 compression and 1 tension */
REAL *sprdir[6]; /* spring direction (0:3 current and 3:6 reference) */
int *sprflg; /* spring flags */
int *sproffset; /* spring stroke offest load curve number */
REAL *sprfric; /* spring tangential friction coefficient */
REAL *sprkskn; /* ratio of normal to tangential spring and dashpot parameters */
REAL *sprsdsp[3]; /* tangential frictional displacement vector */
REAL *stroke0; /* initial spring stroke */
REAL *stroke[3]; /* current stroke: 0 current, 1 total compression, 2 total tension */
REAL *sprfrc[3]; /* total, spring, friction force magnitudes */
int springs_changed; /* spring input data changed flag */
int spring_buffer_size; /* size of the spring constraint buffer */
int spring_lookup_size; /* size of the spring force lookup tables */
int dashpot_lookup_size; /* size of the dashpot force lookup tables */
int unload_lookup_size; /* size of the unload force lookup tables */

int trqsprnum; /* number of torsion springs */
int *trqsprid; /* torsion spring id --> number returned to user */
int *trqsprmap; /* map of torsion spring ids to spring indices */
int *trqsprpart[2]; /* torsion spring particle indices */
REAL *trqzdir0[3]; /* torsion spring reference z direction */
REAL *trqxdir0[3]; /* torsion spring reference x direction */
REAL *krpy[3][2]; /* kroll, kpitch, kyaw spring angle-torque lookup tables */
int *krpyidx[3]; /* spring torque lookup start indexes */
REAL *drpy[3][2]; /* droll, dpitch, dyaw dashpot ang. velocity-torque lookup tables */
int *drpyidx[3]; /* dashpot torque lookup start indexes */
int *trqcone; /* cone constraint flags */
REAL *trqzdir1[3]; /* output: torsion spring z current direction */
REAL *trqxdir1[3]; /* output: torsion spring x current direction */
REAL *trqrpy[3]; /* output: torsion spring angles roll, pitch, yaw */
REAL *trqrpytot[3]; /* output: total moments conjugate with spring angles */
REAL *trqrpyspr[3]; /* output: spring moments wihout damper components */
REAL *trqrefpnt[3]; /* input: reference point coordinates on part[0] or (inf, inf, inf) */
int trqspr_changed; /* torqion spring input changed flag */
int trqspr_buffer_size; /* size of torsion spring constraint buffer */
int krpy_lookup_size[3]; /* size of spring angle-torque lookup tables */
int drpy_lookup_size[3]; /* size of spring ang. velocity-torque lookup tables */

int unsprnum; /* number of unspring definitions */
int *tsprings; /* test springs */
int *tspridx; /* test springs index range */
int *msprings; /* modified springs */
int *mspridx; /* modified springs index range */
REAL *unlim[2];  /* entity limits */
int *unent; /* entity type */
int *unop; /* test springs operator */
int *unabs; /* absolute value flag */
int *nsteps; /* number of steps between checks */
int *nfreq; /* number of nsteps for which tsprings exceed limits before msprings are modified */
int *unaction; /* unloading action; -1: instantaneous unloading, >=0: use lcurve */
int *activate; /* spring activation lists */
int *actidx; /* spring activation index range */
int unspring_buffer_size; /* size of unspring buffer */
int tsprings_buffer_size; /* size of tsprings buffer */
int msprings_buffer_size; /* size of msprings buffer */
int activate_buffer_size; /* size of activate buffer */

int rstnum; /* number of restraint */
int *rstpart; /* restrained particle numbers */
REAL *rstlin[9]; /* restrained linear directions */
REAL *rstang[9]; /* restrained angular directions */
int restrain_buffer_size; /* size of restrained particles buffer */

int jnum; /* number of joints */
int *jpart[2]; /* joint particles */
REAL *jpoint[3]; /* joint points */
REAL *jreac[3]; /* joint reactions */
int joints_changed; /* joints changed flag */
int joints_buffer_size; /* size of joints buffer */

int tmsnum; /* number of time series */
pointer_t *tms; /* time series */
int time_series_buffer_size; /* size of time series buffer */

int lcnum; /* number of load curves */
REAL *lcurve[2]; /* load curve data */
int *lcidx; /* load curve start index */
int lcurve_buffer_size; /* size of load curves buffer */
int lcurve_data_size; /* size of load curves data buffer */

int prsnum; /* number of particles with prescribed motion */
int *prspart; /* prescribed motion particle numbers */
pointer_t *prslin; /* prescribed linear motion time history callbacks */
int *tmslin[3]; /* prescribed linear motion time series */
int *linkind; /* prescribied linear motion signal kind: 0-velocity, 1-acceleration */
pointer_t *prsang; /* prescribed angular motion time history callbacks */
int *tmsang[3]; /* prescribed angular motion time series */
int *angkind; /* prescribied angular motion signal kind: 0-velocity, 1-acceleration */
int prescribe_buffer_size; /* size of prescribed particle motion buffer */

int hisnum; /* number of time histories */
int *hislst; /* history source lists */
int *hisidx; /* history source list start index */
int *hisent; /* history entity */
int *hiskind; /* history kind */
REAL *source[6]; /* source sphere or box definition or optional point */
pointer_t *h5file; /* optional paths to existing .h5 result files */
pointer_t *history; /* Python list storing history */
int history_buffer_size; /* size of history buffer */
int history_list_size; /* size of history particle lists buffer */

int outnum; /* number of output lists */
int *outmode; /* output mode */
int *outpart; /* output particle lists */
int *outidx; /* output particle list start index */
int *outent; /* output entities per output mode */
int outrest[2]; /* 0: default output entities for unlisted particles and, 1: default output mode */
int outformat; /* output format */
int output_buffer_size; /* size of output buffer */
int output_list_size; /* size of output particle lists buffer */

REAL gravity[3]; /* gravity vector */
pointer_t gravfunc[3]; /* gravity callbacks */
int gravtms[3]; /* gravity time series */

REAL damping[6]; /* linear and angular damping */
pointer_t lindamp; /* linead damping callback */
int lindamptms[3]; /* linear damping time series */
pointer_t angdamp; /* angular damping callback */
int angdamptms[3]; /* angular damping time series */

MAP *prescribed_body_forces; /* particle index based map of prescibed body forces */

/* grow integer buffer */
void integer_buffer_grow (int* &src, int num, int size)
{
  int *dst;

  ERRMEM (dst = aligned_int_alloc (size));
  memcpy (dst, src, sizeof (int)*num);
  aligned_int_free (src);
  src = dst;
}

/* grow real buffer */
void real_buffer_grow (REAL* &src, int num, int size)
{
  REAL *dst;

  ERRMEM (dst = aligned_real_alloc (size));
  memcpy (dst, src, sizeof (REAL)*num);
  aligned_real_free (src);
  src = dst;
}

/* grow pointer buffer */
void pointer_buffer_grow (pointer_t* &src, int num, int size)
{
  pointer_t *dst;

  ERRMEM (dst = new pointer_t [size]);
  memcpy (dst, src, sizeof (pointer_t)*num);
  delete [] src;
  src = dst;
}

/* initialize material buffer */
void material_buffer_init ()
{
  material_buffer_size = 64;

  for (int i = 0; i < NMAT; i ++)
  {
    mparam[i] = aligned_real_alloc (material_buffer_size);
  }

  matnum = 0;
}

/* grow material buffer */
void material_buffer_grow ()
{
  material_buffer_size *= 2;

  for (int i = 0; i < NMAT; i ++)
  {
    real_buffer_grow (mparam[i], matnum, material_buffer_size);
  }
}

/* initialize pairing buffer */
void pair_buffer_init ()
{
  pair_buffer_size = 64;

  pairs = aligned_int_alloc (2*pair_buffer_size);
  ikind = aligned_int_alloc (pair_buffer_size);
  for (int i = 0; i < NIPARAM; i ++)
  {
    iparam[i] = aligned_real_alloc (pair_buffer_size);
  }
  uforce = new pointer_t [pair_buffer_size];
}

/* grow pair buffer */
int pair_buffer_grow ()
{
  pair_buffer_size *= 2;

  integer_buffer_grow (pairs, pairnum, 2*pair_buffer_size);
  integer_buffer_grow (ikind, pairnum, pair_buffer_size);
  for (int i = 0; i < NIPARAM; i ++)
  {
    real_buffer_grow (iparam[i], pairnum, pair_buffer_size);
  }
  pointer_buffer_grow (uforce, pairnum, pair_buffer_size);

  return pair_buffer_size;
}

/* reset to default pair */
void pair_reset ()
{
  pairnum = 1; /* default */
  pairs[0] = 0;
  pairs[1] = 0;
  ikind[0] = GRANULAR_FORCE;
  iparam[SPRING][0] = 0.0;
  iparam[DAMPER][0] = 1.0;
  iparam[FRISTAT][0] = 0.0;
  iparam[FRIDYN][0] = 0.0;
  iparam[FRIROL][0] = 0.0;
  iparam[FRIDRIL][0] = 0.0;
  iparam[KSKN][0] = 0.5;
}

/* initialize ellipsoid buffer */
void ellipsoid_buffer_init ()
{
  ellipsoid_buffer_size = 256;

  part = aligned_int_alloc (ellipsoid_buffer_size);
  ellcol = aligned_int_alloc (ellipsoid_buffer_size);
  center[0] = aligned_real_alloc (ellipsoid_buffer_size);
  center[1] = aligned_real_alloc (ellipsoid_buffer_size);
  center[2] = aligned_real_alloc (ellipsoid_buffer_size);
  center[3] = aligned_real_alloc (ellipsoid_buffer_size);
  center[4] = aligned_real_alloc (ellipsoid_buffer_size);
  center[5] = aligned_real_alloc (ellipsoid_buffer_size);
  radii[0] = aligned_real_alloc (ellipsoid_buffer_size);
  radii[1] = aligned_real_alloc (ellipsoid_buffer_size);
  radii[2] = aligned_real_alloc (ellipsoid_buffer_size);
  orient[0] = aligned_real_alloc (ellipsoid_buffer_size);
  orient[1] = aligned_real_alloc (ellipsoid_buffer_size);
  orient[2] = aligned_real_alloc (ellipsoid_buffer_size);
  orient[3] = aligned_real_alloc (ellipsoid_buffer_size);
  orient[4] = aligned_real_alloc (ellipsoid_buffer_size);
  orient[5] = aligned_real_alloc (ellipsoid_buffer_size);
  orient[6] = aligned_real_alloc (ellipsoid_buffer_size);
  orient[7] = aligned_real_alloc (ellipsoid_buffer_size);
  orient[8] = aligned_real_alloc (ellipsoid_buffer_size);
  orient[9] = aligned_real_alloc (ellipsoid_buffer_size);
  orient[10] = aligned_real_alloc (ellipsoid_buffer_size);
  orient[11] = aligned_real_alloc (ellipsoid_buffer_size);
  orient[12] = aligned_real_alloc (ellipsoid_buffer_size);
  orient[13] = aligned_real_alloc (ellipsoid_buffer_size);
  orient[14] = aligned_real_alloc (ellipsoid_buffer_size);
  orient[15] = aligned_real_alloc (ellipsoid_buffer_size);
  orient[16] = aligned_real_alloc (ellipsoid_buffer_size);
  orient[17] = aligned_real_alloc (ellipsoid_buffer_size);
}

/* grow ellipsoid buffer */
int ellipsoid_buffer_grow ()
{
  ellipsoid_buffer_size *= 2;

  integer_buffer_grow (part, ellnum, ellipsoid_buffer_size);
  integer_buffer_grow (ellcol, ellnum, ellipsoid_buffer_size);
  real_buffer_grow (center[0], ellnum, ellipsoid_buffer_size);
  real_buffer_grow (center[1], ellnum, ellipsoid_buffer_size);
  real_buffer_grow (center[2], ellnum, ellipsoid_buffer_size);
  real_buffer_grow (center[3], ellnum, ellipsoid_buffer_size);
  real_buffer_grow (center[4], ellnum, ellipsoid_buffer_size);
  real_buffer_grow (center[5], ellnum, ellipsoid_buffer_size);
  real_buffer_grow (radii[0], ellnum, ellipsoid_buffer_size);
  real_buffer_grow (radii[1], ellnum, ellipsoid_buffer_size);
  real_buffer_grow (radii[2], ellnum, ellipsoid_buffer_size);
  real_buffer_grow (orient[0], ellnum, ellipsoid_buffer_size);
  real_buffer_grow (orient[1], ellnum, ellipsoid_buffer_size);
  real_buffer_grow (orient[2], ellnum, ellipsoid_buffer_size);
  real_buffer_grow (orient[3], ellnum, ellipsoid_buffer_size);
  real_buffer_grow (orient[4], ellnum, ellipsoid_buffer_size);
  real_buffer_grow (orient[5], ellnum, ellipsoid_buffer_size);
  real_buffer_grow (orient[6], ellnum, ellipsoid_buffer_size);
  real_buffer_grow (orient[7], ellnum, ellipsoid_buffer_size);
  real_buffer_grow (orient[8], ellnum, ellipsoid_buffer_size);
  real_buffer_grow (orient[9], ellnum, ellipsoid_buffer_size);
  real_buffer_grow (orient[10], ellnum, ellipsoid_buffer_size);
  real_buffer_grow (orient[11], ellnum, ellipsoid_buffer_size);
  real_buffer_grow (orient[12], ellnum, ellipsoid_buffer_size);
  real_buffer_grow (orient[13], ellnum, ellipsoid_buffer_size);
  real_buffer_grow (orient[14], ellnum, ellipsoid_buffer_size);
  real_buffer_grow (orient[15], ellnum, ellipsoid_buffer_size);
  real_buffer_grow (orient[16], ellnum, ellipsoid_buffer_size);
  real_buffer_grow (orient[17], ellnum, ellipsoid_buffer_size);

  return ellipsoid_buffer_size;
}

/* init particle buffer */
void particle_buffer_init ()
{
  particle_buffer_size = 256;

  parmat = aligned_int_alloc (particle_buffer_size);
  angular[0] = aligned_real_alloc (particle_buffer_size);
  angular[1] = aligned_real_alloc (particle_buffer_size);
  angular[2] = aligned_real_alloc (particle_buffer_size);
  angular[3] = aligned_real_alloc (particle_buffer_size);
  angular[4] = aligned_real_alloc (particle_buffer_size);
  angular[5] = aligned_real_alloc (particle_buffer_size);
  linear[0] = aligned_real_alloc (particle_buffer_size);
  linear[1] = aligned_real_alloc (particle_buffer_size);
  linear[2] = aligned_real_alloc (particle_buffer_size);
  rotation[0] = aligned_real_alloc (particle_buffer_size);
  rotation[1] = aligned_real_alloc (particle_buffer_size);
  rotation[2] = aligned_real_alloc (particle_buffer_size);
  rotation[3] = aligned_real_alloc (particle_buffer_size);
  rotation[4] = aligned_real_alloc (particle_buffer_size);
  rotation[5] = aligned_real_alloc (particle_buffer_size);
  rotation[6] = aligned_real_alloc (particle_buffer_size);
  rotation[7] = aligned_real_alloc (particle_buffer_size);
  rotation[8] = aligned_real_alloc (particle_buffer_size);
  position[0] = aligned_real_alloc (particle_buffer_size);
  position[1] = aligned_real_alloc (particle_buffer_size);
  position[2] = aligned_real_alloc (particle_buffer_size);
  position[3] = aligned_real_alloc (particle_buffer_size);
  position[4] = aligned_real_alloc (particle_buffer_size);
  position[5] = aligned_real_alloc (particle_buffer_size);
  inertia[0] = aligned_real_alloc (particle_buffer_size);
  inertia[1] = aligned_real_alloc (particle_buffer_size);
  inertia[2] = aligned_real_alloc (particle_buffer_size);
  inertia[3] = aligned_real_alloc (particle_buffer_size);
  inertia[4] = aligned_real_alloc (particle_buffer_size);
  inertia[5] = aligned_real_alloc (particle_buffer_size);
  inertia[6] = aligned_real_alloc (particle_buffer_size);
  inertia[7] = aligned_real_alloc (particle_buffer_size);
  inertia[8] = aligned_real_alloc (particle_buffer_size);
  inverse[0] = aligned_real_alloc (particle_buffer_size);
  inverse[1] = aligned_real_alloc (particle_buffer_size);
  inverse[2] = aligned_real_alloc (particle_buffer_size);
  inverse[3] = aligned_real_alloc (particle_buffer_size);
  inverse[4] = aligned_real_alloc (particle_buffer_size);
  inverse[5] = aligned_real_alloc (particle_buffer_size);
  inverse[6] = aligned_real_alloc (particle_buffer_size);
  inverse[7] = aligned_real_alloc (particle_buffer_size);
  inverse[8] = aligned_real_alloc (particle_buffer_size);
  mass = aligned_real_alloc (particle_buffer_size);
  invm = aligned_real_alloc (particle_buffer_size);
  force[0] = aligned_real_alloc (particle_buffer_size);
  force[1] = aligned_real_alloc (particle_buffer_size);
  force[2] = aligned_real_alloc (particle_buffer_size);
  torque[0] = aligned_real_alloc (particle_buffer_size);
  torque[1] = aligned_real_alloc (particle_buffer_size);
  torque[2] = aligned_real_alloc (particle_buffer_size);
  kact[0] = aligned_real_alloc (particle_buffer_size);
  kact[1] = aligned_real_alloc (particle_buffer_size);
  kact[2] = aligned_real_alloc (particle_buffer_size);
  kmax = aligned_real_alloc (particle_buffer_size);
  emax = aligned_real_alloc (particle_buffer_size);
  krot[0] = aligned_real_alloc (particle_buffer_size);
  krot[1] = aligned_real_alloc (particle_buffer_size);
  krot[2] = aligned_real_alloc (particle_buffer_size);
  krot[3] = aligned_real_alloc (particle_buffer_size);
  krot[4] = aligned_real_alloc (particle_buffer_size);
  krot[5] = aligned_real_alloc (particle_buffer_size);
  flags = aligned_int_alloc (particle_buffer_size);
  master = master_alloc (NULL, 0, particle_buffer_size);
  slave = slave_alloc (NULL, 0, particle_buffer_size);

  parnum = 0;
}

/* grow particle buffer */
int particle_buffer_grow ()
{
  particle_buffer_size *= 2;

  integer_buffer_grow (parmat, parnum, particle_buffer_size);
  real_buffer_grow (angular[0], parnum, particle_buffer_size);
  real_buffer_grow (angular[1], parnum, particle_buffer_size);
  real_buffer_grow (angular[2], parnum, particle_buffer_size);
  real_buffer_grow (angular[3], parnum, particle_buffer_size);
  real_buffer_grow (angular[4], parnum, particle_buffer_size);
  real_buffer_grow (angular[5], parnum, particle_buffer_size);
  real_buffer_grow (linear[0], parnum, particle_buffer_size);
  real_buffer_grow (linear[1], parnum, particle_buffer_size);
  real_buffer_grow (linear[2], parnum, particle_buffer_size);
  real_buffer_grow (rotation[0], parnum, particle_buffer_size);
  real_buffer_grow (rotation[1], parnum, particle_buffer_size);
  real_buffer_grow (rotation[2], parnum, particle_buffer_size);
  real_buffer_grow (rotation[3], parnum, particle_buffer_size);
  real_buffer_grow (rotation[4], parnum, particle_buffer_size);
  real_buffer_grow (rotation[5], parnum, particle_buffer_size);
  real_buffer_grow (rotation[6], parnum, particle_buffer_size);
  real_buffer_grow (rotation[7], parnum, particle_buffer_size);
  real_buffer_grow (rotation[8], parnum, particle_buffer_size);
  real_buffer_grow (position[0], parnum, particle_buffer_size);
  real_buffer_grow (position[1], parnum, particle_buffer_size);
  real_buffer_grow (position[2], parnum, particle_buffer_size);
  real_buffer_grow (position[3], parnum, particle_buffer_size);
  real_buffer_grow (position[4], parnum, particle_buffer_size);
  real_buffer_grow (position[5], parnum, particle_buffer_size);
  real_buffer_grow (inertia[0], parnum, particle_buffer_size);
  real_buffer_grow (inertia[1], parnum, particle_buffer_size);
  real_buffer_grow (inertia[2], parnum, particle_buffer_size);
  real_buffer_grow (inertia[3], parnum, particle_buffer_size);
  real_buffer_grow (inertia[4], parnum, particle_buffer_size);
  real_buffer_grow (inertia[5], parnum, particle_buffer_size);
  real_buffer_grow (inertia[6], parnum, particle_buffer_size);
  real_buffer_grow (inertia[7], parnum, particle_buffer_size);
  real_buffer_grow (inertia[8], parnum, particle_buffer_size);
  real_buffer_grow (inverse[0], parnum, particle_buffer_size);
  real_buffer_grow (inverse[1], parnum, particle_buffer_size);
  real_buffer_grow (inverse[2], parnum, particle_buffer_size);
  real_buffer_grow (inverse[3], parnum, particle_buffer_size);
  real_buffer_grow (inverse[4], parnum, particle_buffer_size);
  real_buffer_grow (inverse[5], parnum, particle_buffer_size);
  real_buffer_grow (inverse[6], parnum, particle_buffer_size);
  real_buffer_grow (inverse[7], parnum, particle_buffer_size);
  real_buffer_grow (inverse[8], parnum, particle_buffer_size);
  real_buffer_grow (mass, parnum, particle_buffer_size);
  real_buffer_grow (invm, parnum, particle_buffer_size);
  real_buffer_grow (force[0], parnum, particle_buffer_size);
  real_buffer_grow (force[1], parnum, particle_buffer_size);
  real_buffer_grow (force[2], parnum, particle_buffer_size);
  real_buffer_grow (torque[0], parnum, particle_buffer_size);
  real_buffer_grow (torque[1], parnum, particle_buffer_size);
  real_buffer_grow (torque[2], parnum, particle_buffer_size);
  real_buffer_grow (kact[0], parnum, particle_buffer_size);
  real_buffer_grow (kact[1], parnum, particle_buffer_size);
  real_buffer_grow (kact[2], parnum, particle_buffer_size);
  real_buffer_grow (kmax, parnum, particle_buffer_size);
  real_buffer_grow (emax, parnum, particle_buffer_size);
  real_buffer_grow (krot[0], parnum, particle_buffer_size);
  real_buffer_grow (krot[1], parnum, particle_buffer_size);
  real_buffer_grow (krot[2], parnum, particle_buffer_size);
  real_buffer_grow (krot[3], parnum, particle_buffer_size);
  real_buffer_grow (krot[4], parnum, particle_buffer_size);
  real_buffer_grow (krot[5], parnum, particle_buffer_size);
  integer_buffer_grow (flags, parnum, particle_buffer_size);
  master = master_alloc (master, parnum, particle_buffer_size);
  slave = slave_alloc (slave, parnum, particle_buffer_size);

  return particle_buffer_size;
}

/* init triangle buffer */
void triangle_buffer_init ()
{
  triangle_buffer_size = 256;

  tricol = aligned_int_alloc (triangle_buffer_size);
  triobs = aligned_int_alloc (triangle_buffer_size);
  tri[0][0] = aligned_real_alloc (triangle_buffer_size);
  tri[0][1] = aligned_real_alloc (triangle_buffer_size);
  tri[0][2] = aligned_real_alloc (triangle_buffer_size);
  tri[1][0] = aligned_real_alloc (triangle_buffer_size);
  tri[1][1] = aligned_real_alloc (triangle_buffer_size);
  tri[1][2] = aligned_real_alloc (triangle_buffer_size);
  tri[2][0] = aligned_real_alloc (triangle_buffer_size);
  tri[2][1] = aligned_real_alloc (triangle_buffer_size);
  tri[2][2] = aligned_real_alloc (triangle_buffer_size);

  trinum = 0;
  tricon = 0;
}

/* grow triangle buffer */
int triangle_buffer_grow ()
{
  triangle_buffer_size *= 2;

  integer_buffer_grow (tricol, trinum, triangle_buffer_size);
  integer_buffer_grow (triobs, trinum, triangle_buffer_size);
  real_buffer_grow (tri[0][0], trinum, triangle_buffer_size);
  real_buffer_grow (tri[0][1], trinum, triangle_buffer_size);
  real_buffer_grow (tri[0][2], trinum, triangle_buffer_size);
  real_buffer_grow (tri[1][0], trinum, triangle_buffer_size);
  real_buffer_grow (tri[1][1], trinum, triangle_buffer_size);
  real_buffer_grow (tri[1][2], trinum, triangle_buffer_size);
  real_buffer_grow (tri[2][0], trinum, triangle_buffer_size);
  real_buffer_grow (tri[2][1], trinum, triangle_buffer_size);
  real_buffer_grow (tri[2][2], trinum, triangle_buffer_size);

  return triangle_buffer_size;
}

/* init element buffer */
void element_buffer_init ()
{
  node_buffer_size = 256;
  element_node_buffer_size = 1024;
  element_buffer_size = 256;
  face_buffer_size = 512;
  
  nodes[0] = aligned_real_alloc (node_buffer_size);
  nodes[1] = aligned_real_alloc (node_buffer_size); 
  nodes[2] = aligned_real_alloc (node_buffer_size); 
  nodes[3] = aligned_real_alloc (node_buffer_size); 
  nodes[4] = aligned_real_alloc (node_buffer_size); 
  nodes[5] = aligned_real_alloc (node_buffer_size); 
  nodpart = aligned_int_alloc (node_buffer_size); 
  eletype = aligned_int_alloc (element_buffer_size); 
  elenod = aligned_int_alloc (element_node_buffer_size); 
  eleidx = aligned_int_alloc (element_buffer_size+1); 
  elepart = aligned_int_alloc (element_buffer_size); 
  elemat = aligned_int_alloc (element_buffer_size); 
  facnod[0] = aligned_int_alloc (face_buffer_size); 
  facnod[1] = aligned_int_alloc (face_buffer_size); 
  facnod[2] = aligned_int_alloc (face_buffer_size); 
  facpart = aligned_int_alloc (face_buffer_size); 
  factri = aligned_int_alloc (face_buffer_size); 

  nodnum = 0;
  elenum = 0;
  facnum = 0;
  faccon = 0;
  eleidx[elenum] = 0;
}

/* grow element buffer */
void element_buffer_grow (int node_count, int element_node_count, int element_count, int triangle_count)
{
  if (node_buffer_size < nodnum + node_count)
  {
    node_buffer_size = 2 * (node_buffer_size + node_count);
    real_buffer_grow (nodes[0], nodnum, node_buffer_size);
    real_buffer_grow (nodes[1], nodnum, node_buffer_size);
    real_buffer_grow (nodes[2], nodnum, node_buffer_size);
    real_buffer_grow (nodes[3], nodnum, node_buffer_size);
    real_buffer_grow (nodes[4], nodnum, node_buffer_size);
    real_buffer_grow (nodes[5], nodnum, node_buffer_size);
    integer_buffer_grow (nodpart, nodnum, node_buffer_size);
  }

  if (element_node_buffer_size < eleidx[elenum] + element_node_count)
  {
    element_node_buffer_size = 2 * (eleidx[elenum] + element_node_count);
    integer_buffer_grow (elenod, eleidx[elenum], element_node_buffer_size);
  }

  if (element_buffer_size < elenum + element_count)
  {
    element_buffer_size = 2 * (elenum + element_count);
    integer_buffer_grow (eletype, elenum, element_buffer_size);
    integer_buffer_grow (eleidx, elenum+1, element_buffer_size+1);
    integer_buffer_grow (elepart, elenum, element_buffer_size);
    integer_buffer_grow (elemat, elenum, element_buffer_size);
  }

  if (face_buffer_size < facnum + triangle_count)
  {
    face_buffer_size = 2 * (facnum + triangle_count);
    integer_buffer_grow (facnod[0], facnum, face_buffer_size);
    integer_buffer_grow (facnod[1], facnum, face_buffer_size);
    integer_buffer_grow (facnod[2], facnum, face_buffer_size);
    integer_buffer_grow (facpart, facnum, face_buffer_size);
    integer_buffer_grow (factri, facnum, face_buffer_size);
  }
}

/* init obstacle buffer */
void obstacle_buffer_init ()
{
  obstacle_buffer_size = 256;

  trirng = aligned_int_alloc (2*obstacle_buffer_size);
  obspnt = aligned_real_alloc (3*obstacle_buffer_size);
  obsang = aligned_real_alloc (3*obstacle_buffer_size);
  obslin = aligned_real_alloc (3*obstacle_buffer_size);
  anghis = new pointer_t [obstacle_buffer_size];
  linhis = new pointer_t [obstacle_buffer_size];

  obsnum = 0;
}

/* grow obstacle buffer */
int obstacle_buffer_grow ()
{
  obstacle_buffer_size *= 2;

  integer_buffer_grow (trirng, 2*obsnum, 2*obstacle_buffer_size);
  real_buffer_grow (obspnt, 3*obsnum, 3*obstacle_buffer_size);
  real_buffer_grow (obsang, 3*obsnum, 3*obstacle_buffer_size);
  real_buffer_grow (obslin, 3*obsnum, 3*obstacle_buffer_size);
  pointer_buffer_grow (anghis, obsnum, obstacle_buffer_size);
  pointer_buffer_grow (linhis, obsnum, obstacle_buffer_size);

  return obstacle_buffer_size;
}

/* init spring buffer */
void spring_buffer_init ()
{
  spring_buffer_size = 256;
  spring_lookup_size = 1024;
  dashpot_lookup_size = 1024;
  unload_lookup_size = 1024;

  sprid = aligned_int_alloc (spring_buffer_size);
  sprmap = aligned_int_alloc (spring_buffer_size);
  sprtype = aligned_int_alloc (spring_buffer_size);
  unspring = aligned_int_alloc (spring_buffer_size);
  sprpart[0] = aligned_int_alloc (spring_buffer_size);
  sprpart[1] = aligned_int_alloc (spring_buffer_size);
  sprpnt[0][0] = aligned_real_alloc (spring_buffer_size);
  sprpnt[0][1] = aligned_real_alloc (spring_buffer_size);
  sprpnt[0][2] = aligned_real_alloc (spring_buffer_size);
  sprpnt[0][3] = aligned_real_alloc (spring_buffer_size);
  sprpnt[0][4] = aligned_real_alloc (spring_buffer_size);
  sprpnt[0][5] = aligned_real_alloc (spring_buffer_size);
  sprpnt[1][0] = aligned_real_alloc (spring_buffer_size);
  sprpnt[1][1] = aligned_real_alloc (spring_buffer_size);
  sprpnt[1][2] = aligned_real_alloc (spring_buffer_size);
  sprpnt[1][3] = aligned_real_alloc (spring_buffer_size);
  sprpnt[1][4] = aligned_real_alloc (spring_buffer_size);
  sprpnt[1][5] = aligned_real_alloc (spring_buffer_size);
  spring[0] = aligned_real_alloc (spring_lookup_size);
  spring[1] = aligned_real_alloc (spring_lookup_size);
  spridx = aligned_int_alloc (spring_buffer_size+1);
  dashpot[0] = aligned_real_alloc (dashpot_lookup_size);
  dashpot[1] = aligned_real_alloc (dashpot_lookup_size);
  dashidx = aligned_int_alloc (spring_buffer_size+1);
  unload[0] = aligned_real_alloc (unload_lookup_size);
  unload[1] = aligned_real_alloc (unload_lookup_size);
  unidx = aligned_int_alloc (spring_buffer_size+1);
  yield[0] = aligned_real_alloc (spring_buffer_size);
  yield[1] = aligned_real_alloc (spring_buffer_size);
  sprdir[0] = aligned_real_alloc (spring_buffer_size);
  sprdir[1] = aligned_real_alloc (spring_buffer_size);
  sprdir[2] = aligned_real_alloc (spring_buffer_size);
  sprdir[3] = aligned_real_alloc (spring_buffer_size);
  sprdir[4] = aligned_real_alloc (spring_buffer_size);
  sprdir[5] = aligned_real_alloc (spring_buffer_size);
  sprflg = aligned_int_alloc (spring_buffer_size);
  sproffset = aligned_int_alloc (spring_buffer_size);
  sprfric = aligned_real_alloc (spring_buffer_size);
  sprkskn = aligned_real_alloc (spring_buffer_size);
  sprsdsp[0] = aligned_real_alloc (spring_buffer_size);
  sprsdsp[1] = aligned_real_alloc (spring_buffer_size);
  sprsdsp[2] = aligned_real_alloc (spring_buffer_size);
  stroke0 = aligned_real_alloc (spring_buffer_size);
  stroke[0] = aligned_real_alloc (spring_buffer_size);
  stroke[1] = aligned_real_alloc (spring_buffer_size);
  stroke[2] = aligned_real_alloc (spring_buffer_size);
  sprfrc[0] = aligned_real_alloc (spring_buffer_size);
  sprfrc[1] = aligned_real_alloc (spring_buffer_size);
  sprfrc[2] = aligned_real_alloc (spring_buffer_size);

  sprnum = 0;
  spridx[sprnum] = 0;
  dashidx[sprnum] = 0;
  unidx[sprnum] = 0;
  springs_changed = 0;
}

/* init torsion spring buffer */
void trqspr_buffer_init ()
{
  trqspr_buffer_size = 256;
  krpy_lookup_size[0] = 1024;
  krpy_lookup_size[1] = 1024;
  krpy_lookup_size[2] = 1024;
  drpy_lookup_size[0] = 1024;
  drpy_lookup_size[1] = 1024;
  drpy_lookup_size[2] = 1024;

  trqsprid = aligned_int_alloc (trqspr_buffer_size);
  trqsprmap = aligned_int_alloc (trqspr_buffer_size);
  trqsprpart[0] = aligned_int_alloc (trqspr_buffer_size);
  trqsprpart[1] = aligned_int_alloc (trqspr_buffer_size);
  trqzdir0[0] = aligned_real_alloc (trqspr_buffer_size);
  trqzdir0[1] = aligned_real_alloc (trqspr_buffer_size);
  trqzdir0[2] = aligned_real_alloc (trqspr_buffer_size);
  trqxdir0[0] = aligned_real_alloc (trqspr_buffer_size);
  trqxdir0[1] = aligned_real_alloc (trqspr_buffer_size);
  trqxdir0[2] = aligned_real_alloc (trqspr_buffer_size);
  krpy[0][0] = aligned_real_alloc (krpy_lookup_size[0]);
  krpy[0][1] = aligned_real_alloc (krpy_lookup_size[0]);
  krpy[1][0] = aligned_real_alloc (krpy_lookup_size[1]);
  krpy[1][1] = aligned_real_alloc (krpy_lookup_size[1]);
  krpy[2][0] = aligned_real_alloc (krpy_lookup_size[2]);
  krpy[2][1] = aligned_real_alloc (krpy_lookup_size[2]);
  krpyidx[0] = aligned_int_alloc (trqspr_buffer_size+1);
  krpyidx[1] = aligned_int_alloc (trqspr_buffer_size+1);
  krpyidx[2] = aligned_int_alloc (trqspr_buffer_size+1);
  drpy[0][0] = aligned_real_alloc (drpy_lookup_size[0]);
  drpy[0][1] = aligned_real_alloc (drpy_lookup_size[0]);
  drpy[1][0] = aligned_real_alloc (drpy_lookup_size[1]);
  drpy[1][1] = aligned_real_alloc (drpy_lookup_size[1]);
  drpy[2][0] = aligned_real_alloc (drpy_lookup_size[2]);
  drpy[2][1] = aligned_real_alloc (drpy_lookup_size[2]);
  drpyidx[0] = aligned_int_alloc (trqspr_buffer_size+1);
  drpyidx[1] = aligned_int_alloc (trqspr_buffer_size+1);
  drpyidx[2] = aligned_int_alloc (trqspr_buffer_size+1);
  trqcone = aligned_int_alloc (trqspr_buffer_size);
  trqzdir1[0] = aligned_real_alloc (trqspr_buffer_size);
  trqzdir1[1] = aligned_real_alloc (trqspr_buffer_size);
  trqzdir1[2] = aligned_real_alloc (trqspr_buffer_size);
  trqxdir1[0] = aligned_real_alloc (trqspr_buffer_size);
  trqxdir1[1] = aligned_real_alloc (trqspr_buffer_size);
  trqxdir1[2] = aligned_real_alloc (trqspr_buffer_size);
  trqrpy[0] = aligned_real_alloc (trqspr_buffer_size);
  trqrpy[1] = aligned_real_alloc (trqspr_buffer_size);
  trqrpy[2] = aligned_real_alloc (trqspr_buffer_size);
  trqrpytot[0] = aligned_real_alloc (trqspr_buffer_size);
  trqrpytot[1] = aligned_real_alloc (trqspr_buffer_size);
  trqrpytot[2] = aligned_real_alloc (trqspr_buffer_size);
  trqrpyspr[0] = aligned_real_alloc (trqspr_buffer_size);
  trqrpyspr[1] = aligned_real_alloc (trqspr_buffer_size);
  trqrpyspr[2] = aligned_real_alloc (trqspr_buffer_size);
  trqrefpnt[0] = aligned_real_alloc (trqspr_buffer_size);
  trqrefpnt[1] = aligned_real_alloc (trqspr_buffer_size);
  trqrefpnt[2] = aligned_real_alloc (trqspr_buffer_size);

  trqsprnum = 0;
  krpyidx[0][trqsprnum] = 0;
  krpyidx[1][trqsprnum] = 0;
  krpyidx[2][trqsprnum] = 0;
  drpyidx[0][trqsprnum] = 0;
  drpyidx[1][trqsprnum] = 0;
  drpyidx[2][trqsprnum] = 0;
  trqspr_changed = 0;
}

/* grow torsion spring buffer */
void trqspr_buffer_grow (int krpy_lookup[3], int drpy_lookup[3])
{
  if (trqsprnum+1 >= trqspr_buffer_size)
  {
    trqspr_buffer_size *= 2;

    integer_buffer_grow(trqsprid, trqsprnum, trqspr_buffer_size);
    integer_buffer_grow(trqsprmap, trqsprnum, trqspr_buffer_size);
    integer_buffer_grow(trqsprpart[0], trqsprnum, trqspr_buffer_size);
    integer_buffer_grow(trqsprpart[1], trqsprnum, trqspr_buffer_size);
    real_buffer_grow(trqzdir0[0], trqsprnum, trqspr_buffer_size);
    real_buffer_grow(trqzdir0[1], trqsprnum, trqspr_buffer_size);
    real_buffer_grow(trqzdir0[2], trqsprnum, trqspr_buffer_size);
    real_buffer_grow(trqxdir0[0], trqsprnum, trqspr_buffer_size);
    real_buffer_grow(trqxdir0[1], trqsprnum, trqspr_buffer_size);
    real_buffer_grow(trqxdir0[2], trqsprnum, trqspr_buffer_size);
    integer_buffer_grow(krpyidx[0], trqsprnum+1, trqspr_buffer_size+1);
    integer_buffer_grow(krpyidx[1], trqsprnum+1, trqspr_buffer_size+1);
    integer_buffer_grow(krpyidx[2], trqsprnum+1, trqspr_buffer_size+1);
    integer_buffer_grow(drpyidx[0], trqsprnum+1, trqspr_buffer_size+1);
    integer_buffer_grow(drpyidx[1], trqsprnum+1, trqspr_buffer_size+1);
    integer_buffer_grow(drpyidx[2], trqsprnum+1, trqspr_buffer_size+1);
    integer_buffer_grow(trqcone, trqsprnum, trqspr_buffer_size);
    real_buffer_grow(trqzdir1[0], trqsprnum, trqspr_buffer_size);
    real_buffer_grow(trqzdir1[1], trqsprnum, trqspr_buffer_size);
    real_buffer_grow(trqzdir1[2], trqsprnum, trqspr_buffer_size);
    real_buffer_grow(trqxdir1[0], trqsprnum, trqspr_buffer_size);
    real_buffer_grow(trqxdir1[1], trqsprnum, trqspr_buffer_size);
    real_buffer_grow(trqxdir1[2], trqsprnum, trqspr_buffer_size);
    real_buffer_grow(trqrpy[0], trqsprnum, trqspr_buffer_size);
    real_buffer_grow(trqrpy[1], trqsprnum, trqspr_buffer_size);
    real_buffer_grow(trqrpy[2], trqsprnum, trqspr_buffer_size);
    real_buffer_grow(trqrpytot[0], trqsprnum, trqspr_buffer_size);
    real_buffer_grow(trqrpytot[1], trqsprnum, trqspr_buffer_size);
    real_buffer_grow(trqrpytot[2], trqsprnum, trqspr_buffer_size);
    real_buffer_grow(trqrpyspr[0], trqsprnum, trqspr_buffer_size);
    real_buffer_grow(trqrpyspr[1], trqsprnum, trqspr_buffer_size);
    real_buffer_grow(trqrpyspr[2], trqsprnum, trqspr_buffer_size);
    real_buffer_grow(trqrefpnt[0], trqsprnum, trqspr_buffer_size);
    real_buffer_grow(trqrefpnt[1], trqsprnum, trqspr_buffer_size);
    real_buffer_grow(trqrefpnt[2], trqsprnum, trqspr_buffer_size);
  }

  for (int i = 0; i < 3; i ++)
  {
    if (krpy_lookup_size[i] < krpyidx[i][trqsprnum] + krpy_lookup[i])
    {
      krpy_lookup_size[i] = 2 * (krpyidx[i][trqsprnum] + krpy_lookup[i]);
      real_buffer_grow (krpy[i][0], krpyidx[i][trqsprnum], krpy_lookup_size[i]);
      real_buffer_grow (krpy[i][1], krpyidx[i][trqsprnum], krpy_lookup_size[i]);
    }

    if (drpy_lookup_size[i] < drpyidx[i][trqsprnum] + drpy_lookup[i])
    {
      drpy_lookup_size[i] = 2 * (drpyidx[i][trqsprnum] + drpy_lookup[i]);
      real_buffer_grow (drpy[i][0], drpyidx[i][trqsprnum], drpy_lookup_size[i]);
      real_buffer_grow (drpy[i][1], drpyidx[i][trqsprnum], drpy_lookup_size[i]);
    }
  }
}

/* grow spring buffer */
void spring_buffer_grow (int spring_lookup, int dashpot_lookup, int unload_lookup)
{
  if (sprnum+1 >= spring_buffer_size)
  {
    spring_buffer_size *= 2;

    integer_buffer_grow(sprid, sprnum, spring_buffer_size);
    integer_buffer_grow(sprmap, sprnum, spring_buffer_size);
    integer_buffer_grow(sprtype, sprnum, spring_buffer_size);
    integer_buffer_grow(unspring, sprnum, spring_buffer_size);
    integer_buffer_grow(sprpart[0], sprnum, spring_buffer_size);
    integer_buffer_grow(sprpart[1], sprnum, spring_buffer_size);
    real_buffer_grow(sprpnt[0][0], sprnum, spring_buffer_size);
    real_buffer_grow(sprpnt[0][1], sprnum, spring_buffer_size);
    real_buffer_grow(sprpnt[0][2], sprnum, spring_buffer_size);
    real_buffer_grow(sprpnt[0][3], sprnum, spring_buffer_size);
    real_buffer_grow(sprpnt[0][4], sprnum, spring_buffer_size);
    real_buffer_grow(sprpnt[0][5], sprnum, spring_buffer_size);
    real_buffer_grow(sprpnt[1][0], sprnum, spring_buffer_size);
    real_buffer_grow(sprpnt[1][1], sprnum, spring_buffer_size);
    real_buffer_grow(sprpnt[1][2], sprnum, spring_buffer_size);
    real_buffer_grow(sprpnt[1][3], sprnum, spring_buffer_size);
    real_buffer_grow(sprpnt[1][4], sprnum, spring_buffer_size);
    real_buffer_grow(sprpnt[1][5], sprnum, spring_buffer_size);
    integer_buffer_grow(spridx, sprnum+1, spring_buffer_size+1);
    integer_buffer_grow(dashidx, sprnum+1, spring_buffer_size+1);
    integer_buffer_grow(unidx, sprnum+1, spring_buffer_size+1);
    real_buffer_grow(yield[0], sprnum, spring_buffer_size);
    real_buffer_grow(yield[1], sprnum, spring_buffer_size);
    real_buffer_grow(sprdir[0], sprnum, spring_buffer_size);
    real_buffer_grow(sprdir[1], sprnum, spring_buffer_size);
    real_buffer_grow(sprdir[2], sprnum, spring_buffer_size);
    real_buffer_grow(sprdir[3], sprnum, spring_buffer_size);
    real_buffer_grow(sprdir[4], sprnum, spring_buffer_size);
    real_buffer_grow(sprdir[5], sprnum, spring_buffer_size);
    integer_buffer_grow(sprflg, sprnum, spring_buffer_size);
    integer_buffer_grow(sproffset, sprnum, spring_buffer_size);
    real_buffer_grow(sprfric, sprnum, spring_buffer_size);
    real_buffer_grow(sprkskn, sprnum, spring_buffer_size);
    real_buffer_grow(sprsdsp[0], sprnum, spring_buffer_size);
    real_buffer_grow(sprsdsp[1], sprnum, spring_buffer_size);
    real_buffer_grow(sprsdsp[2], sprnum, spring_buffer_size);
    real_buffer_grow(stroke0, sprnum, spring_buffer_size);
    real_buffer_grow(stroke[0], sprnum, spring_buffer_size);
    real_buffer_grow(stroke[1], sprnum, spring_buffer_size);
    real_buffer_grow(stroke[2], sprnum, spring_buffer_size);
    real_buffer_grow(sprfrc[0], sprnum, spring_buffer_size);
    real_buffer_grow(sprfrc[1], sprnum, spring_buffer_size);
    real_buffer_grow(sprfrc[2], sprnum, spring_buffer_size);
  }

  if (spring_lookup_size < spridx[sprnum] + spring_lookup)
  {
    spring_lookup_size = 2 * (spridx[sprnum] + spring_lookup);
    real_buffer_grow (spring[0], spridx[sprnum], spring_lookup_size);
    real_buffer_grow (spring[1], spridx[sprnum], spring_lookup_size);
  }

  if (dashpot_lookup_size < dashidx[sprnum] + dashpot_lookup)
  {
    dashpot_lookup_size = 2 * (dashidx[sprnum] + dashpot_lookup);
    real_buffer_grow (dashpot[0], dashidx[sprnum], dashpot_lookup_size);
    real_buffer_grow (dashpot[1], dashidx[sprnum], dashpot_lookup_size);
  }

  if (unload_lookup_size < unidx[sprnum] + unload_lookup)
  {
    unload_lookup_size = 2 * (unidx[sprnum] + unload_lookup);
    real_buffer_grow (unload[0], unidx[sprnum], unload_lookup_size);
    real_buffer_grow (unload[1], unidx[sprnum], unload_lookup_size);
  }
}

/* init unspring buffer */
void unspring_buffer_init ()
{
  unspring_buffer_size = 256;
  tsprings_buffer_size = 1024;
  msprings_buffer_size = 1024;
  activate_buffer_size = 1024;

  tsprings = aligned_int_alloc (tsprings_buffer_size);
  tspridx = aligned_int_alloc (unspring_buffer_size+1);
  msprings = aligned_int_alloc (msprings_buffer_size);
  mspridx = aligned_int_alloc (unspring_buffer_size+1);
  unlim[0] = aligned_real_alloc (unspring_buffer_size);
  unlim[1] = aligned_real_alloc (unspring_buffer_size);
  unent = aligned_int_alloc (tsprings_buffer_size);
  unop = aligned_int_alloc (tsprings_buffer_size);
  unabs = aligned_int_alloc (tsprings_buffer_size);
  nsteps = aligned_int_alloc (tsprings_buffer_size);
  nfreq = aligned_int_alloc (tsprings_buffer_size);
  unaction = aligned_int_alloc (tsprings_buffer_size);
  activate = aligned_int_alloc (activate_buffer_size);
  actidx = aligned_int_alloc (unspring_buffer_size+1);

  unsprnum = 0;
  tspridx[unsprnum] = 0;
  mspridx[unsprnum] = 0;
  actidx[unsprnum] = 0;
}

/* grow unspring buffer */
void unspring_buffer_grow (int tsprings_increment, int msprings_increment, int activate_increment)
{
  if (unsprnum+1 >= unspring_buffer_size)
  {
    unspring_buffer_size *= 2;

    integer_buffer_grow(tspridx, unsprnum+1, unspring_buffer_size+1);
    integer_buffer_grow(mspridx, unsprnum+1, unspring_buffer_size+1);
    real_buffer_grow(unlim[0], unsprnum, unspring_buffer_size);
    real_buffer_grow(unlim[1], unsprnum, unspring_buffer_size);
    integer_buffer_grow(unent, unsprnum, unspring_buffer_size);
    integer_buffer_grow(unop, unsprnum, unspring_buffer_size);
    integer_buffer_grow(unabs, unsprnum, unspring_buffer_size);
    integer_buffer_grow(nsteps, unsprnum, unspring_buffer_size);
    integer_buffer_grow(nfreq, unsprnum, unspring_buffer_size);
    integer_buffer_grow(unaction, unsprnum, unspring_buffer_size);
    integer_buffer_grow(actidx, unsprnum+1, unspring_buffer_size+1);
  }

  if (tsprings_buffer_size < tspridx[unsprnum] + tsprings_increment)
  {
    tsprings_buffer_size = 2 * (tspridx[unsprnum] + tsprings_increment);
    integer_buffer_grow(tsprings, tspridx[unsprnum], tsprings_buffer_size);
  }

  if (msprings_buffer_size < mspridx[unsprnum] + msprings_increment)
  {
    msprings_buffer_size = 2 * (mspridx[unsprnum] + msprings_increment);
    integer_buffer_grow(msprings, mspridx[unsprnum], msprings_buffer_size);
  }

  if (activate_buffer_size < actidx[unsprnum] + activate_increment)
  {
    activate_buffer_size = 2 * (actidx[unsprnum] + activate_increment);
    integer_buffer_grow(activate, actidx[unsprnum], activate_buffer_size);
  }
}

/* init restrained particles buffer */
void restrain_buffer_init ()
{
  restrain_buffer_size = 256;

  rstpart = aligned_int_alloc (restrain_buffer_size);
  rstlin[0] = aligned_real_alloc (restrain_buffer_size);
  rstlin[1] = aligned_real_alloc (restrain_buffer_size);
  rstlin[2] = aligned_real_alloc (restrain_buffer_size);
  rstlin[3] = aligned_real_alloc (restrain_buffer_size);
  rstlin[4] = aligned_real_alloc (restrain_buffer_size);
  rstlin[5] = aligned_real_alloc (restrain_buffer_size);
  rstlin[6] = aligned_real_alloc (restrain_buffer_size);
  rstlin[7] = aligned_real_alloc (restrain_buffer_size);
  rstlin[8] = aligned_real_alloc (restrain_buffer_size);
  rstang[0] = aligned_real_alloc (restrain_buffer_size);
  rstang[1] = aligned_real_alloc (restrain_buffer_size);
  rstang[2] = aligned_real_alloc (restrain_buffer_size);
  rstang[3] = aligned_real_alloc (restrain_buffer_size);
  rstang[4] = aligned_real_alloc (restrain_buffer_size);
  rstang[5] = aligned_real_alloc (restrain_buffer_size);
  rstang[6] = aligned_real_alloc (restrain_buffer_size);
  rstang[7] = aligned_real_alloc (restrain_buffer_size);
  rstang[8] = aligned_real_alloc (restrain_buffer_size);

  rstnum = 0;
}

/* grow restrained particles buffer */
int restrain_buffer_grow ()
{
  restrain_buffer_size *= 2;

  integer_buffer_grow (rstpart, rstnum, restrain_buffer_size);
  real_buffer_grow (rstlin[0], rstnum, restrain_buffer_size);
  real_buffer_grow (rstlin[1], rstnum, restrain_buffer_size);
  real_buffer_grow (rstlin[2], rstnum, restrain_buffer_size);
  real_buffer_grow (rstlin[3], rstnum, restrain_buffer_size);
  real_buffer_grow (rstlin[4], rstnum, restrain_buffer_size);
  real_buffer_grow (rstlin[5], rstnum, restrain_buffer_size);
  real_buffer_grow (rstlin[6], rstnum, restrain_buffer_size);
  real_buffer_grow (rstlin[7], rstnum, restrain_buffer_size);
  real_buffer_grow (rstlin[8], rstnum, restrain_buffer_size);
  real_buffer_grow (rstang[0], rstnum, restrain_buffer_size);
  real_buffer_grow (rstang[1], rstnum, restrain_buffer_size);
  real_buffer_grow (rstang[2], rstnum, restrain_buffer_size);
  real_buffer_grow (rstang[3], rstnum, restrain_buffer_size);
  real_buffer_grow (rstang[4], rstnum, restrain_buffer_size);
  real_buffer_grow (rstang[5], rstnum, restrain_buffer_size);
  real_buffer_grow (rstang[6], rstnum, restrain_buffer_size);
  real_buffer_grow (rstang[7], rstnum, restrain_buffer_size);
  real_buffer_grow (rstang[8], rstnum, restrain_buffer_size);

  return restrain_buffer_size;
}

/* init joints buffer */
void joints_buffer_init ()
{
  joints_buffer_size = 256;

  jpart[0] = aligned_int_alloc (joints_buffer_size);
  jpart[1] = aligned_int_alloc (joints_buffer_size);
  jpoint[0] = aligned_real_alloc (joints_buffer_size);
  jpoint[1] = aligned_real_alloc (joints_buffer_size);
  jpoint[2] = aligned_real_alloc (joints_buffer_size);
  jreac[0] = aligned_real_alloc (joints_buffer_size);
  jreac[1] = aligned_real_alloc (joints_buffer_size);
  jreac[2] = aligned_real_alloc (joints_buffer_size);

  jnum = 0;
  joints_changed = 0;
}

/* grow joints buffer */
int joints_buffer_grow ()
{
  joints_buffer_size *= 2;

  integer_buffer_grow (jpart[0], jnum, joints_buffer_size);
  integer_buffer_grow (jpart[1], jnum, joints_buffer_size);
  real_buffer_grow (jpoint[0], jnum, joints_buffer_size);
  real_buffer_grow (jpoint[1], jnum, joints_buffer_size);
  real_buffer_grow (jpoint[2], jnum, joints_buffer_size);
  real_buffer_grow (jreac[0], jnum, joints_buffer_size);
  real_buffer_grow (jreac[1], jnum, joints_buffer_size);
  real_buffer_grow (jreac[2], jnum, joints_buffer_size);

  return joints_buffer_size;
}

/* init time series buffer */
void time_series_buffer_init ()
{
  time_series_buffer_size = 256;

  tms = new pointer_t [time_series_buffer_size];

  tmsnum = 0;
}

/* grow time series buffer */
int time_series_buffer_grow ()
{
  time_series_buffer_size *= 2;

  pointer_buffer_grow (tms, tmsnum, time_series_buffer_size);

  return time_series_buffer_size;
}

/* init load curve buffer */
void lcurve_buffer_init ()
{
  lcurve_buffer_size = 256;
  lcurve_data_size = 1024;

  lcurve[0] = aligned_real_alloc (lcurve_data_size);
  lcurve[1] = aligned_real_alloc (lcurve_data_size);
  lcidx = aligned_int_alloc (lcurve_buffer_size+1);

  lcnum = 0;
  lcidx[lcnum] = 0;
}

/* grow load curve buffer */
void lcurve_buffer_grow (int increment)
{
  if (lcnum+1 >= lcurve_buffer_size)
  {
    lcurve_buffer_size *= 2;

    integer_buffer_grow(lcidx, lcnum+1, lcurve_buffer_size+1);
  }

  if (lcurve_data_size < lcidx[lcnum] + increment)
  {
    lcurve_data_size = 2 * (lcidx[lcnum] + increment);
    real_buffer_grow (lcurve[0], lcidx[lcnum], lcurve_data_size);
    real_buffer_grow (lcurve[1], lcidx[lcnum], lcurve_data_size);
  }
}

/* get load curve from time series */
int lcurve_from_time_series (int ts)
{
  if (ts >= 0 && ts < tmsnum)
  {
    TMS *x = (TMS*)tms[ts];

    int j = x->lc;

    if (j < 0)
    {
      j = lcnum;

      x->lc = j;

      lcurve_buffer_grow (x->size);

      int i = ++lcnum;

      lcidx[i] = lcidx[i-1];
      for (int k = 0; k < x->size; k ++)
      {
	lcurve[0][lcidx[i]] = x->points[k][0];
	lcurve[1][lcidx[i]] = x->points[k][1];
	lcidx[i] ++;
      }
    }

    return j;
  }
  else return -1;
}

/* init prescribed particles motion buffer */
void prescribe_buffer_init ()
{
  prescribe_buffer_size = 256;

  prspart = aligned_int_alloc (prescribe_buffer_size);
  prslin = new pointer_t [prescribe_buffer_size];
  tmslin[0] = aligned_int_alloc (prescribe_buffer_size);
  tmslin[1] = aligned_int_alloc (prescribe_buffer_size);
  tmslin[2] = aligned_int_alloc (prescribe_buffer_size);
  linkind = aligned_int_alloc (prescribe_buffer_size);
  prsang = new pointer_t [prescribe_buffer_size];
  tmsang[0] = aligned_int_alloc (prescribe_buffer_size);
  tmsang[1] = aligned_int_alloc (prescribe_buffer_size);
  tmsang[2] = aligned_int_alloc (prescribe_buffer_size);
  angkind = aligned_int_alloc (prescribe_buffer_size);

  prsnum = 0;
}

/* grow prescribed particles motion buffer */
int prescribe_buffer_grow ()
{
  prescribe_buffer_size *= 2;

  integer_buffer_grow (prspart, prsnum, prescribe_buffer_size);
  pointer_buffer_grow (prslin, prsnum, prescribe_buffer_size);
  integer_buffer_grow (tmslin[0], prsnum, prescribe_buffer_size);
  integer_buffer_grow (tmslin[1], prsnum, prescribe_buffer_size);
  integer_buffer_grow (tmslin[2], prsnum, prescribe_buffer_size);
  integer_buffer_grow (linkind, prsnum, prescribe_buffer_size);
  pointer_buffer_grow (prsang, prsnum, prescribe_buffer_size);
  integer_buffer_grow (tmsang[0], prsnum, prescribe_buffer_size);
  integer_buffer_grow (tmsang[1], prsnum, prescribe_buffer_size);
  integer_buffer_grow (tmsang[2], prsnum, prescribe_buffer_size);
  integer_buffer_grow (angkind, prsnum, prescribe_buffer_size);

  return prescribe_buffer_size;
}

/* init history buffer */
void history_buffer_init ()
{
  history_buffer_size = 256;
  history_list_size = 1024;

  hislst = aligned_int_alloc (history_list_size);
  hisidx = aligned_int_alloc (history_buffer_size+1);
  hisent = aligned_int_alloc (history_buffer_size);
  hiskind = aligned_int_alloc (history_buffer_size);
  source[0] = aligned_real_alloc (history_buffer_size);
  source[1] = aligned_real_alloc (history_buffer_size);
  source[2] = aligned_real_alloc (history_buffer_size);
  source[3] = aligned_real_alloc (history_buffer_size);
  source[4] = aligned_real_alloc (history_buffer_size);
  source[5] = aligned_real_alloc (history_buffer_size);
  h5file = new pointer_t [history_buffer_size];
  history = new pointer_t [history_buffer_size];

  hisnum = 0;
  hisidx[hisnum] = 0;
}

/* grow history buffer */
void history_buffer_grow (int list_size)
{
  if (hisnum+1 >= history_buffer_size)
  {
    history_buffer_size *= 2;

    integer_buffer_grow (hisidx, hisnum, history_buffer_size+1);
    integer_buffer_grow (hisent, hisnum, history_buffer_size);
    integer_buffer_grow (hiskind, hisnum, history_buffer_size);
    real_buffer_grow (source[0], hisnum, history_buffer_size);
    real_buffer_grow (source[1], hisnum, history_buffer_size);
    real_buffer_grow (source[2], hisnum, history_buffer_size);
    real_buffer_grow (source[3], hisnum, history_buffer_size);
    real_buffer_grow (source[4], hisnum, history_buffer_size);
    real_buffer_grow (source[5], hisnum, history_buffer_size);
    pointer_buffer_grow (h5file, hisnum, history_buffer_size);
    pointer_buffer_grow (history, hisnum, history_buffer_size);
  }

  if (history_list_size < hisidx[hisnum] + list_size)
  {
    history_list_size = 2 * (hisidx[hisnum] + list_size);
    integer_buffer_grow (hislst, hisidx[hisnum], history_list_size);
  }
}

/* init output buffer */
void output_buffer_init ()
{
  output_buffer_size = 256;
  output_list_size = 1024;

  outmode = aligned_int_alloc (output_buffer_size);
  outpart = aligned_int_alloc (output_list_size);
  outidx = aligned_int_alloc (output_buffer_size+1);
  outent = aligned_int_alloc (output_buffer_size);
  outrest[0] = OUT_NUMBER|OUT_COLOR|OUT_DISPL|OUT_LENGTH|OUT_ORIENT|OUT_ORIENT1|OUT_ORIENT2|OUT_ORIENT3|
               OUT_LINVEL|OUT_ANGVEL|OUT_FORCE|OUT_TORQUE|OUT_F|OUT_FN|OUT_FT|OUT_SF|OUT_AREA|OUT_PAIR|OUT_SS|OUT_FF|
	       OUT_XDIR|OUT_YDIR|OUT_ZDIR|OUT_TRQROT|OUT_TRQTOT|OUT_TRQSPR|OUT_JREAC;
  outrest[1] = OUT_MODE_SPH|OUT_MODE_MESH|OUT_MODE_RB|OUT_MODE_CD|OUT_MODE_SL|OUT_MODE_ST|OUT_MODE_JT;
  outformat = OUT_FORMAT_XDMF;

  outnum = 0;
  outidx[outnum] = 0;
}

/* grow output buffer */
void output_buffer_grow (int list_size)
{
  if (outnum+1 >= output_buffer_size)
  {
    output_buffer_size *= 2;

    integer_buffer_grow (outmode, outnum, output_buffer_size);
    integer_buffer_grow (outidx, outnum, output_buffer_size+1);
    integer_buffer_grow (outent, outnum, output_buffer_size);
  }

  if (output_list_size < outidx[outnum] + list_size)
  {
    output_list_size = 2 * (outidx[outnum] + list_size);
    integer_buffer_grow (outpart, outidx[outnum], output_list_size);
  }
}

/* try shuffle ellipsoids */
static int shuffle_ellipsoids (int k)
{
  int start = ellcon, end = ellnum, mid = (start+end)/2, i, j;

  /* particles are inserted in ascending index order
   * therefore we can apply binary search */
  while (start < end)
  {
    i = part[mid];
    if (i < k) start = mid;
    else if (i > k) end = mid;
    else break;
    mid = (start+end)/2;
  }

  if (mid == ellnum || part[mid] != k) return 0; /* not found */

  /* find the range of ellipsoids used by this particle */
  for (start = mid; start > ellcon && part[start] == k; start--);
  for (end = mid; end < ellnum && part[end] == k; end++);

  /* swap initial ellipsoids with those from the found range */
  for (i = ellcon, j = start; i < start; i ++)
  {
    int co = ellcol[i];
    int pa = part[i];
    REAL c[6] = {center[0][i], center[1][i], center[2][i], center[3][i], center[4][i], center[5][i]};
    REAL r[3] = {radii[0][i], radii[1][i], radii[2][i]};
    REAL o[18] = {orient[0][i], orient[1][i], orient[2][i], orient[3][i], orient[4][i], orient[5][i],
                  orient[6][i], orient[7][i], orient[8][i], orient[9][i], orient[10][i], orient[11][i],
	          orient[12][i], orient[13][i], orient[14][i], orient[15][i], orient[16][i], orient[17][i]};

    ellcol[i] = ellcol[j];
    part[i] = part[j];
    center[0][i] = center[0][j];
    center[1][i] = center[1][j];
    center[2][i] = center[2][j];
    center[3][i] = center[3][j];
    center[4][i] = center[4][j];
    center[5][i] = center[5][j];
    radii[0][i] = radii[0][j];
    radii[1][i] = radii[1][j];
    radii[2][i] = radii[2][j];
    orient[0][i] = orient[0][j];
    orient[1][i] = orient[1][j];
    orient[2][i] = orient[2][j];
    orient[3][i] = orient[3][j];
    orient[4][i] = orient[4][j];
    orient[5][i] = orient[5][j];
    orient[6][i] = orient[6][j];
    orient[7][i] = orient[7][j];
    orient[8][i] = orient[8][j];
    orient[9][i] = orient[9][j];
    orient[10][i] = orient[10][j];
    orient[11][i] = orient[11][j];
    orient[12][i] = orient[12][j];
    orient[13][i] = orient[13][j];
    orient[14][i] = orient[14][j];
    orient[15][i] = orient[15][j];
    orient[16][i] = orient[16][j];
    orient[17][i] = orient[17][j];

    ellcol[j] = co;
    part[j] = pa;
    center[0][j] = c[0];
    center[1][j] = c[1];
    center[2][j] = c[2];
    center[3][j] = c[3];
    center[4][j] = c[4];
    center[5][j] = c[5];
    radii[0][j] = r[0];
    radii[1][j] = r[1];
    radii[2][j] = r[2];
    orient[0][j] = o[0];
    orient[1][j] = o[1];
    orient[2][j] = o[2];
    orient[3][j] = o[3];
    orient[4][j] = o[4];
    orient[5][j] = o[5];
    orient[6][j] = o[6];
    orient[7][j] = o[7];
    orient[8][j] = o[8];
    orient[9][j] = o[9];
    orient[10][j] = o[10];
    orient[11][j] = o[11];
    orient[12][j] = o[12];
    orient[13][j] = o[13];
    orient[14][j] = o[14];
    orient[15][j] = o[15];
    orient[16][j] = o[16];
    orient[17][j] = o[17];
  }

  ellcon += end-start; /* update start index of ellipsoids used in contact detection */

  return 1;
}

/* try shuffle faces */
static void shuffle_faces (int k)
{
  int start = faccon, end = facnum, mid = (start+end)/2, i, j, l;

  /* particles are inserted in ascending index order
   * therefore we can apply binary search */
  while (start < end)
  {
    i = facpart[mid];
    if (i < k) start = mid;
    else if (i > k) end = mid;
    else break;
    mid = (start+end)/2;
  }

  if (mid == facnum || facpart[mid] != k) return; /* not found */

  /* find the range of faces used by this particle */
  for (start = mid; start > faccon && facpart[start] == k; start--);
  for (end = mid; end < facnum && facpart[end] == k; end++);

  /* swap initial faces with those from the found range */
  for (i = faccon, j = start; i < start; i ++)
  {
    int no[3] = {facnod[0][i], facnod[1][i], facnod[2][i]};
    int pa = facpart[i];
    int tr = factri[i];

    facnod[0][i] = facnod[0][j];
    facnod[1][i] = facnod[1][j];
    facnod[2][i] = facnod[2][j];
    facpart[i] = facpart[j];
    factri[i] = factri[j];

    facnod[0][j] = no[0];
    facnod[1][j] = no[1];
    facnod[2][j] = no[2];
    facpart[j] = pa;
    factri[j] = tr;
  }

  /* shuffle triangles */
  for (i = faccon, j = i+end-start; i < j; i ++)
  {
    k = tricon++; /* shift start index of triangles used in contact */
    l = factri[i];

    if (k != l)
    {
      int co = tricol[k];
      int ob = triobs[k];
      REAL t[3][3] = {{tri[0][0][k], tri[0][1][k], tri[0][2][k]},
                      {tri[1][0][k], tri[1][1][k], tri[1][2][k]},
		      {tri[2][0][k], tri[2][1][k], tri[2][2][k]}};

      tricol[k] = tricol[l];
      triobs[k] = triobs[l];
      tri[0][0][k] = tri[0][0][l];
      tri[0][1][k] = tri[0][1][l];
      tri[0][2][k] = tri[0][2][l];
      tri[1][0][k] = tri[1][0][l];
      tri[1][1][k] = tri[1][1][l];
      tri[1][2][k] = tri[1][2][l];
      tri[2][0][k] = tri[2][0][l];
      tri[2][1][k] = tri[2][1][l];
      tri[2][2][k] = tri[2][2][l];

      tricol[l] = co;
      triobs[l] = ob;
      tri[0][0][l] = t[0][0];
      tri[0][1][l] = t[0][1];
      tri[0][2][l] = t[0][2];
      tri[1][0][l] = t[1][0];
      tri[1][1][l] = t[1][1];
      tri[1][2][l] = t[1][2];
      tri[2][0][l] = t[2][0];
      tri[2][1][l] = t[2][1];
      tri[2][2][l] = t[2][2];
    }
  }

  faccon += end-start; /* update start index of faces used in contact detection */
}

/* declare particle 'k' analytical */
void declare_analytical (int k)
{
  /* try shuffle ellipsoids */
  int ok = shuffle_ellipsoids (k);

  /* if no ellipsoids use this particle, try shuffle faces (and triangles) */
  if (!ok) shuffle_faces (k);
}

/* progress bar by Ross Hemsley;
 * http://www.rosshemsley.co.uk/2011/02/creating-a-progress-bar-in-c-or-any-other-console-app/ */
static void progressbar (unsigned int x, unsigned int n, unsigned int w = 50)
{
  if (n < 100)
  {
    x *= 100/n;
    n = 100;
  }

  if ((x != n) && (x % (n/100) != 0)) return;

  using namespace std;
       
  float ratio  =  x/(float)n;
  int c =  ratio * w;
	        
  cout << setw(3) << (int)(ratio*100) << "% [";
  for (int x=0; x<c; x++) cout << "=";
  for (int x=c; x<w; x++) cout << " ";
  cout << "]\r" << flush;
}

/* temporary surface pairing */
struct pair
{
  int color1;
  int color2;
  REAL iparam[NIPARAM];
};

/* material comparison by color pair */
struct cmp_pair
{
  bool operator() (const pair& a, const pair& b)
  {
    if (a.color1 == b.color1) return a.color2 < b.color2;
    else return a.color1 < b.color1;
  }
};

/* sort materials according to surface color pairs */
static void sort_materials ()
{
  std::vector<pair> v;

  v.reserve (pairnum);

  for (int i = 0; i < pairnum; i++)
  {
    pair x;

    x.color1 = pairs[2*i];
    x.color2 = pairs[2*i+1];
    for (int j = 0; j < NIPARAM; j ++) x.iparam[j] = iparam[j][i];

    v.push_back (x);
  }

  std::sort (v.begin(), v.end(), cmp_pair());

  int i = 0;

  for (std::vector<pair>::const_iterator p = v.begin(); p != v.end(); ++p, ++i)
  {
    pairs[2*i] = p->color1;
    pairs[2*i+1] = p->color2;
    for (int j = 0; j < NIPARAM; j ++) iparam[j][i] = p->iparam[j];
  }
}

/* temporary spring data structure */
struct spring_data
{
  int part[2];
  int number;
};

/* spring comparison by particle index */
struct cmp_spring
{
  bool operator() (const spring_data & a, const spring_data & b)
  {
    if (a.part[0] == b.part[0]) return a.part[1] < b.part[1];
    else return a.part[0] < b.part[0];
  }
};

/* sort springs according to particle indices */
static void sort_springs ()
{
  std::vector<spring_data> v;

  v.reserve (sprnum);

  for (int i = 0; i < sprnum; i++)
  {
    spring_data x;

    x.part[0] = sprpart[0][i];
    x.part[1] = sprpart[1][i];
    x.number = i;

    v.push_back (x);
  }

  std::sort (v.begin(), v.end(), cmp_spring());

  int *sprid; /* spring id --> number returned to user */
  int *sprtype; /* spring type */
  int *unspring; /* unspring action */
  int *sprpart[2]; /* spring constraint particle numbers */
  REAL *sprpnt[2][6]; /* spring constraint current and reference points */
  REAL *spring[2]; /* spring force lookup tables */
  int *spridx; /* spring force lookup start index */
  REAL *dashpot[2]; /* dashpot force lookup tables */
  int *dashidx; /* dashpot force lookup start index */
  REAL *unload[2]; /* unload force lookup tables */
  int *unidx; /* unload force lookup start index */
  REAL *yield[2]; /* spring yield limits */
  REAL *sprdir[6]; /* spring direction */
  int *sprflg; /* spring flags */
  int *sproffset; /* spring stroke offset time series number */
  REAL *sprfric; /* spring tangential friction coefficient */
  REAL *sprkskn; /* ratio of normal to tangential spring and dashpot parameters */
  REAL *sprsdsp[3]; /* tangential frictional displacement vector */
  REAL *stroke0; /* initial spring stroke */
  REAL *stroke[3]; /* current stroke */
  REAL *sprfrc[3]; /* spring force */

  sprid = aligned_int_alloc (spring_buffer_size);
  sprtype = aligned_int_alloc (spring_buffer_size);
  unspring = aligned_int_alloc (spring_buffer_size);
  sprpart[0] = aligned_int_alloc (spring_buffer_size);
  sprpart[1] = aligned_int_alloc (spring_buffer_size);
  sprpnt[0][0] = aligned_real_alloc (spring_buffer_size);
  sprpnt[0][1] = aligned_real_alloc (spring_buffer_size);
  sprpnt[0][2] = aligned_real_alloc (spring_buffer_size);
  sprpnt[0][3] = aligned_real_alloc (spring_buffer_size);
  sprpnt[0][4] = aligned_real_alloc (spring_buffer_size);
  sprpnt[0][5] = aligned_real_alloc (spring_buffer_size);
  sprpnt[1][0] = aligned_real_alloc (spring_buffer_size);
  sprpnt[1][1] = aligned_real_alloc (spring_buffer_size);
  sprpnt[1][2] = aligned_real_alloc (spring_buffer_size);
  sprpnt[1][3] = aligned_real_alloc (spring_buffer_size);
  sprpnt[1][4] = aligned_real_alloc (spring_buffer_size);
  sprpnt[1][5] = aligned_real_alloc (spring_buffer_size);
  spring[0] = aligned_real_alloc (spring_lookup_size);
  spring[1] = aligned_real_alloc (spring_lookup_size);
  spridx = aligned_int_alloc (spring_buffer_size+1);
  dashpot[0] = aligned_real_alloc (dashpot_lookup_size);
  dashpot[1] = aligned_real_alloc (dashpot_lookup_size);
  dashidx = aligned_int_alloc (spring_buffer_size+1);
  unload[0] = aligned_real_alloc (unload_lookup_size);
  unload[1] = aligned_real_alloc (unload_lookup_size);
  unidx = aligned_int_alloc (spring_buffer_size+1);
  yield[0] = aligned_real_alloc (spring_buffer_size);
  yield[1] = aligned_real_alloc (spring_buffer_size);
  sprdir[0] = aligned_real_alloc (spring_buffer_size);
  sprdir[1] = aligned_real_alloc (spring_buffer_size);
  sprdir[2] = aligned_real_alloc (spring_buffer_size);
  sprdir[3] = aligned_real_alloc (spring_buffer_size);
  sprdir[4] = aligned_real_alloc (spring_buffer_size);
  sprdir[5] = aligned_real_alloc (spring_buffer_size);
  sprflg = aligned_int_alloc (spring_buffer_size);
  sproffset = aligned_int_alloc (spring_buffer_size);
  sprfric = aligned_real_alloc (spring_buffer_size);
  sprkskn = aligned_real_alloc (spring_buffer_size);
  sprsdsp[0] = aligned_real_alloc (spring_buffer_size);
  sprsdsp[1] = aligned_real_alloc (spring_buffer_size);
  sprsdsp[2] = aligned_real_alloc (spring_buffer_size);
  stroke0 = aligned_real_alloc (spring_buffer_size);
  stroke[0] = aligned_real_alloc (spring_buffer_size);
  stroke[1] = aligned_real_alloc (spring_buffer_size);
  stroke[2] = aligned_real_alloc (spring_buffer_size);
  sprfrc[0] = aligned_real_alloc (spring_buffer_size);
  sprfrc[1] = aligned_real_alloc (spring_buffer_size);
  sprfrc[2] = aligned_real_alloc (spring_buffer_size);
 
  int i = 0;

  spridx[0] = dashidx[0] = unidx[0] = 0;

  for (std::vector<spring_data>::const_iterator x = v.begin(); x != v.end(); ++x, ++i)
  {
    sprid[i] = parmec::sprid[x->number];
    parmec::sprmap[sprid[i]] = i; /* reverse mapping --> id to current index */
    sprtype[i] = parmec::sprtype[x->number];
    unspring[i] = parmec::unspring[x->number];
    sprpart[0][i] = parmec::sprpart[0][x->number];
    sprpart[1][i] = parmec::sprpart[1][x->number];
    sprpnt[0][0][i] = parmec::sprpnt[0][0][x->number];
    sprpnt[0][1][i] = parmec::sprpnt[0][1][x->number];
    sprpnt[0][2][i] = parmec::sprpnt[0][2][x->number];
    sprpnt[0][3][i] = parmec::sprpnt[0][3][x->number];
    sprpnt[0][4][i] = parmec::sprpnt[0][4][x->number];
    sprpnt[0][5][i] = parmec::sprpnt[0][5][x->number];
    sprpnt[1][0][i] = parmec::sprpnt[1][0][x->number];
    sprpnt[1][1][i] = parmec::sprpnt[1][1][x->number];
    sprpnt[1][2][i] = parmec::sprpnt[1][2][x->number];
    sprpnt[1][3][i] = parmec::sprpnt[1][3][x->number];
    sprpnt[1][4][i] = parmec::sprpnt[1][4][x->number];
    sprpnt[1][5][i] = parmec::sprpnt[1][5][x->number];
    yield[0][i] = parmec::yield[0][x->number];
    yield[1][i] = parmec::yield[1][x->number];
    sprdir[0][i] = parmec::sprdir[0][x->number];
    sprdir[1][i] = parmec::sprdir[1][x->number];
    sprdir[2][i] = parmec::sprdir[2][x->number];
    sprdir[3][i] = parmec::sprdir[3][x->number];
    sprdir[4][i] = parmec::sprdir[4][x->number];
    sprdir[5][i] = parmec::sprdir[5][x->number];
    sprflg[i] = parmec::sprflg[x->number];
    sproffset[i] = parmec::sproffset[x->number];
    sprfric[i] = parmec::sprfric[x->number];
    sprkskn[i] = parmec::sprkskn[x->number];
    sprsdsp[0][i] = parmec::sprsdsp[0][x->number];
    sprsdsp[1][i] = parmec::sprsdsp[1][x->number];
    sprsdsp[2][i] = parmec::sprsdsp[2][x->number];
    stroke0[i] = parmec::stroke0[x->number];
    stroke[0][i] = parmec::stroke[0][x->number];
    stroke[1][i] = parmec::stroke[1][x->number];
    stroke[2][i] = parmec::stroke[2][x->number];
    sprfrc[0][i] = parmec::sprfrc[0][x->number];
    sprfrc[1][i] = parmec::sprfrc[1][x->number];
    sprfrc[2][i] = parmec::sprfrc[2][x->number];

    spridx[i+1] = spridx[i] + parmec::spridx[x->number+1]-parmec::spridx[x->number];
    for (int j = spridx[i], k = parmec::spridx[x->number]; j < spridx[i+1]; j ++, k ++)
    {
      spring[0][j] = parmec::spring[0][k];
      spring[1][j] = parmec::spring[1][k];
    }

    dashidx[i+1] = dashidx[i] + parmec::dashidx[x->number+1]-parmec::dashidx[x->number];
    for (int j = dashidx[i], k = parmec::dashidx[x->number]; j < dashidx[i+1]; j ++, k ++)
    {
      dashpot[0][j] = parmec::dashpot[0][k];
      dashpot[1][j] = parmec::dashpot[1][k];
    }

    unidx[i+1] = unidx[i] + parmec::unidx[x->number+1]-parmec::unidx[x->number];
    for (int j = unidx[i], k = parmec::unidx[x->number]; j < unidx[i+1]; j ++, k ++)
    {
      unload[0][j] = parmec::unload[0][k];
      unload[1][j] = parmec::unload[1][k];
    }
  }

  aligned_int_free (parmec::sprid);
  aligned_int_free (parmec::sprtype);
  aligned_int_free (parmec::unspring);
  aligned_int_free (parmec::sprpart[0]);
  aligned_int_free (parmec::sprpart[1]);
  aligned_real_free (parmec::sprpnt[0][0]);
  aligned_real_free (parmec::sprpnt[0][1]);
  aligned_real_free (parmec::sprpnt[0][2]);
  aligned_real_free (parmec::sprpnt[0][3]);
  aligned_real_free (parmec::sprpnt[0][4]);
  aligned_real_free (parmec::sprpnt[0][5]);
  aligned_real_free (parmec::sprpnt[1][0]);
  aligned_real_free (parmec::sprpnt[1][1]);
  aligned_real_free (parmec::sprpnt[1][2]);
  aligned_real_free (parmec::sprpnt[1][3]);
  aligned_real_free (parmec::sprpnt[1][4]);
  aligned_real_free (parmec::sprpnt[1][5]);
  aligned_real_free (parmec::spring[0]);
  aligned_real_free (parmec::spring[1]);
  aligned_int_free (parmec::spridx);
  aligned_real_free (parmec::dashpot[0]);
  aligned_real_free (parmec::dashpot[1]);
  aligned_int_free (parmec::dashidx);
  aligned_real_free (parmec::unload[0]);
  aligned_real_free (parmec::unload[1]);
  aligned_int_free (parmec::unidx);
  aligned_real_free (parmec::yield[0]);
  aligned_real_free (parmec::yield[1]);
  aligned_real_free (parmec::sprdir[0]);
  aligned_real_free (parmec::sprdir[1]);
  aligned_real_free (parmec::sprdir[2]);
  aligned_real_free (parmec::sprdir[3]);
  aligned_real_free (parmec::sprdir[4]);
  aligned_real_free (parmec::sprdir[5]);
  aligned_int_free (parmec::sprflg);
  aligned_int_free (parmec::sproffset);
  aligned_real_free (parmec::sprfric);
  aligned_real_free (parmec::sprkskn);
  aligned_real_free (parmec::sprsdsp[0]);
  aligned_real_free (parmec::sprsdsp[1]);
  aligned_real_free (parmec::sprsdsp[2]);
  aligned_real_free (parmec::stroke0);
  aligned_real_free (parmec::stroke[0]);
  aligned_real_free (parmec::stroke[1]);
  aligned_real_free (parmec::stroke[2]);
  aligned_real_free (parmec::sprfrc[0]);
  aligned_real_free (parmec::sprfrc[1]);
  aligned_real_free (parmec::sprfrc[2]);

  parmec::sprid = sprid;
  parmec::sprtype = sprtype;
  parmec::unspring = unspring;
  parmec::sprpart[0] = sprpart[0];
  parmec::sprpart[1] = sprpart[1];
  parmec::sprpnt[0][0] = sprpnt[0][0];
  parmec::sprpnt[0][1] = sprpnt[0][1];
  parmec::sprpnt[0][2] = sprpnt[0][2];
  parmec::sprpnt[0][3] = sprpnt[0][3];
  parmec::sprpnt[0][4] = sprpnt[0][4];
  parmec::sprpnt[0][5] = sprpnt[0][5];
  parmec::sprpnt[1][0] = sprpnt[1][0];
  parmec::sprpnt[1][1] = sprpnt[1][1];
  parmec::sprpnt[1][2] = sprpnt[1][2];
  parmec::sprpnt[1][3] = sprpnt[1][3];
  parmec::sprpnt[1][4] = sprpnt[1][4];
  parmec::sprpnt[1][5] = sprpnt[1][5];
  parmec::spring[0] = spring[0];
  parmec::spring[1] = spring[1];
  parmec::spridx = spridx;
  parmec::dashpot[0] = dashpot[0];
  parmec::dashpot[1] = dashpot[1];
  parmec::dashidx = dashidx;
  parmec::unload[0] = unload[0];
  parmec::unload[1] = unload[1];
  parmec::unidx = unidx;
  parmec::yield[0] = yield[0];
  parmec::yield[1] = yield[1];
  parmec::sprdir[0] = sprdir[0];
  parmec::sprdir[1] = sprdir[1];
  parmec::sprdir[2] = sprdir[2];
  parmec::sprdir[3] = sprdir[3];
  parmec::sprdir[4] = sprdir[4];
  parmec::sprdir[5] = sprdir[5];
  parmec::sprflg = sprflg;
  parmec::sproffset = sproffset;
  parmec::sprfric = sprfric;
  parmec::sprkskn = sprkskn;
  parmec::sprsdsp[0] = sprsdsp[0];
  parmec::sprsdsp[1] = sprsdsp[1];
  parmec::sprsdsp[2] = sprsdsp[2];
  parmec::stroke0 = stroke0;
  parmec::stroke[0] = stroke[0];
  parmec::stroke[1] = stroke[1];
  parmec::stroke[2] = stroke[2];
  parmec::sprfrc[0] = sprfrc[0];
  parmec::sprfrc[1] = sprfrc[1];
  parmec::sprfrc[2] = sprfrc[2];

#if 0 /* print spring statistics */
  int j_avg = 0, n_j = 1;
  int k_avg = 0, n_k = 1;

  for (i = 0; i < parmec::sprnum; i ++)
  {
    int j = 0, k = 0;

    if (i)
    {
      if (parmec::sprpart[0][i-1] == parmec::sprpart[0][i]) j_avg ++;
      else n_j ++;

      if (parmec::sprpart[1][i-1] == parmec::sprpart[1][i]) k_avg ++;
      else n_k ++;
    }

    /* printf ("spring %d particles = (%d, %d)\n", i, parmec::sprpart[0][i], parmec::sprpart[1][i]); */
  }

  j_avg /= n_j;
  k_avg /= n_k;

  printf ("Average number of particles per spring = (%d, %d)\n", j_avg, k_avg);
#endif
}

/* sort torsion springs according to particle indices */
static void sort_trqspr ()
{
  std::vector<spring_data> v;

  v.reserve (trqsprnum);

  for (int i = 0; i < trqsprnum; i++)
  {
    spring_data x;

    x.part[0] = trqsprpart[0][i];
    x.part[1] = trqsprpart[1][i];
    x.number = i;

    v.push_back (x);
  }

  std::sort (v.begin(), v.end(), cmp_spring());

  int *trqsprid; /* torsion spring id --> number returned to user */
  int *trqsprmap; /* map of torsion spring ids to spring indices */
  int *trqsprpart[2]; /* torsion spring particle indices */
  REAL *trqzdir0[3]; /* torsion spring z reference direction */
  REAL *trqxdir0[3]; /* torsion spring x reference direction */
  REAL *krpy[3][2]; /* kroll, kpitch, kyaw spring angle-torque lookup tables */
  int *krpyidx[3]; /* spring torque lookup start indexes */
  REAL *drpy[3][2]; /* droll, dpitch, dyaw dashpot ang. velocity-torque lookup tables */
  int *drpyidx[3]; /* dashpot torque lookup start indexes */
  int *trqcone; /* cone constraint flags */
  REAL *trqzdir1[3]; /* output: torsion spring z current direction */
  REAL *trqxdir1[3]; /* output: torsion spring x current direction */
  REAL *trqrpy[3]; /* output: torsion spring angles roll, pitch, yaw */
  REAL *trqrpytot[3]; /* output: total moments conjugate with spring angles */
  REAL *trqrpyspr[3]; /* output: spring moments wihout damper components */
  REAL *trqrefpnt[3]; /* input: reference point coordinates on part[0] or (inf, inf, inf) */

  trqsprid = aligned_int_alloc (trqspr_buffer_size);
  trqsprpart[0] = aligned_int_alloc (trqspr_buffer_size);
  trqsprpart[1] = aligned_int_alloc (trqspr_buffer_size);
  trqzdir0[0] = aligned_real_alloc (trqspr_buffer_size);
  trqzdir0[1] = aligned_real_alloc (trqspr_buffer_size);
  trqzdir0[2] = aligned_real_alloc (trqspr_buffer_size);
  trqxdir0[0] = aligned_real_alloc (trqspr_buffer_size);
  trqxdir0[1] = aligned_real_alloc (trqspr_buffer_size);
  trqxdir0[2] = aligned_real_alloc (trqspr_buffer_size);
  krpy[0][0] = aligned_real_alloc (krpy_lookup_size[0]);
  krpy[0][1] = aligned_real_alloc (krpy_lookup_size[0]);
  krpy[1][0] = aligned_real_alloc (krpy_lookup_size[1]);
  krpy[1][1] = aligned_real_alloc (krpy_lookup_size[1]);
  krpy[2][0] = aligned_real_alloc (krpy_lookup_size[2]);
  krpy[2][1] = aligned_real_alloc (krpy_lookup_size[2]);
  krpyidx[0] = aligned_int_alloc (trqspr_buffer_size+1);
  krpyidx[1] = aligned_int_alloc (trqspr_buffer_size+1);
  krpyidx[2] = aligned_int_alloc (trqspr_buffer_size+1);
  drpy[0][0] = aligned_real_alloc (drpy_lookup_size[0]);
  drpy[0][1] = aligned_real_alloc (drpy_lookup_size[0]);
  drpy[1][0] = aligned_real_alloc (drpy_lookup_size[1]);
  drpy[1][1] = aligned_real_alloc (drpy_lookup_size[1]);
  drpy[2][0] = aligned_real_alloc (drpy_lookup_size[2]);
  drpy[2][1] = aligned_real_alloc (drpy_lookup_size[2]);
  drpyidx[0] = aligned_int_alloc (trqspr_buffer_size+1);
  drpyidx[1] = aligned_int_alloc (trqspr_buffer_size+1);
  drpyidx[2] = aligned_int_alloc (trqspr_buffer_size+1);
  trqcone = aligned_int_alloc (trqspr_buffer_size);
  trqzdir1[0] = aligned_real_alloc (trqspr_buffer_size);
  trqzdir1[1] = aligned_real_alloc (trqspr_buffer_size);
  trqzdir1[2] = aligned_real_alloc (trqspr_buffer_size);
  trqxdir1[0] = aligned_real_alloc (trqspr_buffer_size);
  trqxdir1[1] = aligned_real_alloc (trqspr_buffer_size);
  trqxdir1[2] = aligned_real_alloc (trqspr_buffer_size);
  trqrpy[0] = aligned_real_alloc (trqspr_buffer_size);
  trqrpy[1] = aligned_real_alloc (trqspr_buffer_size);
  trqrpy[2] = aligned_real_alloc (trqspr_buffer_size);
  trqrpytot[0] = aligned_real_alloc (trqspr_buffer_size);
  trqrpytot[1] = aligned_real_alloc (trqspr_buffer_size);
  trqrpytot[2] = aligned_real_alloc (trqspr_buffer_size);
  trqrpyspr[0] = aligned_real_alloc (trqspr_buffer_size);
  trqrpyspr[1] = aligned_real_alloc (trqspr_buffer_size);
  trqrpyspr[2] = aligned_real_alloc (trqspr_buffer_size);
  trqrefpnt[0] = aligned_real_alloc (trqspr_buffer_size);
  trqrefpnt[1] = aligned_real_alloc (trqspr_buffer_size);
  trqrefpnt[2] = aligned_real_alloc (trqspr_buffer_size);

  int i = 0;

  krpyidx[0][0] = krpyidx[1][0] = krpyidx[2][0] =
  drpyidx[0][0] = drpyidx[1][0] = drpyidx[2][0] = 0;

  for (std::vector<spring_data>::const_iterator x = v.begin(); x != v.end(); ++x, ++i)
  {
    trqsprid[i] = parmec::trqsprid[x->number];
    parmec::trqsprmap[trqsprid[i]] = i; /* reverse mapping --> id to current index */
    trqsprpart[0][i] = parmec::trqsprpart[0][x->number];
    trqsprpart[1][i] = parmec::trqsprpart[1][x->number];
    trqzdir0[0][i] = parmec::trqzdir0[0][x->number];
    trqzdir0[1][i] = parmec::trqzdir0[1][x->number];
    trqzdir0[2][i] = parmec::trqzdir0[2][x->number];
    trqxdir0[0][i] = parmec::trqxdir0[0][x->number];
    trqxdir0[1][i] = parmec::trqxdir0[1][x->number];
    trqxdir0[2][i] = parmec::trqxdir0[2][x->number];
    trqzdir1[0][i] = parmec::trqzdir1[0][x->number];
    trqzdir1[1][i] = parmec::trqzdir1[1][x->number];
    trqzdir1[2][i] = parmec::trqzdir1[2][x->number];
    trqxdir1[0][i] = parmec::trqxdir1[0][x->number];
    trqxdir1[1][i] = parmec::trqxdir1[1][x->number];
    trqxdir1[2][i] = parmec::trqxdir1[2][x->number];
    trqrpy[0][i] = parmec::trqrpy[0][x->number];
    trqrpy[1][i] = parmec::trqrpy[1][x->number];
    trqrpy[2][i] = parmec::trqrpy[2][x->number];
    trqrpytot[0][i] = parmec::trqrpytot[0][x->number];
    trqrpytot[1][i] = parmec::trqrpytot[1][x->number];
    trqrpytot[2][i] = parmec::trqrpytot[2][x->number];
    trqrpyspr[0][i] = parmec::trqrpyspr[0][x->number];
    trqrpyspr[1][i] = parmec::trqrpyspr[1][x->number];
    trqrpyspr[2][i] = parmec::trqrpyspr[2][x->number];
    trqrefpnt[0][i] = parmec::trqrefpnt[0][x->number];
    trqrefpnt[1][i] = parmec::trqrefpnt[1][x->number];
    trqrefpnt[2][i] = parmec::trqrefpnt[2][x->number];
    trqcone[i] = parmec::trqcone[x->number];

    for (int l = 0; l < 3; l ++)
    {
      krpyidx[l][i+1] = krpyidx[l][i] + parmec::krpyidx[l][x->number+1]-parmec::krpyidx[l][x->number];
      for (int j = krpyidx[l][i], k = parmec::krpyidx[l][x->number]; j < krpyidx[l][i+1]; j ++, k ++)
      {
	krpy[l][0][j] = parmec::krpy[l][0][k];
	krpy[l][1][j] = parmec::krpy[l][1][k];
      }

      drpyidx[l][i+1] = drpyidx[l][i] + parmec::drpyidx[l][x->number+1]-parmec::drpyidx[l][x->number];
      for (int j = drpyidx[l][i], k = parmec::drpyidx[l][x->number]; j < drpyidx[l][i+1]; j ++, k ++)
      {
	drpy[l][0][j] = parmec::drpy[l][0][k];
	drpy[l][1][j] = parmec::drpy[l][1][k];
      }
    }
  }

  aligned_int_free (parmec::trqsprid);
  aligned_int_free (parmec::trqsprpart[0]);
  aligned_int_free (parmec::trqsprpart[1]);
  aligned_real_free (parmec::trqzdir0[0]);
  aligned_real_free (parmec::trqzdir0[1]);
  aligned_real_free (parmec::trqzdir0[2]);
  aligned_real_free (parmec::trqxdir0[0]);
  aligned_real_free (parmec::trqxdir0[1]);
  aligned_real_free (parmec::trqxdir0[2]);
  aligned_real_free (parmec::krpy[0][0]);
  aligned_real_free (parmec::krpy[0][1]);
  aligned_real_free (parmec::krpy[1][0]);
  aligned_real_free (parmec::krpy[1][1]);
  aligned_real_free (parmec::krpy[2][0]);
  aligned_real_free (parmec::krpy[2][1]);
  aligned_int_free (parmec::krpyidx[0]);
  aligned_int_free (parmec::krpyidx[1]);
  aligned_int_free (parmec::krpyidx[2]);
  aligned_real_free (parmec::drpy[0][0]);
  aligned_real_free (parmec::drpy[0][1]);
  aligned_real_free (parmec::drpy[1][0]);
  aligned_real_free (parmec::drpy[1][1]);
  aligned_real_free (parmec::drpy[2][0]);
  aligned_real_free (parmec::drpy[2][1]);
  aligned_int_free (parmec::drpyidx[0]);
  aligned_int_free (parmec::drpyidx[1]);
  aligned_int_free (parmec::drpyidx[2]);
  aligned_int_free (parmec::trqcone);
  aligned_real_free (parmec::trqzdir1[0]);
  aligned_real_free (parmec::trqzdir1[1]);
  aligned_real_free (parmec::trqzdir1[2]);
  aligned_real_free (parmec::trqxdir1[0]);
  aligned_real_free (parmec::trqxdir1[1]);
  aligned_real_free (parmec::trqxdir1[2]);
  aligned_real_free (parmec::trqrpy[0]);
  aligned_real_free (parmec::trqrpy[1]);
  aligned_real_free (parmec::trqrpy[2]);
  aligned_real_free (parmec::trqrpytot[0]);
  aligned_real_free (parmec::trqrpytot[1]);
  aligned_real_free (parmec::trqrpytot[2]);
  aligned_real_free (parmec::trqrpyspr[0]);
  aligned_real_free (parmec::trqrpyspr[1]);
  aligned_real_free (parmec::trqrpyspr[2]);
  aligned_real_free (parmec::trqrefpnt[0]);
  aligned_real_free (parmec::trqrefpnt[1]);
  aligned_real_free (parmec::trqrefpnt[2]);

  parmec::trqsprid = trqsprid;
  parmec::trqsprpart[0] = trqsprpart[0];
  parmec::trqsprpart[1] = trqsprpart[1];
  parmec::trqzdir0[0] = trqzdir0[0];
  parmec::trqzdir0[1] = trqzdir0[1];
  parmec::trqzdir0[2] = trqzdir0[2];
  parmec::trqxdir0[0] = trqxdir0[0];
  parmec::trqxdir0[1] = trqxdir0[1];
  parmec::trqxdir0[2] = trqxdir0[2];
  parmec::krpy[0][0] = krpy[0][0];
  parmec::krpy[0][1] = krpy[0][1];
  parmec::krpy[1][0] = krpy[1][0];
  parmec::krpy[1][1] = krpy[1][1];
  parmec::krpy[2][0] = krpy[2][0];
  parmec::krpy[2][1] = krpy[2][1];
  parmec::krpyidx[0] = krpyidx[0];
  parmec::krpyidx[1] = krpyidx[1];
  parmec::krpyidx[2] = krpyidx[2];
  parmec::drpy[0][0] = drpy[0][0];
  parmec::drpy[0][1] = drpy[0][1];
  parmec::drpy[1][0] = drpy[1][0];
  parmec::drpy[1][1] = drpy[1][1];
  parmec::drpy[2][0] = drpy[2][0];
  parmec::drpy[2][1] = drpy[2][1];
  parmec::drpyidx[0] = drpyidx[0];
  parmec::drpyidx[1] = drpyidx[1];
  parmec::drpyidx[2] = drpyidx[2];
  parmec::trqcone = trqcone;
  parmec::trqzdir1[0] = trqzdir1[0];
  parmec::trqzdir1[1] = trqzdir1[1];
  parmec::trqzdir1[2] = trqzdir1[2];
  parmec::trqxdir1[0] = trqxdir1[0];
  parmec::trqxdir1[1] = trqxdir1[1];
  parmec::trqxdir1[2] = trqxdir1[2];
  parmec::trqrpy[0] = trqrpy[0];
  parmec::trqrpy[1] = trqrpy[1];
  parmec::trqrpy[2] = trqrpy[2];
  parmec::trqrpytot[0] = trqrpytot[0];
  parmec::trqrpytot[1] = trqrpytot[1];
  parmec::trqrpytot[2] = trqrpytot[2];
  parmec::trqrpyspr[0] = trqrpyspr[0];
  parmec::trqrpyspr[1] = trqrpyspr[1];
  parmec::trqrpyspr[2] = trqrpyspr[2];
  parmec::trqrefpnt[0] = trqrefpnt[0];
  parmec::trqrefpnt[1] = trqrefpnt[1];
  parmec::trqrefpnt[2] = trqrefpnt[2];
}

/* add up prescribed body forces */
static void prescribe_body_forces (MAP *prescribed_body_forces, REAL *force[3], REAL *toruqe[3])
{
  for (MAP *item = MAP_First (prescribed_body_forces); item; item = MAP_Next (item))
  {
    struct prescribed_body_force *p = (struct prescribed_body_force*) item->data;

    p->inner_force[0] = force[0][p->particle]; /* export forces --> overwrite with the latest value */
    p->inner_force[1] = force[1][p->particle];
    p->inner_force[2] = force[2][p->particle];

    p->inner_torque[0] = torque[0][p->particle];
    p->inner_torque[1] = torque[1][p->particle];
    p->inner_torque[2] = torque[2][p->particle];

    force[0][p->particle] += p->outer_force[0]; /* import forces --> add up as force/torque vectors are zeroed at every time step */
    force[1][p->particle] += p->outer_force[1];
    force[2][p->particle] += p->outer_force[2];

    torque[0][p->particle] += p->outer_torque[0];
    torque[1][p->particle] += p->outer_torque[1];
    torque[2][p->particle] += p->outer_torque[2];
  }
}

/* init parmec library */
void init()
{
  material_buffer_init ();
  pair_buffer_init ();
  ellipsoid_buffer_init ();
  particle_buffer_init ();
  triangle_buffer_init ();
  element_buffer_init ();
  obstacle_buffer_init ();
  spring_buffer_init ();
  trqspr_buffer_init ();
  unspring_buffer_init ();
  joints_buffer_init ();
  restrain_buffer_init ();
  time_series_buffer_init ();
  lcurve_buffer_init ();
  prescribe_buffer_init ();
  history_buffer_init ();
  output_buffer_init ();

  reset ();

  ntasks = ispc_num_cores();
}

/* reset all data */
void reset ()
{
  output_frame = 0;

  stepnum = 0;
  curtime = 0.0;
  curstep = 0.0;
  curtime_output = 0.0;
  curtime_history = 0.0;

  matnum = 0;
  pairnum = 0;
  ellnum = 0;
  ellcon = 0;
  parnum = 0;
  trinum = 0;
  tricon = 0;
  nodnum = 0;
  elenum = 0;
  facnum = 0;
  faccon = 0;
  obsnum = 0;
  sprnum = 0;
  trqsprnum = 0;
  rstnum = 0;
  jnum = 0;
  tmsnum = 0;
  prsnum = 0;
  hisnum = 0;
  outnum = 0;

  springs_changed = 0; /* unset linear springs changed flag */
  trqspr_changed = 0; /* unset torsion springs changed flag */
  joints_changed = 0; /* unset joints changed flag */

  /* unselected particles default output flags */
  outrest[0] = OUT_NUMBER|OUT_COLOR|OUT_DISPL|OUT_LENGTH|OUT_ORIENT|OUT_ORIENT1|OUT_ORIENT2|OUT_ORIENT3|
               OUT_LINVEL|OUT_ANGVEL|OUT_FORCE|OUT_TORQUE|OUT_F|OUT_FN|OUT_FT|OUT_SF|OUT_AREA|OUT_PAIR|OUT_SS|OUT_FF|
	       OUT_XDIR|OUT_YDIR|OUT_ZDIR|OUT_TRQROT|OUT_TRQTOT|OUT_TRQSPR|OUT_JREAC;
  outrest[1] = OUT_MODE_SPH|OUT_MODE_MESH|OUT_MODE_RB|OUT_MODE_CD|OUT_MODE_SL|OUT_MODE_ST|OUT_MODE_JT;
  outformat = OUT_FORMAT_XDMF;

  /* zero global damping by default */
  damping[0] = damping[1] = damping[2] = damping[3] = damping[4] = damping[5] = 0.0;
  lindamp = angdamp = NULL;
  lindamptms[0] = lindamptms[1] = lindamptms[2] = -1;
  angdamptms[0] = angdamptms[2] = angdamptms[2] = -1;

  /* zero gravity by default */
  gravity[0] = gravity[1] = gravity[2] = 0.0;
  gravfunc[0] = gravfunc[1] = gravfunc[2] = NULL;
  gravtms[0] = gravtms[1] = gravtms[2] = -1;

  /* no prescribed body forces by default */
  prescribed_body_forces = NULL;

  pair_reset();

  output_reset();
}

/* run DEM simulation */
REAL dem (REAL duration, REAL step, REAL *interval, pointer_t *interval_func, int *interval_tms, char *prefix, int verbose, double adaptive)
{
  REAL time, dt, step0, step1;
  REAL auxiliary_interval[2];
  timing tt;

  timerstart (&tt);

  if ((interval_func || interval_tms) && !interval)
  {
    interval = auxiliary_interval;
  }

  if (prefix)
  {
    static char prefix0[1024] = {'\0'};
    ASSERT (!(strcmp(prefix0, prefix) != 0 && curtime > 0.0),
      "ERROR: output file prefix has changed for time > 0.0; "
      "INFO: the output prefix can only change at time 0.0 or after RESET()\n");
    strncpy (prefix0, prefix, 1024);

    char *out;
    int len, i;
    len = strlen (output_path);
    for (i = len-1; i >= 0; i --)
    {
      if (output_path[i] == '/' || /* unix */
	  output_path[i] == '\\') break; /* windows */
    }
    len = (i+1) + strlen (prefix) + 1;
    ERRMEM (out = new char [len]);
    output_path[i+1] = '\0';
    strcpy (out, output_path);
    strcpy (out+i+1, prefix);
    delete output_path;
    output_path = out;
  }

  REAL *icenter[6] = {center[0]+ellcon, center[1]+ellcon, center[2]+ellcon, center[3]+ellcon, center[4]+ellcon, center[5]+ellcon};
  REAL *iradii[3] = {radii[0]+ellcon, radii[1]+ellcon, radii[2]+ellcon};
  REAL *iorient[18] = {orient[0]+ellcon, orient[1]+ellcon, orient[2]+ellcon, orient[3]+ellcon, orient[4]+ellcon, orient[5]+ellcon,
		       orient[6]+ellcon, orient[7]+ellcon, orient[8]+ellcon, orient[9]+ellcon, orient[10]+ellcon, orient[11]+ellcon,
		       orient[12]+ellcon, orient[13]+ellcon, orient[14]+ellcon, orient[15]+ellcon, orient[16]+ellcon, orient[17]+ellcon};
  REAL *itri[3][3] = {{tri[0][0]+tricon, tri[0][1]+tricon, tri[0][2]+tricon},
		      {tri[1][0]+tricon, tri[1][1]+tricon, tri[1][2]+tricon},
		      {tri[2][0]+tricon, tri[2][1]+tricon, tri[2][2]+tricon}};
  int *ifacnod[3] = {facnod[0]+faccon, facnod[1]+faccon, facnod[2]+faccon};

  sort_materials ();

  if (springs_changed)
  {
    sort_springs ();

    springs_changed = 0;
  }

  if (trqspr_changed)
  {
    sort_trqspr ();

    trqspr_changed = 0;
  }

  if (joints_changed)
  {
    reset_joints_symbolic (jnum, jpart);

    joints_changed = 0;
  }

  if (curtime == 0.0)
  {
    step0 = step;

    if (interval)
    {
      output_files ();

      output_history ();
    }

    euler (ntasks, parnum, angular, linear, rotation, position, 0.5*step0);

    shapes (ntasks, ellnum, part, center, radii, orient, nodnum, nodes,
            nodpart, NULL, facnum, facnod, factri, tri, rotation, position);

    obstaclev (obsnum, obsang, obslin, anghis, linhis, 0.5*step0);

    obstacles (obsnum, trirng, obspnt, obsang, obslin, tri, step0);
  }
  else
  {
    step0 = curstep;
  }

  invert_inertia (ntasks, parnum, inertia, inverse, mass, invm);

  restrain_velocities (ntasks, rstnum, rstpart, rstlin, rstang, linear, angular, rotation);

  partitioning *tree = partitioning_create (ntasks, ellnum-ellcon, icenter);

  /* time stepping */
  for (time = 0.0; time < duration; time += 0.5*(step0+step1), curtime += 0.5*(step0+step1), step0 = step1)
  {
    if (partitioning_store (ntasks, tree, ellnum-ellcon, ellcol+ellcon, part+ellcon, icenter, iradii, iorient) > 0)
    {
      partitioning_destroy (tree);

      tree = partitioning_create (ntasks, ellnum-ellcon, icenter);

      ASSERT (partitioning_store (ntasks, tree, ellnum-ellcon, ellcol+ellcon, part+ellcon, icenter, iradii, iorient) == 0, "Repartitioning failed");
    }

    condet (ntasks, tree, master, parnum, ellnum-ellcon, ellcol+ellcon, part+ellcon,
            icenter, iradii, iorient, trinum-tricon, tricol+tricon, triobs+tricon, itri);

    read_gravity_and_damping (curtime, tms, gravfunc, gravtms, gravity, lindamp, lindamptms, angdamp, angdamptms, damping);

    forces (ntasks, master, slave, parnum, angular, linear, rotation, position, inertia, inverse, mass, invm, obspnt, obslin,
            obsang, parmat, mparam, pairnum, pairs, ikind, iparam, step0, sprnum, sprtype, unspring, sprmap, sprpart, sprpnt,
	    spring, spridx, dashpot, dashidx, unload, unidx, yield, sprdir, sprflg, sproffset, sprfric, sprkskn, sprsdsp, stroke0,
	    stroke, sprfrc, lcurve, lcidx, gravity, trqsprnum, trqsprpart, trqzdir0, trqxdir0, krpy, krpyidx, drpy, drpyidx, trqcone,
	    trqzdir1, trqxdir1, trqrpy, trqrpytot, trqrpyspr, force, torque, kact, kmax, emax, krot, (adaptive > 0.0 && adaptive <= 1.0),
	    unsprnum, tsprings, tspridx, msprings, mspridx, unlim, unent, unop, unabs, nsteps, nfreq, unaction, activate, actidx,
	    stepnum, curtime);

    prescribe_body_forces (prescribed_body_forces, force, torque);

    if (adaptive > 0.0 && adaptive <= 1.0)
    {
      step1 = adaptive_timestep (ntasks, parnum, mass, inertia, kact, kmax, emax, krot, step0, adaptive);
    }
    else
    {
      step1 = step0;
    }

    solve_joints (jnum, jpart, jpoint, jreac, parnum, position, rotation, inertia,
      inverse, mass, invm, damping, linear, angular, force, torque, step0, step1);

    restrain_forces (ntasks, rstnum, rstpart, rstlin, rstang, force, torque);

    prescribe_acceleration (prsnum, tms, prspart, prslin, tmslin, linkind, prsang,
                            tmsang, angkind, curtime, mass, inertia, force, torque);

    dynamics (ntasks, master, slave, parnum, angular, linear, rotation, position,
      inertia, inverse, mass, invm, damping, force, torque, flags, step0, step1);

    prescribe_velocity (prsnum, tms, prspart, prslin, tmslin, linkind, prsang,
                        tmsang, angkind, curtime, rotation, linear, angular);

    if (interval && interval_func && interval_func[0])
    {
      interval[0] = current_interval(interval_func[0], curtime);
    }
    else if (interval && interval_tms && interval_tms[0] >= 0 && interval_tms[0] < tmsnum)
    {
      interval[0] = TMS_Value ((TMS*)tms[interval_tms[0]], curtime);
    }

    if (interval && interval_func && interval_func[1])
    {
      interval[1] = current_interval(interval_func[1], curtime);
    }
    else if (interval && interval_tms && interval_tms[1] >= 0 && interval_tms[1] < tmsnum)
    {
      interval[1] = TMS_Value ((TMS*)tms[interval_tms[1]], curtime);
    }

    if (interval && curtime >= curtime_output + interval[0]) /* full update, due to output */
    {
      shapes (ntasks, ellnum, part, center, radii, orient,
	      nodnum, nodes, nodpart, NULL, facnum, facnod,
	      factri, tri, rotation, position);
    }
    else /* partial update, skipping analytical particles */
    {
      shapes (ntasks, ellnum-ellcon, part+ellcon, icenter, iradii, iorient,
	      nodnum, nodes, nodpart, flags, facnum-faccon, ifacnod,
	      factri+faccon, tri, rotation, position);
    }

    obstaclev (obsnum, obsang, obslin, anghis, linhis, curtime+step0);

    obstacles (obsnum, trirng, obspnt, obsang, obslin, tri, 0.5*(step0+step1));

    if (interval && curtime >= curtime_output + interval[0])
    {
      output_files ();

      curtime_output += interval[0];
    }

    if (interval && curtime >= curtime_history + interval[1])
    {
      output_history ();

      curtime_history += interval[1];
    }

    if (verbose) progressbar (2.0*time/(step0+step1), 2.0*duration/(step0+step1));
  }

  partitioning_destroy (tree);

  curstep = step1;

  stepnum ++;

  dt = timerend (&tt);

  if (verbose) printf("[ ===             %10.3f sec                    === ]\n", dt);

  return dt;
}

#ifdef __cplusplus
} /* namespace */
#endif

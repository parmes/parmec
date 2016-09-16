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
#include "macros.h"
#include "parmec.h"
#include "parmec_ispc.h"
#include "condet_ispc.h"

using namespace parmec;
using namespace ispc; /* ISPC calls are used below */

namespace parmec
{
char *outpath; /* output path */

int threads; /* number of hardware threads */

REAL curtime; /* current time */

REAL gravity[3]; /* gravity vector */

int matnum; /* number of materials */
REAL *mparam[NMAT]; /* material parameters */
int material_buffer_size; /* size of the buffer */

int pairnum; /* number of pairings */
int *pairs; /* color pairs */
int *ikind; /* interaction kind */
REAL *iparam[NIPARAM]; /* interaction parameters */
callback_t *uforce; /* user force callbacks */
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
int *flags; /* particle flags */
ispc::master_conpnt *master; /* master contact points */
ispc::slave_conpnt *slave; /* slave contact points */
int particle_buffer_size; /* size of the buffer */

int trinum; /* number of triangles */
int tricon; /* index of the first triangle used in contact detection */
int *tricol; /* triangle color */
int *triobs; /* triangle obstacle */
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
callback_t *anghis; /* angular velocity history */
callback_t *linhis; /* linear velocity history */
int obstacle_buffer_size; /* size of the buffer */

int sprnum; /* number of spring constraints */
int *sprpart[2]; /* spring constraint particle numbers */
REAL *sprpnt[2][6]; /* spring constraint current and reference points */
REAL *spring[2]; /* spring force lookup tables */
int *spridx; /* spring force lookup start index */
REAL *dashpot[2]; /* dashpot force lookup tables */
int *dashidx; /* dashpot force lookup start index */
REAL *sprdir[3]; /* spring direction */
int *sprdirup; /* spring direction update flag */
REAL *stroke0; /* initial spring stroke */
int spring_buffer_size; /* size of the spring constraint buffer */
int spring_lookup_size; /* size of the spring force lookup tables */
int dashpot_lookup_size; /* size of the dashpot force lookup tables */

int cnsnum; /* number of constraints */
int *cnspart; /* constrained particle numbers */
REAL *cnslin[9]; /* constrained linear directions */
REAL *cnsang[9]; /* constrained angular directions */
int constrain_buffer_size; /* size of constrained particles buffer */

int prsnum; /* number of particles with prescribed motion */
int *prspart; /* prescribed motion particle numbers */
callback_t *prslin; /* prescribed linear motion time history callbacks */
int *linkind; /* prescribied linear motion signal kind: 0-velocity, 1-acceleration */
callback_t *prsang; /* prescribed angular motion time history callbacks */
int *angkind; /* prescribied angular motion signal kind: 0-velocity, 1-acceleration */
int prescribe_buffer_size; /* size of prescribed particle motion buffer */

int hisnum; /* number of time histories */
int *hispart; /* history particle lists */
int *hisidx; /* history particle list start index */
int *hisent; /* history entity */
int *hiskind; /* history kind */
REAL *source[6]; /* source sphere or box definition or optional point */
callback_t *history; /* Python list storing history */
int history_buffer_size; /* size of history buffer */
int history_list_size; /* size of history particle lists buffer */

int outnum; /* number of output lists */
int *outpart; /* output particle lists */
int *outidx; /* output particle list start index */
int *outent; /* output entities */
int outrest; /* default output entities for unlisted particles */
int output_buffer_size; /* size of output buffer */
int output_list_size; /* size of output particle lists buffer */

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

/* grow callback buffer */
void callback_buffer_grow (callback_t* &src, int num, int size)
{
  callback_t *dst;

  ERRMEM (dst = new callback_t [size]);
  memcpy (dst, src, sizeof (callback_t)*num);
  delete src;
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
int material_buffer_grow ()
{
  material_buffer_size *= 2;

  for (int i = 0; i < NMAT; i ++)
  {
    real_buffer_grow (mparam[i], matnum, material_buffer_size);
  }
}

/* initialize pairing buffer */
int pair_buffer_init ()
{
  pair_buffer_size = 64;

  pairs = aligned_int_alloc (2*pair_buffer_size);
  ikind = aligned_int_alloc (pair_buffer_size);
  for (int i = 0; i < NIPARAM; i ++)
  {
    iparam[i] = aligned_real_alloc (pair_buffer_size);
  }
  uforce = new callback_t [pair_buffer_size];
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
  callback_buffer_grow (uforce, pairnum, pair_buffer_size);

  return pair_buffer_size;
}

/* reset to default pair */
void pair_reset ()
{
  pairnum = 1; /* default */
  pairs[0] = 0;
  pairs[1] = 0;
  ikind[0] = GRANULAR_FORCE;
  iparam[SPRING][0] = 1E6;
  iparam[DAMPER][0] = 1.0;
  iparam[FRISTAT][0] = 0.0;
  iparam[FRIDYN][0] = 0.0;
  iparam[FRIROL][0] = 0.0;
  iparam[FRIDRIL][0] = 0.0;
  iparam[KSKN][0] = 0.5;
}

/* initialize ellipsoid buffer */
int ellipsoid_buffer_init ()
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
  integer_buffer_grow (flags, parnum, particle_buffer_size);
  master = master_alloc (master, parnum, particle_buffer_size);
  slave = slave_alloc (slave, parnum, particle_buffer_size);

  return particle_buffer_size;
}

/* init triangle buffer */
int triangle_buffer_init ()
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
int element_buffer_init ()
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
  eleidx = aligned_int_alloc (element_buffer_size); 
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
int obstacle_buffer_init ()
{
  obstacle_buffer_size = 256;

  trirng = aligned_int_alloc (2*obstacle_buffer_size);
  obspnt = aligned_real_alloc (3*obstacle_buffer_size);
  obsang = aligned_real_alloc (3*obstacle_buffer_size);
  obslin = aligned_real_alloc (3*obstacle_buffer_size);
  anghis = new callback_t [obstacle_buffer_size];
  linhis = new callback_t [obstacle_buffer_size];

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
  callback_buffer_grow (anghis, obsnum, obstacle_buffer_size);
  callback_buffer_grow (linhis, obsnum, obstacle_buffer_size);

  return obstacle_buffer_size;
}

/* init spring buffer */
int spring_buffer_init ()
{
  spring_buffer_size = 256;
  spring_lookup_size = 1024;
  dashpot_lookup_size = 1024;

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
  sprdir[0] = aligned_real_alloc (spring_buffer_size);
  sprdir[1] = aligned_real_alloc (spring_buffer_size);
  sprdir[2] = aligned_real_alloc (spring_buffer_size);
  sprdirup = aligned_int_alloc (spring_buffer_size);
  stroke0 = aligned_real_alloc (spring_buffer_size);

  sprnum = 0;
  spridx[sprnum] = 0;
  dashidx[sprnum] = 0;
}

/* grow spring buffer */
void spring_buffer_grow (int spring_lookup, int dashpot_lookup)
{
  if (sprnum+1 >= spring_buffer_size)
  {
    spring_buffer_size *= 2;

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
    real_buffer_grow(sprdir[0], sprnum, spring_buffer_size);
    real_buffer_grow(sprdir[1], sprnum, spring_buffer_size);
    real_buffer_grow(sprdir[2], sprnum, spring_buffer_size);
    integer_buffer_grow(sprdirup, sprnum, spring_buffer_size);
    real_buffer_grow(stroke0, sprnum, spring_buffer_size);
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
}

/* init constrained particles buffer */
int constrain_buffer_init ()
{
  constrain_buffer_size = 256;

  cnspart = aligned_int_alloc (constrain_buffer_size);
  cnslin[0] = aligned_real_alloc (constrain_buffer_size);
  cnslin[1] = aligned_real_alloc (constrain_buffer_size);
  cnslin[2] = aligned_real_alloc (constrain_buffer_size);
  cnslin[3] = aligned_real_alloc (constrain_buffer_size);
  cnslin[4] = aligned_real_alloc (constrain_buffer_size);
  cnslin[5] = aligned_real_alloc (constrain_buffer_size);
  cnslin[6] = aligned_real_alloc (constrain_buffer_size);
  cnslin[7] = aligned_real_alloc (constrain_buffer_size);
  cnslin[8] = aligned_real_alloc (constrain_buffer_size);
  cnsang[0] = aligned_real_alloc (constrain_buffer_size);
  cnsang[1] = aligned_real_alloc (constrain_buffer_size);
  cnsang[2] = aligned_real_alloc (constrain_buffer_size);
  cnsang[3] = aligned_real_alloc (constrain_buffer_size);
  cnsang[4] = aligned_real_alloc (constrain_buffer_size);
  cnsang[5] = aligned_real_alloc (constrain_buffer_size);
  cnsang[6] = aligned_real_alloc (constrain_buffer_size);
  cnsang[7] = aligned_real_alloc (constrain_buffer_size);
  cnsang[8] = aligned_real_alloc (constrain_buffer_size);

  cnsnum = 0;
}

/* grow constrained particles buffer */
int constrain_buffer_grow ()
{
  constrain_buffer_size *= 2;

  integer_buffer_grow (cnspart, cnsnum, constrain_buffer_size);
  real_buffer_grow (cnslin[0], cnsnum, constrain_buffer_size);
  real_buffer_grow (cnslin[1], cnsnum, constrain_buffer_size);
  real_buffer_grow (cnslin[2], cnsnum, constrain_buffer_size);
  real_buffer_grow (cnslin[3], cnsnum, constrain_buffer_size);
  real_buffer_grow (cnslin[4], cnsnum, constrain_buffer_size);
  real_buffer_grow (cnslin[5], cnsnum, constrain_buffer_size);
  real_buffer_grow (cnslin[6], cnsnum, constrain_buffer_size);
  real_buffer_grow (cnslin[7], cnsnum, constrain_buffer_size);
  real_buffer_grow (cnslin[8], cnsnum, constrain_buffer_size);
  real_buffer_grow (cnsang[0], cnsnum, constrain_buffer_size);
  real_buffer_grow (cnsang[1], cnsnum, constrain_buffer_size);
  real_buffer_grow (cnsang[2], cnsnum, constrain_buffer_size);
  real_buffer_grow (cnsang[3], cnsnum, constrain_buffer_size);
  real_buffer_grow (cnsang[4], cnsnum, constrain_buffer_size);
  real_buffer_grow (cnsang[5], cnsnum, constrain_buffer_size);
  real_buffer_grow (cnsang[6], cnsnum, constrain_buffer_size);
  real_buffer_grow (cnsang[7], cnsnum, constrain_buffer_size);
  real_buffer_grow (cnsang[8], cnsnum, constrain_buffer_size);

  return constrain_buffer_size;
}

/* init prescribed particles motion buffer */
int prescribe_buffer_init ()
{
  prescribe_buffer_size = 256;

  prspart = aligned_int_alloc (prescribe_buffer_size);
  prslin = new callback_t [prescribe_buffer_size];
  linkind = aligned_int_alloc (prescribe_buffer_size);
  prsang = new callback_t [prescribe_buffer_size];
  angkind = aligned_int_alloc (prescribe_buffer_size);

  prsnum = 0;
}

/* grow prescribed particles motion buffer */
int prescribe_buffer_grow ()
{
  prescribe_buffer_size *= 2;

  integer_buffer_grow (prspart, prsnum, prescribe_buffer_size);
  callback_buffer_grow (prslin, prsnum, prescribe_buffer_size);
  integer_buffer_grow (linkind, prsnum, prescribe_buffer_size);
  callback_buffer_grow (prsang, prsnum, prescribe_buffer_size);
  integer_buffer_grow (angkind, prsnum, prescribe_buffer_size);

  return prescribe_buffer_size;
}

/* init history buffer */
int history_buffer_init ()
{
  history_buffer_size = 256;
  history_list_size = 1024;

  hispart = aligned_int_alloc (history_list_size);
  hisidx = aligned_int_alloc (history_buffer_size+1);
  hisent = aligned_int_alloc (history_buffer_size);
  hiskind = aligned_int_alloc (history_buffer_size);
  source[0] = aligned_real_alloc (history_buffer_size);
  source[1] = aligned_real_alloc (history_buffer_size);
  source[2] = aligned_real_alloc (history_buffer_size);
  source[3] = aligned_real_alloc (history_buffer_size);
  source[4] = aligned_real_alloc (history_buffer_size);
  source[5] = aligned_real_alloc (history_buffer_size);
  history = new callback_t [history_buffer_size];

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
    callback_buffer_grow (history, hisnum, history_buffer_size);
  }

  if (history_list_size < hisidx[hisnum] + list_size)
  {
    history_list_size = 2 * (hisidx[hisnum] + list_size);
    integer_buffer_grow (hispart, hisidx[hisnum], history_list_size);
  }
}

/* reset all data */
void reset_all_data ()
{
  curtime = 0.0;

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
  cnsnum = 0;
  prsnum = 0;
  hisnum = 0;
  outnum = 0;
  outrest = OUT_COLOR|OUT_DISP|OUT_LINVEL|OUT_ANGVEL|OUT_FORCE|OUT_TORQUE;

  pair_reset();
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
} /* end of namespace */

int main (int argc, char *argv[])
{
  if (argc == 1)
  {
    printf ("SYNOPSIS: parmec [-threads n] path/to/file.py\n");
    printf ("         -threads n: number of threads (default: hardware supported maximum)\n");
    return 1;
  }
  else
  {
    material_buffer_init ();
    pair_buffer_init ();
    ellipsoid_buffer_init ();
    particle_buffer_init ();
    triangle_buffer_init ();
    element_buffer_init ();
    obstacle_buffer_init ();
    spring_buffer_init ();
    constrain_buffer_init ();
    prescribe_buffer_init ();
    reset_all_data ();

    if (strcmp (argv[1], "-threads") == 0 && argc > 2)
    {
      threads = atoi (argv[2]);
    }
    else threads = ispc_num_cores();

    int len = strlen (argv[argc-1]);
    ERRMEM (outpath = new char [len+1]);
    strcpy (outpath, argv[argc-1]);
   
    if (outpath[len-3] != '.' ||
        outpath[len-2] != 'p' ||
	outpath[len-1] != 'y')
    {
      fprintf (stderr, "ERROR: input file does not have '.py' extension!\n");
      return 2;
    }
    else outpath[len-3] = '\0';

    input (argv[argc-1]);
  }

  return 0;
}

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
int *elenod; /* element nodes */
int *eleidx; /* element nodes index */
int *elepart; /* element particle index */
int numfac; /* number of faces (triangulated) */
int *facnod; /* face nodes */
int *factri; /* face to triangle mapping */

int obsnum; /* number of obstacles */
int *trirng; /* triangles range */
REAL *obspnt; /* obstacle spatial points */
REAL *obsang; /* obstacle angular velocities at t and t+h */
REAL *obslin; /* obstacle linear velocities at t and t+h */
callback_t *anghis; /* angular velocity history */
callback_t *linhis; /* linear velocity history */
int obstacle_buffer_size; /* size of the obstacle buffer */

/* ISPC calls are used below */
using namespace ispc;

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

/* grow triangle buffer */
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
  obsnum = 0;

  pair_reset();
}

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
    obstacle_buffer_init ();
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

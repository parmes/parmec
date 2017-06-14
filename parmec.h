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

#include "condet_ispc.h"
#include "constants.h"
#include "map.h"

#ifndef __parmec__
#define __parmec__

#ifdef __cplusplus
namespace parmec { /* namespace */
#endif

typedef void* pointer_t; /* generall Python callback */

extern char **argv; /* input arguments */
extern int argc; /* input arguments count */

extern int ntasks; /* number of tasks */

extern char *output_path; /* output path */
extern int output_frame; /* output files frame */

extern REAL curtime; /* current time */
extern REAL curstep; /* current step */
extern REAL curtime_output; /* current output time */
extern REAL curtime_history; /* current history time */

extern int matnum; /* number of materials */
extern REAL *mparam[NMAT]; /* material parameters */
extern int material_buffer_size; /* size of the buffer */
extern int material_buffer_grow (); /* grow buffer */

extern int pairnum; /* number of pairings */
extern int *pairs; /* color pairs */
extern int *ikind; /* interaction kind */
extern REAL *iparam[NIPARAM]; /* interaction parameters */
extern pointer_t *uforce; /* user force callbacks */
extern int pair_buffer_size; /* size of the buffer */
extern int pair_buffer_grow (); /* grow buffer */
extern void pair_reset (); /* reset pairing defaults */

extern int ellnum; /* number of ellipsoids */
extern int ellcon; /* index of the first ellipsoid used in contact detection */
extern int *ellcol; /* ellipsoid color */
extern int *part; /* ellipsoid to particle map */
extern REAL *center[6]; /* ellipsoid current and reference centers */
extern REAL *radii[3]; /* ellipsoid radii; radii[1][] < 0 indicates sphere */
extern REAL *orient[18]; /* ellipsoid current and reference orientations */
extern int ellipsoid_buffer_size; /* size of the buffer */
extern int ellipsoid_buffer_grow (); /* grow buffer */

extern int parnum; /* particle count */
extern int *parmat; /* particle material */
extern REAL *angular[6]; /* angular velocities (referencial, spatial) */
extern REAL *linear[3]; /* linear velocities */
extern REAL *rotation[9]; /* rotation operators */
extern REAL *position[6]; /* mass center current and reference positions */
extern REAL *inertia[9]; /* inertia tensors */
extern REAL *inverse[9]; /* inverse inertia tensors */
extern REAL *mass; /* scalar mass */
extern REAL *invm; /* inverse scalar mass */
extern REAL *force[3]; /* total spatial force */
extern REAL *torque[3]; /* total spatial torque */
extern REAL *kact[3]; /* time step control --> number of active constraints per particle */
extern REAL *kmax; /* time step control --> maximum stiffness coefficient per particle */
extern REAL *emax; /* time step control --> maximum damper coefficient per particle */
extern REAL *krot[6]; /* time step control --> symmetric rotational unit stiffness matrix per particle */
extern int *flags; /* particle flags */
extern ispc::master_conpnt *master; /* master contact points */
extern ispc::slave_conpnt *slave; /* slave contact points */
extern int particle_buffer_size; /* size of the buffer */
extern int particle_buffer_grow (); /* grow buffer */

extern int trinum; /* number of triangles */
extern int tricon; /* index of the first triangle used in contact detection */
extern int *tricol; /* triangle color */
extern int *triobs; /* triangle obstacle --> <0 - moving obstalce (-index-2), -1 - static obstacle, >= 0 - triangulated particle */
extern REAL *tri[3][3]; /* triangle vertices */
extern int triangle_buffer_size; /* size of the buffer */
extern int triangle_buffer_grow (); /* grow buffer */

extern int nodnum; /* number of nodes */
extern REAL *nodes[6]; /* current and reference nodes */
extern int *nodpart; /* node particle index */
extern int elenum; /* number of elements */
extern int *eletype; /* element type (4, 5, 6, 8) */
extern int *elenod; /* element nodes */
extern int *eleidx; /* element nodes start index */
extern int *elepart; /* element particle index */
extern int *elemat; /* element material index */
extern int facnum; /* number of faces (triangulated) */
extern int faccon; /* index of the first face used in contact detection */
extern int *facnod[3]; /* face nodes */
extern int *facpart; /* face particle index */
extern int *factri; /* face to triangle mapping */
extern int node_buffer_size; /* size of the nodes buffer */
extern int element_node_buffer_size; /* size of the element nodes buffer */
extern int element_buffer_size; /* size of the element buffers */
extern int face_buffer_size; /* size of the face buffer */
extern void element_buffer_grow (int node_count, int element_node_count, int element_count, int triangle_count); /* grow buffer */

extern int obsnum; /* number of obstacles */
extern int *trirng; /* triangles range */
extern REAL *obspnt; /* obstacle spatial points */
extern REAL *obslin; /* obstacle linear velocities at t and t+h */
extern REAL *obsang; /* obstacle angular velocities at t and t+h */
extern pointer_t *linhis; /* linear velocity history */
extern pointer_t *anghis; /* angular velocity history */
extern int obstacle_buffer_size; /* size of the buffer */
extern int obstacle_buffer_grow (); /* grow buffer */

extern int sprnum; /* number of spring constraints */
extern int *sprid; /* spring id --> number returned to user */
extern int *sprmap; /* map of spring ids to spring indices */
extern int *sprtype; /* spring type */
extern int *sprpart[2]; /* spring constraint particle numbers */
extern REAL *sprpnt[2][6]; /* spring constraint current and reference points */
extern REAL *spring[2]; /* spring force lookup tables */
extern int *spridx; /* spring force lookup start index */
extern REAL *dashpot[2]; /* dashpot force lookup tables */
extern int *dashidx; /* dashpot force lookup start index */
extern REAL *unload[2]; /* spring unloading lookup tables */
extern int *unidx; /* spring unloading lookup start index */
extern REAL *yield[2]; /* spring yield limits: 0 compression and 1 tension */
extern REAL *sprdir[3]; /* spring direction */
extern int *sprflg; /* spring flags */
extern REAL *stroke0; /* initial spring stroke */
extern REAL *stroke[3]; /* current stroke: 0 current, 1 total compression, 2 total tension */
extern REAL *sprfrc[2]; /* total and spring force magnitude */
extern int springs_changed; /* spring input data changed flag */
extern int spring_buffer_size; /* size of the spring constraint buffer */
extern int spring_lookup_size; /* size of the spring force lookup tables */
extern int dashpot_lookup_size; /* size of the dashpot force lookup tables */
extern int unload_lookup_size; /* size of the unload force lookup tables */
extern void spring_buffer_grow (int spring_lookup, int dashpot_lookup, int unload_lookup); /* grow buffer */

extern int cnsnum; /* number of constraints */
extern int *cnspart; /* constrained particle numbers */
extern REAL *cnslin[9]; /* constrained linear directions */
extern REAL *cnsang[9]; /* constrained angular directions */
extern int constrain_buffer_size; /* size of constrained particles buffer */
extern int constrain_buffer_grow (); /* grow buffer */

extern int tmsnum; /* number of time series */
extern pointer_t *tms; /* time series */
extern int time_series_buffer_size; /* size of time series buffer */
extern int time_series_buffer_grow (); /* grow buffer */

extern int prsnum; /* number of particles with prescribed motion */
extern int *prspart; /* prescribed motion particle numbers */
extern pointer_t *prslin; /* prescribed linear motion time history callbacks */
extern int *tmslin[3]; /* prescribed linear motion time series */
extern int *linkind; /* prescribied linear motion signal kind: 0-velocity, 1-acceleration */
extern pointer_t *prsang; /* prescribed angular motion time history callbacks */
extern int *tmsang[3]; /* prescribed angular motion time series */
extern int *angkind; /* prescribied angular motion signal kind: 0-velocity, 1-acceleration */
extern int prescribe_buffer_size; /* size of prescribed particle motion buffer */
extern int prescribe_buffer_grow (); /* grow buffer */

extern int hisnum; /* number of time histories */
extern int *hislst; /* history source lists */
extern int *hisidx; /* history source list start index */
extern int *hisent; /* history entity */
extern int *hiskind; /* history kind */
extern REAL *source[6]; /* source sphere or box definition or optional point */
extern pointer_t *history; /* Python list storing history */
extern int history_buffer_size; /* size of history buffer */
extern int history_list_size; /* size of history particle lists buffer */
extern void history_buffer_grow (int list_size); /* grow buffer */

extern int outnum; /* number of output lists */
extern int *outmode; /* output mode */
extern int *outpart; /* output particle lists */
extern int *outidx; /* output particle list start index */
extern int *outent; /* output entities per output mode */
extern int outrest[2]; /* 0: default output entities for unlisted particles and, 1: default output mode */
extern int outformat; /* output format */
extern int output_buffer_size; /* size of output buffer */
extern int output_list_size; /* size of output particle lists buffer */
extern void output_buffer_grow (int list_size); /* grow buffer */

extern REAL gravity[3]; /* gravity vector */
extern pointer_t gravfunc[3]; /* gravity callbacks */
extern int gravtms[3]; /* gravity time series */

extern REAL damping[6]; /* linear and angular damping */
extern pointer_t lindamp; /* linead damping callback */
extern int lindamptms[3]; /* linear damping time series */
extern pointer_t angdamp; /* angular damping callback */
extern int angdamptms[3]; /* angular damping time series */

struct prescribed_body_force /* externally prescribed body force */
{
  int particle;
  REAL inner_force[3];
  REAL inner_torque[3];
  REAL outer_force[3];
  REAL outer_torque[3];
};

extern MAP *prescribed_body_forces; /* particle index based map of prescibed body forces */

extern void declare_analytical (int k); /* declare particle 'k' analytical */

/**************** library interface ****************/

void init (); /* init memory */

void reset (); /* reset all simulation data */

int input (const char *path); /* interpret an input file (return 0 on success) */

/* run DEM simulation (return timed duration) */
REAL dem (REAL duration, /* simulation duration */
          REAL step, /* time step */
	  REAL *interval, /* optional constant output intervals (two); when passed with interval_func or interval_tms it stores current values */
	  pointer_t *interval_func, /* optional output interval callbacks (two) */
          int *interval_tms, /* optional output interval TSERIES numbers (two) */
	  char *prefix, /* optional output directory prefix */
	  int verbose, /* verbosity flag; 0 disables verbose output */
	  double adaptive); /* adaptive time stepping ratio; 0.0 disables adaptive time stepping */

#ifdef __cplusplus
} /* namespace */
#endif
#endif

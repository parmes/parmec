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

#include "macros.h"
#include "condet.h"
#include "constants.h"

/* zero force and torque */
task void zero_force_and_torque (uniform int span, uniform int parnum, uniform REAL * uniform force[3], uniform REAL * uniform torque[3])
{
  uniform int start = taskIndex*span;
  uniform int end = taskIndex == taskCount-1 ? parnum: start+span;

  foreach (i = start ... end)
  {
    force[0][i] = 0.0;
    force[1][i] = 0.0;
    force[2][i] = 0.0;
    torque[0][i] = 0.0;
    torque[1][i] = 0.0;
    torque[2][i] = 0.0;
  }
}

/* clear slave contact points */
task void clear_slaves (uniform int span, uniform slave_conpnt slave[], uniform int parnum)
{
  uniform int start = taskIndex*span;
  uniform int end = taskIndex == taskCount-1 ? parnum: start+span;

  for (uniform int i = start; i < end; i ++)
  {
    uniform slave_conpnt * uniform ptr = slave[i].next;
    while (ptr)
    {
      uniform slave_conpnt * uniform next = ptr->next;
      delete ptr;
      ptr = next;
    }

    slave[i].size = 0;
    slave[i].next = NULL;
  }
}

/* allocate new slave contact point that can be written to */
static uniform slave_conpnt * uniform newcon (uniform slave_conpnt * uniform slave, uniform int taskindex, uniform int *k)
{
  while (slave->lock != taskindex)
    atomic_compare_exchange_global (&slave->lock, -1, taskindex); /* lock access */

  uniform slave_conpnt * uniform con = slave;

  while (con->size == CONBUF) con = con->next; /* rewind to the end */

  *k = con->size ++;

  if (con->size == CONBUF)
  {
    uniform slave_conpnt * uniform ptr = uniform new uniform slave_conpnt;
    ptr->size = 0;
    ptr->next = NULL;
    con->next = ptr; /* append new item at the end */
  }

  slave->lock = -1; /* unlock */

  return con;
}

/* return pairing index based on (i,j) pairing of colors */
static inline int pairing (uniform int nummat, uniform int pairs[], int i, int j)
{
  int p[2] = {i, j}, start = 1, end = nummat;

  if (i > j)
  {
    p[0] = j;
    p[1] = i;
  }

  while (start < end)
  {
    int mid = ((end-start)>>1) + start;
    if (p[0] > pairs[2*mid]) start = mid;
    else if(p[0] == pairs[2*mid])
    {
      if (p[1] > pairs[2*mid+1]) start = mid;
      else if (p[1] == pairs[2*mid+1]) return mid;
      else end = mid;
    }
    else end = mid;
  }

  return 0; /* default material */
}

/* granular material force */
inline static int granural_force (REAL n[3], REAL vij[3], REAL oij[3], REAL depth, uniform int i, int j,
  uniform REAL mass[], uniform REAL invm[], uniform REAL * uniform iparam[NIPARAM], int ij, REAL state[NSTATE],
  uniform REAL step, REAL f[3])
{
  REAL ma;

  if (j >= 0) ma = 1.0 / (invm[i] + invm[j]);
  else ma = mass[i];

  /* normal force */

  REAL kn = iparam[SPRING][ij];
  REAL en = iparam[DAMPER][ij] * 2.0 * sqrt(kn*ma);
  REAL vn = DOT(vij,n);
  REAL fn = kn*depth + en*vn;

  f[0] = fn*n[0];
  f[1] = fn*n[1];
  f[2] = fn*n[2];

  /* sliding friction */

  REAL kt = kn * iparam[KSKN][ij];
  REAL et = en * iparam[KSKN][ij];
  REAL mus = iparam[FRISTAT][ij];
  REAL mud = iparam[FRIDYN][ij];

  REAL vt[3], ft[3];

  vt[0] = vij[0] - vn*n[0];
  vt[1] = vij[1] - vn*n[1];
  vt[2] = vij[2] - vn*n[2];

  state[0] += step*vt[0];
  state[1] += step*vt[1];
  state[2] += step*vt[2];

  ft[0] = kt*state[0] + et*vt[0];
  ft[1] = kt*state[1] + et*vt[1];
  ft[2] = kt*state[2] + et*vt[2];

  REAL len = LEN(ft);
  REAL lim = mus*fn;
  if (len > lim)
  {
    lim = mud*fn/len;
    ft[0] *= lim;
    ft[1] *= lim;
    ft[2] *= lim;
  }

  f[0] += ft[0];
  f[1] += ft[1];
  f[2] += ft[2];

  /* rolling and drilling friction */

  REAL mur = iparam[FRIROL][ij];
  REAL mul = iparam[FRIDRIL][ij];

  /* TODO */

  return depth < 0.0 ? 1 : 0;
}

/* contact forces task */
task void contacts_task (uniform int span, uniform master_conpnt master[], uniform slave_conpnt slave[],
  uniform int parnum, uniform REAL * uniform angular[6], uniform REAL * uniform linear[3],
  uniform REAL * uniform rotation[9], uniform REAL * uniform position[3], uniform REAL * uniform inertia[9],
  uniform REAL * uniform inverse[9], uniform REAL mass[], uniform REAL invm[], uniform REAL obspnt[],
  uniform REAL obslin[], uniform REAL obsang[], uniform int parmat[], uniform REAL * uniform mparam[NMAT],
  uniform int pairnum, uniform int pairs[], uniform int ikind[], uniform REAL * uniform iparam[NIPARAM],
  uniform REAL step)
{
  uniform int start = taskIndex*span;
  uniform int end = taskIndex == taskCount-1 ? parnum: start+span;

  for (uniform int i = start; i < end; i ++)
  {
    uniform REAL oi[3], v[3], x[3];

    oi[0] = angular[3][i];
    oi[1] = angular[4][i];
    oi[2] = angular[5][i];

    v[0] = linear[0][i];
    v[1] = linear[1][i];
    v[2] = linear[2][i];

    x[0] = position[0][i];
    x[1] = position[1][i];
    x[2] = position[2][i];

    /* update contact forces */
    for (uniform master_conpnt * uniform con = &master[i]; con; con = con->next)
    {
      uniform int gone[CONBUF];

      foreach (k = 0 ... con->size)
      {
	REAL p[3], n[3], z[3], vi[3], vj[3], oj[3], vij[3], oij[3];

	p[0] = con->point[0][k];
	p[1] = con->point[1][k];
	p[2] = con->point[2][k];

	n[0] = con->normal[0][k];
	n[1] = con->normal[1][k];
	n[2] = con->normal[2][k];

	z[0] = p[0]-x[0];
	z[1] = p[1]-x[1];
	z[2] = p[2]-x[2];

	vi[0] = oi[1]*z[2]-oi[2]*z[1] + v[0];
	vi[1] = oi[2]*z[0]-oi[0]*z[2] + v[1];
	vi[2] = oi[0]*z[1]-oi[1]*z[0] + v[2];

	int j = con->slave[0][k];

	if (j >= 0) /* particle-particle */
	{
	  z[0] = p[0]-position[0][j];
	  z[1] = p[1]-position[1][j];
	  z[2] = p[2]-position[2][j];

	  oj[0] = angular[3][j];
	  oj[1] = angular[4][j];
	  oj[2] = angular[5][j];

	  vj[0] = oj[1]*z[2]-oj[2]*z[1] + linear[0][j];
	  vj[1] = oj[2]*z[0]-oj[0]*z[2] + linear[1][j];
	  vj[2] = oj[0]*z[1]-oj[1]*z[0] + linear[2][j];
	}
	else /* particle-obstacle */
	{
	  int l = -j-2; /* obstacle index, see condet.ispc:drop_triangle and input.cpp:OBSTACLE */

	  if (l >= 0) /* moving obstacle */
	  {
	    z[0] = p[0]-obspnt[3*l];
	    z[1] = p[1]-obspnt[3*l+1];
	    z[2] = p[2]-obspnt[3*l+2];

	    oj[0] = obsang[3*l];
	    oj[1] = obsang[3*l+1];
	    oj[2] = obsang[3*l+2];

	    vj[0] = oj[1]*z[2]-oj[2]*z[1] + obslin[3*l];
	    vj[1] = oj[2]*z[0]-oj[0]*z[2] + obslin[3*l+1];
	    vj[2] = oj[0]*z[1]-oj[1]*z[0] + obslin[3*l+2];
	  }
	  else /* static boundary */
	  {
	    SET (oj, 0.0);
	    SET (vj, 0.0);
	  }
	}

	SUB (vj, vi, vij); /* relative linear velocity */
	SUB (oj, oi, oij); /* relative angular velocity */

	int ij = pairing (pairnum, pairs, con->color[0][k], con->color[1][k]);

	REAL state[NSTATE], f[3];

	for (uniform int l = 0; l < NSTATE; l ++)
	{
	  state[l] = con->state[l][k];
	}

	switch (ikind[ij])
	{
	case GRANULAR_FORCE:
	  gone[k] = granural_force (n, vij, oij, con->depth[k], i, j, mass, invm, iparam, ij, state, step, f);
	break;
	case USER_FORCE:
	  assert (0);
	  /* TODO */
	break;
	default:
	  print ("ERROR: invalid pairing kind");
	  assert (0);
	break;
	}

	con->force[0][k] = f[0];
	con->force[1][k] = f[1];
	con->force[2][k] = f[2];

	for (uniform int l = 0; l < NSTATE; l ++)
	{
	  con->state[l][k] = state[l];
	}
      }

      uniform int ngone = 0;

      for (uniform int k = 0; k < con->size; k ++)
      {
	if (gone[k] != 0)
	{
	  uniform int j = k+1;

	  while (j < con->size && gone[j] != 0) j ++;

	  if (j < con->size)
	  {
	    con->master[k] = con->master[j];
	    con->slave[0][k] = con->slave[0][j];
	    con->slave[1][k] = con->slave[1][j];
	    con->color[0][k] = con->color[0][j];
	    con->color[1][k] = con->color[1][j];
	    con->point[0][k] = con->point[0][j];
	    con->point[1][k] = con->point[1][j];
	    con->point[2][k] = con->point[2][j];
	    con->normal[0][k] = con->normal[0][j];
	    con->normal[1][k] = con->normal[1][j];
	    con->normal[2][k] = con->normal[2][j];
	    con->depth[k] = con->depth[j];
	    con->force[0][k] = con->force[0][j];
	    con->force[1][k] = con->force[1][j];
	    con->force[2][k] = con->force[2][j];

	    gone[j] = -1; /* not to be used again */
	  }

	  if (gone[k] > 0) ngone ++;
	}
      }

      con->size -= ngone;
    }

    uniform master_conpnt * uniform con = master[i].next;

    while (con && con->next) /* delete empty items */
    {
      uniform master_conpnt * uniform next = con->next;

      if (next->size == 0)
      {
	con->next = next->next;

	delete next;
      }

      con = con->next;
    }

    /* symmetrical copy into slave contact points */
    for (uniform master_conpnt * uniform con = &master[i]; con; con = con->next)
    {
      for (uniform int j = 0; j < con->size; j ++)
      {
	uniform slave_conpnt *ptr;
	uniform int k;

	if (con->slave[0][j] >= 0) /* particle-particle contact */
	{
	  ptr = newcon (&slave[con->slave[0][j]], taskIndex, &k);

	  ptr->master[0][k] = i;
	  ptr->master[1][k] = con->master[j];
	  ptr->point[0][k] = con->point[0][j];
	  ptr->point[1][k] = con->point[1][j];
	  ptr->point[2][k] = con->point[2][j];
	  ptr->force[0][k] = -con->force[0][j];
	  ptr->force[1][k] = -con->force[1][j];
	  ptr->force[2][k] = -con->force[2][j];
	}
      }
    }
  }
}

/* look up table based linear spline interpolation */
static inline uniform REAL look_up (uniform REAL * uniform tab[2], uniform int lo, uniform int hi, uniform REAL arg)
{
  while (hi > lo+2)
  {
    uniform int mid = lo + (hi-lo)/2;
    if (arg > tab[0][mid]) lo = mid;
    else hi = mid;
  }

  return tab[1][lo] + (tab[1][lo+1]-tab[1][lo])*(arg-tab[0][lo])/(tab[0][lo+1]-tab[0][lo]);
}

/* look up table based linear spline slope based on (x,y)-search (i == 0 -> x; i == 1 -> y) */
static inline uniform REAL look_up_slope (uniform int i, uniform REAL * uniform tab[2], uniform int lo, uniform int hi, uniform REAL arg)
{
  while (hi > lo+2)
  {
    uniform int mid = lo + (hi-lo)/2;
    if (arg > tab[i][mid]) lo = mid;
    else hi = mid;
  }

  return (tab[1][lo+1]-tab[1][lo])/(tab[0][lo+1]-tab[0][lo]);
}

/* contact (and gravity) forces accumulation task */
task void contacts_acc_task (uniform int span, uniform master_conpnt master[], uniform slave_conpnt slave[],
  uniform int parnum, uniform REAL * uniform position[6], uniform REAL mass[], uniform REAL gravity[3],
  uniform REAL * uniform force[3], uniform REAL * uniform torque[3])
{
  uniform int start = taskIndex*span;
  uniform int end = taskIndex == taskCount-1 ? parnum: start+span;

  for (uniform int i = start; i < end; i ++) /* force accumulation */
  {
    REAL f[3], a[3], fs[3], ts[3];
    uniform REAL po[3], ma;

    po[0] = position[0][i];
    po[1] = position[1][i];
    po[2] = position[2][i];

    ma = mass[i];

    SET (fs, 0.0);
    SET (ts, 0.0);

    for (uniform master_conpnt * uniform m = &master[i]; m; m = m->next)
    {
      foreach (j = 0 ... m->size)
      {
	f[0] = m->force[0][j];
	f[1] = m->force[1][j];
	f[2] = m->force[2][j];

	a[0] = m->point[0][j]-po[0];
	a[1] = m->point[1][j]-po[1];
	a[2] = m->point[2][j]-po[2];

	ACC (f, fs);
	PRODUCTADD (a, f, ts);
      }
    }

    for (uniform slave_conpnt * uniform s = &slave[i]; s; s = s->next)
    {
      foreach (j = 0 ... s->size)
      {
	f[0] = s->force[0][j];
	f[1] = s->force[1][j];
	f[2] = s->force[2][j];

	a[0] = s->point[0][j]-po[0];
	a[1] = s->point[1][j]-po[1];
	a[2] = s->point[2][j]-po[2];

	ACC (f, fs);
	PRODUCTADD (a, f, ts);
      }
    }

    force[0][i] = reduce_add (fs[0]) + ma * gravity[0];
    force[1][i] = reduce_add (fs[1]) + ma * gravity[1];
    force[2][i] = reduce_add (fs[2]) + ma * gravity[2];

    torque[0][i] = reduce_add (ts[0]);
    torque[1][i] = reduce_add (ts[1]);
    torque[2][i] = reduce_add (ts[2]);
  }
}

/* update spring foces */
task void springs_task (uniform int span, uniform int sprnum, uniform int sprtype[], uniform int * uniform sprpart[2],
  uniform REAL * uniform sprpnt[2][6], uniform REAL * uniform spring[2], uniform int spridx[],
  uniform REAL * uniform dashpot[2], uniform int dashidx[], uniform REAL * uniform unload[2], uniform int unidx[],
  uniform REAL * uniform yield[2], uniform REAL * uniform sprdir[3],
  uniform int sprflg[], uniform REAL stroke0[], uniform REAL * uniform stroke[2], uniform REAL * uniform sprfrc[2],
  uniform REAL * uniform angular[6], uniform REAL * uniform linear[3], uniform REAL * uniform rotation[9],
  uniform REAL * uniform position[6], uniform REAL * uniform force[3], uniform REAL * uniform torque[3],
  uniform int lock[])
{
  uniform int start = taskIndex*span;
  uniform int end = taskIndex == taskCount-1 ? sprnum: start+span;

  uniform REAL oj[3], vj[3], xj[3], Xj[3], Lj[9], ok[3], vk[3], xk[3], Xk[3], Lk[9];
  uniform int j0 = -1, k0 = -1;

  if (start > 0) /* start at sprpart[0][] change --> avoid atomic accumulation of force[][j] and torque[][j] below */
  {
    for (; start < end && sprpart[0][start] == sprpart[0][start-1]; start ++);
  }

  if (end < sprnum) /* end at sprpart[0][] change --> avoid atomic accumulation of force[][j] and torque[][j] below */
  {
    for (; end < sprnum-1 && sprpart[0][end] == sprpart[0][end+1]; end ++);
  }

  uniform REAL fj[3], tj[3], fk[3], tk[3]; /* local force and torque accumulation buffers */

  ZERO (fj);
  ZERO (tj);
  ZERO (fk);
  ZERO (tk);

  /* update spring points and spring forces; update global force and torque */
  for (uniform int i = start; i < end; i ++) /* XXX --> try foreach */
  {
    uniform REAL len, inv, velocity, spring_force, dashpot_force, total_force;
    uniform REAL Y[3], p[3], q[3], z[3], w[3], u[3], a[3], b[3];
    uniform int j = sprpart[0][i], k = sprpart[1][i];

    if (j0 != j) /* springs are sorted (j,k)-wise --> avoid repeated reads */
    {
      oj[0] = angular[3][j];
      oj[1] = angular[4][j];
      oj[2] = angular[5][j];

      vj[0] = linear[0][j];
      vj[1] = linear[1][j];
      vj[2] = linear[2][j];

      xj[0] = position[0][j];
      xj[1] = position[1][j];
      xj[2] = position[2][j];

      Xj[0] = position[3][j];
      Xj[1] = position[4][j];
      Xj[2] = position[5][j];

      Lj[0] = rotation[0][j];
      Lj[1] = rotation[1][j];
      Lj[2] = rotation[2][j];
      Lj[3] = rotation[3][j];
      Lj[4] = rotation[4][j];
      Lj[5] = rotation[5][j];
      Lj[6] = rotation[6][j];
      Lj[7] = rotation[7][j];
      Lj[8] = rotation[8][j];
    }

    Y[0] = sprpnt[0][3][i];
    Y[1] = sprpnt[0][4][i];
    Y[2] = sprpnt[0][5][i];

    SUB (Y, Xj, Y);
    NVADDMUL (xj, Lj, Y, p);

    sprpnt[0][0][i] = p[0];
    sprpnt[0][1][i] = p[1];
    sprpnt[0][2][i] = p[2];

    COPY (vj, u);
    SCALE(u, -1.0);
    SUB (p, xj, a);
    PRODUCTSUB (oj, a, u);

    if (k >= 0)
    {
      if (k0 != k)
      {
	ok[0] = angular[3][k];
	ok[1] = angular[4][k];
	ok[2] = angular[5][k];

	vk[0] = linear[0][k];
	vk[1] = linear[1][k];
	vk[2] = linear[2][k];

	xk[0] = position[0][k];
	xk[1] = position[1][k];
	xk[2] = position[2][k];

	Xk[0] = position[3][k];
	Xk[1] = position[4][k];
	Xk[2] = position[5][k];

	Lk[0] = rotation[0][k];
	Lk[1] = rotation[1][k];
	Lk[2] = rotation[2][k];
	Lk[3] = rotation[3][k];
	Lk[4] = rotation[4][k];
	Lk[5] = rotation[5][k];
	Lk[6] = rotation[6][k];
	Lk[7] = rotation[7][k];
	Lk[8] = rotation[8][k];
      }

      Y[0] = sprpnt[1][3][i];
      Y[1] = sprpnt[1][4][i];
      Y[2] = sprpnt[1][5][i];

      SUB (Y, Xk, Y);
      NVADDMUL (xk, Lk, Y, q);

      sprpnt[1][0][i] = q[0];
      sprpnt[1][1][i] = q[1];
      sprpnt[1][2][i] = q[2];

      ACC (vk, u);
      SUB (q, xk, b);
      PRODUCTADD (ok, b, u);
    }
    else
    {
      q[0] = sprpnt[1][0][i];
      q[1] = sprpnt[1][1][i];
      q[2] = sprpnt[1][2][i];
    }

    SUB (q, p, w);

    switch (sprflg[i]&(SPRDIR_FOLLOWER|SPRDIR_CONSTANT|SPRDIR_PLANAR))
    {
    case SPRDIR_FOLLOWER:

      len = LEN(w);
      inv = len > 0.0 ? 1.0/len : 0.0; /* zero force in case of zero length via (temporary) zero direction */

      z[0] = inv*w[0];
      z[1] = inv*w[1];
      z[2] = inv*w[2];

      sprdir[0][i] = z[0]; /* store current direction */
      sprdir[1][i] = z[1];
      sprdir[2][i] = z[2];
      break;

    case SPRDIR_CONSTANT:

      z[0] = sprdir[0][i];
      z[1] = sprdir[1][i];
      z[2] = sprdir[2][i];

      len = DOT (w, z);
      break;

    case SPRDIR_PLANAR:

      inv = w[0]*sprdir[0][i]+w[1]*sprdir[1][i]+w[2]*sprdir[2][i]; /* dot of (q-p) and constant direction */

      w[0] -= inv*sprdir[0][i];
      w[1] -= inv*sprdir[1][i];
      w[2] -= inv*sprdir[2][i];

      len = LEN(w);
      inv = len > 0.0 ? 1.0/len : 0.0; /* zero force in case of zero length via (temporary) zero direction */

      z[0] = inv*w[0];
      z[1] = inv*w[1];
      z[2] = inv*w[2];
      break;
    }

    switch (sprtype[i])
    {
    case SPRING_NONLINEAR_ELASTIC:
    {
      stroke[0][i] = len - stroke0[i];

      spring_force = look_up (spring, spridx[i], spridx[i+1], stroke[0][i]);
    }
    break;
    case SPRING_GENERAL_NONLINEAR:
    {
      if ((sprflg[i]&SPRING_YIELDED) == 0) /* elastic */
      {
	stroke[0][i] = len - stroke0[i];

	spring_force = look_up (spring, spridx[i], spridx[i+1], stroke[0][i]);

	if (spring_force < yield[0][i] || spring_force > yield[1][i])
	{
	  sprflg[i] |= SPRING_YIELDED;

	  stroke[1][i] = 0.0; /* plastic stroke is zero --> yielded for the first time */
	}
      }
      else if (sprflg[i]&SPRING_KINEMATIC) /* yielded --> kinematic hardening */
      {
	uniform REAL s0 = stroke[0][i];

	uniform REAL s1 = len - stroke0[i];

	uniform REAL ds = abs(s1)-abs(s0);

	uniform REAL f0 = sprfrc[1][i];

	stroke[0][i] = s1;

	uniform REAL df = 0.0;

	if (ds > 0) /* loading or re-loading */
	{
	  uniform REAL f1 = look_up (spring, spridx[i], spridx[i+1], s1-stroke[1][i]);

	  if (sprflg[i]&SPRING_UNLOADING)
	  {
            uniform REAL dfds = look_up_slope (1, unload, unidx[i], unidx[i+1], f0);

	    uniform REAL f2 = f0 + dfds*ds;

	    if (f2 > f1) /* back to loading */
	    {
	      df = f1-f0;

	      sprflg[i] &= ~SPRING_UNLOADING;
	    }
	    else /* re-loading */
	    {
	      df = f2-f0;
	    }
	  }
	}
	else /* unloading */
	{
          uniform REAL dfds = look_up_slope (1, unload, unidx[i], unidx[i+1], f0);

	  df = dfds*ds;

	  sprflg[i] |= SPRING_UNLOADING;
	}

	spring_force = f0 + df;

	if (spring_force * f0 <= 0.0) /* transition from compression into tension or vice versa */
	{
	  stroke[1][i] = s0 - f0*(s1 - s0)/df; /* record plastic stroke at zero crossing of (s0,f0)---(s1,f1) */
	}

	/* TODO --> completely untested */
      }
      else if (sprflg[i]&SPRING_ISOTROPIC) /* yielded --> isotropic hardening */
      {
	print ("ERROR: isotropic hardening has not been implemented --> using elastic spring\n");

	stroke[0][i] = len - stroke0[i];

	spring_force = look_up (spring, spridx[i], spridx[i+1], stroke[0][i]);

	/* TODO --> implement */
      }
      else /* error */
      {
	print ("ERROR: hardening rule not set for a general nonlinear spring --> using elastic spring\n");

	stroke[0][i] = len - stroke0[i];

	spring_force = look_up (spring, spridx[i], spridx[i+1], stroke[0][i]);
      }
    }
    break;
    }

    velocity = DOT (z, u);

    dashpot_force = look_up (dashpot, dashidx[i], dashidx[i+1], velocity);

    total_force = spring_force + (spring_force == 0.0 ? 0.0 : dashpot_force);

    SCALE (z, total_force);

    sprfrc[0][i] = total_force;
    sprfrc[1][i] = spring_force;

    PRODUCT (a, z, w);

    if (j0 != j && j0 >= 0) /* global acc previous and local reset current */
    {
      force[0][j0] += fj[0];
      force[1][j0] += fj[1];
      force[2][j0] += fj[2];
      torque[0][j0] += tj[0];
      torque[1][j0] += tj[1];
      torque[2][j0] += tj[2];

      COPY (z, fj);
      COPY (w, tj);
    }
    else /* local acc current */
    {
      ACC (z, fj);
      ACC (w, tj);
    }

    if (k >= 0)
    {
      PRODUCT (b, z, w);

      if (k0 != k && k0 >= 0) /* global acc previous and local reset current */
      {
	while (lock[k0] != taskIndex)
	  atomic_compare_exchange_global (&lock[k0], -1, taskIndex); /* lock access to k0 force */

	force[0][k0] -= fk[0];
	force[1][k0] -= fk[1];
	force[2][k0] -= fk[2];
	torque[0][k0] -= tk[0];
	torque[1][k0] -= tk[1];
	torque[2][k0] -= tk[2];

	lock[k0] = -1; /* unlock */

	COPY (z, fk);
	COPY (w, tk);
      }
      else /* local acc current */
      {
	ACC (z, fk);
	ACC (w, tk);
      }

      k0 = k;
    }

    j0 = j;
  }

  if (j0 >= 0) /* final global acc */
  {
    force[0][j0] += fj[0];
    force[1][j0] += fj[1];
    force[2][j0] += fj[2];
    torque[0][j0] += tj[0];
    torque[1][j0] += tj[1];
    torque[2][j0] += tj[2];
  }

  if (k0 >= 0) /* final global acc */
  {
    while (lock[k0] != taskIndex)
      atomic_compare_exchange_global (&lock[k0], -1, taskIndex); /* lock access to k0 force */

    force[0][k0] -= fk[0];
    force[1][k0] -= fk[1];
    force[2][k0] -= fk[2];
    torque[0][k0] -= tk[0];
    torque[1][k0] -= tk[1];
    torque[2][k0] -= tk[2];

    lock[k0] = -1; /* unlock */
  }
}

/* update forces */
export void forces (uniform int threads, uniform master_conpnt master[], uniform slave_conpnt slave[],
  uniform int parnum, uniform REAL * uniform angular[6], uniform REAL * uniform linear[3],
  uniform REAL * uniform rotation[9], uniform REAL * uniform position[3], uniform REAL * uniform inertia[9],
  uniform REAL * uniform inverse[9], uniform REAL mass[], uniform REAL invm[], uniform REAL obspnt[],
  uniform REAL obslin[], uniform REAL obsang[], uniform int parmat[], uniform REAL * uniform mparam[NMAT],
  uniform int pairnum, uniform int pairs[], uniform int ikind[], uniform REAL * uniform iparam[NIPARAM],
  uniform REAL step, uniform int sprnum, uniform int sprtype[], uniform int * uniform sprpart[2],
  uniform REAL * uniform sprpnt[2][6], uniform REAL * uniform spring[2], uniform int spridx[],
  uniform REAL * uniform dashpot[2], uniform int dashidx[], uniform REAL * uniform unload[2], uniform int unidx[], 
  uniform REAL * uniform yield[2], uniform REAL * uniform sprdir[3], uniform int sprflg[], uniform REAL stroke0[],
  uniform REAL * uniform stroke[2], uniform REAL * uniform sprfrc[2],
  uniform REAL gravity[3], uniform REAL * uniform force[3], uniform REAL * uniform torque[3])
{
  launch [threads] zero_force_and_torque (parnum/threads, parnum, force, torque);

  launch [threads] clear_slaves (parnum/threads, slave, parnum);
  sync;

  launch [threads] contacts_task (parnum/threads, master, slave, parnum, angular, linear, rotation, position, inertia, inverse,
                                  mass, invm, obspnt, obslin, obsang, parmat, mparam, pairnum, pairs, ikind, iparam, step);
  sync;

  launch [threads] contacts_acc_task (parnum/threads, master, slave, parnum, position, mass, gravity, force, torque);
  sync;

  uniform int * uniform lock = uniform new uniform int [parnum];
  foreach (i = 0 ... parnum) lock[i] = -1;

  launch [threads] springs_task (sprnum/threads, sprnum, sprtype, sprpart, sprpnt, spring, spridx, dashpot, dashidx, unload, unidx, yield,
                                 sprdir, sprflg, stroke0, stroke, sprfrc, angular, linear, rotation, position, force, torque, lock);
  sync;

  delete lock;
}

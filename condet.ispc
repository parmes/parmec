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
#include "partition.h"

/* sphere-ellipsoid contact */
inline static REAL sphere_ellipsoid (uniform REAL p[3], uniform REAL r, REAL center[3],
  REAL radii[3], REAL orient[9], uniform REAL point[3][LSIZE], uniform REAL normal[3][LSIZE], int j)
{
  assert (0); /* TODO */
  return 0.0;
}
inline static REAL sphere_ellipsoid (REAL p[3], REAL r, uniform REAL center[3],
  uniform REAL radii[3], uniform REAL orient[9], uniform REAL point[3][LSIZE], uniform REAL normal[3][LSIZE], int j)
{
  assert (0); /* TODO */
  return 0.0;
}

/* ellipsoid-ellipsoid contact */
inline static REAL ellipsoid_ellipsoid (uniform REAL p1[3], uniform REAL r1[3], uniform REAL o1[9],
  REAL p2[3], REAL r2[3], REAL o2[9], uniform REAL point[3][LSIZE], uniform REAL normal[3][LSIZE], int j)
{
  assert (0); /* TODO */
  return 0.0;
}

/* allocate new master contact point that can be written to */
static uniform master_conpnt * uniform newcon (uniform master_conpnt * uniform master, uniform int taskindex, uniform int *k)
{
  if (taskindex >= 0)
  {
    while (master->lock != taskindex)
      atomic_compare_exchange_global (&master->lock, -1, taskindex); /* lock access */
  }

  uniform master_conpnt * uniform con = master;

  while (con->size == CONBUF && con->next != NULL) con = con->next; /* find available item or rewind to end */

  if (con->size < CONBUF)
  {
    *k = con->size ++;
  }
  else
  {
    uniform master_conpnt * uniform ptr = uniform new uniform master_conpnt;
    ptr->size = 1; /* index zero already used up */
    ptr->next = NULL;
    con->next = ptr; /* append new item at the end */
    con = ptr; /* return new item */
    *k = 0;
  }

  if (taskindex >= 0) master->lock = -1; /* unlock */

  return con;
}

/* maximum of three numbers */
inline static uniform REAL max (uniform REAL a, uniform REAL b, uniform REAL c)
{
  uniform REAL d = (a > b ? a : b);
  return d > c ? d : c;
}

/* minimum of three numbers */
inline static uniform REAL min (uniform REAL a, uniform REAL b, uniform REAL c)
{
  uniform REAL d = (a < b ? a : b);
  return d < c ? d : c;
}

/* drop ellipsoid down the partitioning tree */
static void drop_ellipsoid (uniform partitioning tree[], uniform int node,
  uniform REAL lo[3], uniform REAL hi[3], uniform REAL rx,
  uniform REAL p[3], uniform REAL r[3], uniform REAL or[9],
  uniform int color, uniform int part, uniform int i,
  uniform master_conpnt master[], uniform int taskindex)
{
  uniform int d = tree[node].dimension;

  if (d >= 0) /* node */
  {
    if (lo[d] <= tree[node].coord)
      drop_ellipsoid (tree, tree[node].left, lo, hi, rx, p, r, or, color, part, i, master, taskindex);

    if (hi[d] >= tree[node].coord)
      drop_ellipsoid (tree, tree[node].right, lo, hi, rx, p, r, or, color, part, i, master, taskindex);
  }
  else /* leaf */
  {
    uniform leaf_data * uniform l = tree[node].data;

    uniform REAL point[3][LSIZE];
    uniform REAL normal[3][LSIZE];
    uniform REAL depth[LSIZE];

    if (r[1] < 0.) /* sphere- */
    {
      foreach (j = 0 ... l->size)
      {
	cif (l->radii[1][j] < 0.) /* sphere-sphere */
	{
	  REAL q[3], c[3], len, ilen;

	  c[0] = l->center[0][j];
	  c[1] = l->center[1][j];
	  c[2] = l->center[2][j];
	  q[0] = p[0]-c[0];
	  q[1] = p[1]-c[1];
	  q[2] = p[2]-c[2];
	  len = LEN(q);
	  ilen = len > 0.0 ? 1.0/len : 1.0; /* test with self is possible */
	  point[0][j] = 0.5*(p[0]+c[0]);
	  point[1][j] = 0.5*(p[1]+c[1]); 
	  point[2][j] = 0.5*(p[2]+c[2]); 
	  normal[0][j] = ilen*q[0];
	  normal[1][j] = ilen*q[1];
	  normal[2][j] = ilen*q[2];
	  depth[j] = rx+l->radii[0][j] - len;
	}
	else /* sphere-ellipsoid */
	{
	  REAL center[3] = {l->center[0][j], l->center[1][j], l->center[2][j]};
	  REAL radii[3] = {l->radii[0][j], l->radii[1][j], l->radii[2][j]};
	  REAL orient[9] = {l->orient[0][j], l->orient[1][j], l->orient[2][j],
			    l->orient[3][j], l->orient[4][j], l->orient[5][j],
			    l->orient[6][j], l->orient[7][j], l->orient[8][j]};

	  depth[j] = sphere_ellipsoid (p, rx, center, radii, orient, point, normal, j);
	}
      }
    }
    else /* ellipsoid- */
    {
      foreach (j = 0 ... l->size)
      {
	cif (l->radii[1][j] > 0.) /* ellipsoid-ellipsoid */
	{
	  REAL center[3] = {l->center[0][j], l->center[1][j], l->center[2][j]};
	  REAL radii[3] = {l->radii[0][j], l->radii[1][j], l->radii[2][j]};
	  REAL orient[9] = {l->orient[0][j], l->orient[1][j], l->orient[2][j],
			    l->orient[3][j], l->orient[4][j], l->orient[5][j],
		            l->orient[6][j], l->orient[7][j], l->orient[8][j]};

	  depth[j] = ellipsoid_ellipsoid (p, r, or, center, radii, orient, point, normal, j);
	}
	else /* ellipsoid-sphere */
	{
	  REAL center[3] = {l->center[0][j], l->center[1][j], l->center[2][j]};

	  depth[j] = sphere_ellipsoid (center, l->radii[0][j], p, r, or, point, normal, j);
	}
      }
    }

    for (uniform int j = 0; j < l->size; j ++)
    {
      if (depth[j] > 0.0 && part != l->part[j])
      {
        uniform int found = 0;

	uniform int index = min(part, l->part[j]); /* master contact points are indexed by smaller particle index */

	uniform master_conpnt * uniform mlist = &master[index];

	while (mlist->lock != taskindex)
	  atomic_compare_exchange_global (&mlist->lock, -1, taskindex); /* lock access */

	if (part == index) /* part is master */
	{
	  for (uniform master_conpnt * uniform con = mlist; con; con = con->next)
	  {
	    for (uniform int k = 0; k < con->size; k ++)
	    {
	      if (con->master[k] == i && con->slave[0][k] == l->part[j] && con->slave[1][k] == l->ell[j]) /* found existing contact point */
	      {
		found = 1;
		goto out;
	      }
	    }
	  }
	}
	else /* l->part is master */
	{
          for (uniform master_conpnt * uniform con = mlist; con; con = con->next)
	  {
	    for (uniform int k = 0; k < con->size; k ++)
	    {
	      if (con->master[k] == l->part[j] && con->slave[0][k] == part && con->slave[1][k] == i) /* found existing contact point */
	      {
		found = 1;
		goto out;
	      }
	    }
	  }
	}

out:
	if (found == 0) /* not found above so create new contact point */
	{
	  uniform master_conpnt * uniform con;
	  uniform int k;

	  con = newcon (mlist, -1, &k); /* non-locking alloc */

	  if (part == index) /* part is master */
	  {
	    con->master[k] = i;
	    con->slave[0][k] = l->part[j];
	    con->slave[1][k] = l->ell[j];
	    con->color[0][k] = color;
	    con->color[1][k] = l->color[j];
	    con->point[0][k] = point[0][j];
	    con->point[1][k] = point[1][j];
	    con->point[2][k] = point[2][j];
	    con->normal[0][k] = normal[0][j];
	    con->normal[1][k] = normal[1][j];
	    con->normal[2][k] = normal[2][j];
	    con->depth[k] = depth[j];
	  }
	  else /* l->part is master */
	  {
            con->master[k] = l->part[j];
	    con->slave[0][k] = part;
	    con->slave[1][k] = i;
	    con->color[0][k] = l->color[j];
	    con->color[1][k] = color;
	    con->point[0][k] = point[0][j];
	    con->point[1][k] = point[1][j];
	    con->point[2][k] = point[2][j];
	    con->normal[0][k] = -normal[0][j];
	    con->normal[1][k] = -normal[1][j];
	    con->normal[2][k] = -normal[2][j];
	    con->depth[k] = depth[j];
	  }
	}

	mlist->lock = -1; /* unlock */
      }
    }
  }
}

/* test ellipsoids against those stored in the tree */
task void test_ellipsoids (uniform int span, uniform partitioning tree[], uniform int ellnum, uniform int ellcol[],
  uniform int part[], uniform REAL * uniform center[6], uniform REAL * uniform radii[3], uniform REAL * uniform orient[18],
  uniform master_conpnt master[])
{
  uniform int start = taskIndex*span;
  uniform int end = taskIndex == taskCount-1 ? ellnum: start+span;

  for (uniform int i = start; i < end; i ++)
  {
    uniform REAL p[3] = {center[0][i], center[1][i], center[2][i]};
    uniform REAL r[3] = {radii[0][i], radii[1][i], radii[2][i]};
    uniform REAL rx = max (r[0], r[1], r[2]);
    uniform REAL lo[3] = {p[0]-rx, p[1]-rx, p[2]-rx};
    uniform REAL hi[3] = {p[0]+rx, p[1]+rx, p[2]+rx};
    uniform REAL or[9] = {orient[0][i], orient[1][i], orient[2][i],
                          orient[3][i], orient[4][i], orient[5][i],
			  orient[6][i], orient[7][i], orient[8][i]};

    drop_ellipsoid (tree, 0, lo, hi, rx, p, r, or, ellcol[i], part[i], i, master, taskIndex);
  }
}

#define clamp(a,b,c) ((a) < (b) ? (b) : (a) > (c) ? (c) : (a))
/* triangle-sphere contact based on http://www.gamedev.net/topic/552906-closest-point-on-triangle */
inline static REAL triangle_sphere (uniform REAL ax, uniform REAL ay, uniform REAL az,
  uniform REAL bx, uniform REAL by, uniform REAL bz, uniform REAL cx, uniform REAL cy, uniform REAL cz,
  REAL px, REAL py, REAL pz, REAL r, uniform REAL point[3][LSIZE], uniform REAL normal[3][LSIZE], int j)
{
  uniform REAL edge0[3] = {bx-ax, by-ay, bz-az};
  uniform REAL edge1[3] = {cx-ax, cy-ay, cz-az};
  REAL v0[3] = {ax-px, ay-py, az-pz};

  uniform REAL a = DOT(edge0,edge0);
  uniform REAL b = DOT(edge0,edge1);
  uniform REAL c = DOT(edge1,edge1);
  REAL d = DOT(edge0,v0);
  REAL e = DOT(edge1,v0);

  uniform REAL det = a*c - b*b;
  REAL s = b*e - c*d;
  REAL t = b*d - a*e;

  if (s + t < det)
  {
    if (s < 0.f)
    {
      if (t < 0.f)
      {
	if (d < 0.f)
	{
	  s = clamp(-d/a, 0.f, 1.f);
	  t = 0.f;
	}
	else
	{
	  s = 0.f;
	  t = clamp(-e/c, 0.f, 1.f);
	}
      }
      else
      {
	s = 0.f;
	t = clamp(-e/c, 0.f, 1.f);
      }
    }
    else if (t < 0.f)
    {
      s = clamp(-d/a, 0.f, 1.f);
      t = 0.f;
    }
    else
    {
      REAL invDet = 1.f / det;
      s *= invDet;
      t *= invDet;
    }
  }
  else
  {
    if (s < 0.f)
    {
      REAL tmp0 = b+d;
      REAL tmp1 = c+e;
      if (tmp1 > tmp0)
      {
	REAL numer = tmp1 - tmp0;
	REAL denom = a-2*b+c;
	s = clamp(numer/denom, 0.f, 1.f);
	t = 1-s;
      }
      else
      {
	t = clamp(-e/c, 0.f, 1.f);
	s = 0.f;
      }
    }
    else if (t < 0.f)
    {
      if (a+d > b+e)
      {
	REAL numer = c+e-b-d;
	REAL denom = a-2*b+c;
	s = clamp(numer/denom, 0.f, 1.f);
	t = 1-s;
      }
      else
      {
	s = clamp(-e/c, 0.f, 1.f);
	t = 0.f;
      }
    }
    else
    {
      REAL numer = c+e-b-d;
      REAL denom = a-2*b+c;
      s = clamp(numer/denom, 0.f, 1.f);
      t = 1.f - s;
    }
  }

  point[0][j] = ax + s*edge0[0] + t*edge1[0];
  point[1][j] = ay + s*edge0[1] + t*edge1[1];
  point[2][j] = az + s*edge0[2] + t*edge1[2];
  REAL q[3] = {px-point[0][j], py-point[1][j], pz-point[2][j]};
  REAL len = LEN (q);
  REAL ilen = 1.0/len;
  normal[0][j] = ilen*q[0];
  normal[1][j] = ilen*q[1];
  normal[2][j] = ilen*q[2];
  return r - len;
}

/* triangle-ellipsoid contact */
inline static REAL triangle_ellipsoid (uniform REAL ax, uniform REAL ay, uniform REAL az,
  uniform REAL bx, uniform REAL by, uniform REAL bz, uniform REAL cx, uniform REAL cy, uniform REAL cz,
  REAL p[3], REAL r[3], REAL o[9], uniform REAL point[3][LSIZE], uniform REAL normal[3][LSIZE], int j)
{
  /* TODO */
  return 0.0;
}

/* drop triangle down the partitioning tree */
static void drop_triangle (uniform partitioning tree[], uniform int node, uniform REAL lo[3], uniform REAL hi[3],
  uniform REAL ax, uniform REAL ay, uniform REAL az, uniform REAL bx, uniform REAL by, uniform REAL bz, uniform REAL cx,
  uniform REAL cy, uniform REAL cz, uniform int color, uniform int triobs, uniform int i, uniform master_conpnt master[],
  uniform int taskindex)
{
  uniform int d = tree[node].dimension;

  if (d >= 0) /* node */
  {
    if (lo[d] <= tree[node].coord)
      drop_triangle (tree, tree[node].left, lo, hi, ax, ay, az, bx, by, bz, cx, cy, cz, color, triobs, i, master, taskindex);

    if (hi[d] >= tree[node].coord)
      drop_triangle (tree, tree[node].right, lo, hi, ax, ay, az, bx, by, bz, cx, cy, cz, color, triobs, i, master, taskindex);
  }
  else /* leaf */
  {
    uniform leaf_data * uniform l = tree[node].data;

    uniform REAL point[3][LSIZE];
    uniform REAL normal[3][LSIZE];
    uniform REAL depth[LSIZE];

    foreach (j = 0 ... l->size)
    {
      cif (l->radii[1][j] < 0.) /* triangle-sphere */
      {
	depth[j] = triangle_sphere (ax, ay, az, bx, by, bz, cx, cy, cz,
	  l->center[0][j], l->center[1][j], l->center[2][j], l->radii[0][j], point, normal, j);
      }
      else /* triangle-ellipsoid */
      {
	REAL center[3] = {l->center[0][j], l->center[1][j], l->center[2][j]};
	REAL radii[3] = {l->radii[0][j], l->radii[1][j], l->radii[2][j]};
	REAL orient[9] = {l->orient[0][j], l->orient[1][j], l->orient[2][j],
			  l->orient[3][j], l->orient[4][j], l->orient[5][j],
			  l->orient[6][j], l->orient[7][j], l->orient[8][j]};

	depth[j] = triangle_ellipsoid (ax, ay, az, bx, by, bz, cx, cy, cz, center, radii, orient, point, normal, j);
      }
    }

    for (uniform int j = 0; j < l->size; j ++)
    {
      if (depth[j] > 0.0)
      {
        uniform int found = 0;

	for (uniform master_conpnt * uniform con = &master[l->part[j]]; con; con = con->next)
	{
	  for (uniform int k = 0; k < con->size; k ++)
	  {
	    if (con->master[k] == l->ell[j] && con->slave[1][k] == -(i+1)) /* found existing contact point */
	    {
	      found = 1;
	      goto out;
	    }
	  }
	}
out:
	if (found == 0) /* create new contact point */
	{
	  uniform master_conpnt * uniform con;
	  uniform int k;

	  con = newcon (&master[l->part[j]], taskindex, &k);

	  con->master[k] = l->ell[j];
	  con->slave[0][k] = triobs; /* see input.cpp:OBSTACLE */
	  con->slave[1][k] = -(i+1); /* particle-triangle */ 
	  con->color[0][k] = l->color[j]; 
	  con->color[1][k] = color;
	  con->point[0][k] = point[0][j];
	  con->point[1][k] = point[1][j];
	  con->point[2][k] = point[2][j];
	  con->normal[0][k] = normal[0][j];
	  con->normal[1][k] = normal[1][j];
	  con->normal[2][k] = normal[2][j];
	  con->depth[k] = depth[j];
	}
      }
    }
  }
}

/* test triangles against ellipsoids stored in the tree */
task void test_triangles (uniform int span, uniform partitioning tree[],
  uniform int trinum, uniform int tricol[], uniform int triobs[],
  uniform REAL * uniform tri[3][3], uniform master_conpnt master[])
{
  uniform int start = taskIndex*span;
  uniform int end = taskIndex == taskCount-1 ? trinum: start+span;

  for (uniform int i = start; i < end; i ++)
  {
    uniform REAL ax=tri[0][0][i], ay=tri[0][1][i], az=tri[0][2][i],
                 bx=tri[1][0][i], by=tri[1][1][i], bz=tri[1][2][i],
		 cx=tri[2][0][i], cy=tri[2][1][i], cz=tri[2][2][i];
    uniform REAL lo[3] = {min(ax,bx,cx), min(ay,by,cy), min(az,bz,cz)};
    uniform REAL hi[3] = {max(ax,bx,cx), max(ay,by,cy), max(az,bz,cz)};

    drop_triangle (tree, 0, lo, hi, ax, ay, az, bx, by, bz, cx, cy, cz, tricol[i], triobs[i], i, master,
                   taskCount+taskIndex); /* shift task index not to verlap with drop_ellipsoids */
  }
}

/* triangle-sphere contact update based on http://www.gamedev.net/topic/552906-closest-point-on-triangle */
inline static REAL triangle_sphere (REAL ax, REAL ay, REAL az, REAL bx, REAL by, REAL bz, REAL cx, REAL cy, REAL cz,
  REAL px, REAL py, REAL pz, REAL r, uniform REAL point[3][CONBUF], uniform REAL normal[3][CONBUF], int j)
{
  REAL edge0[3] = {bx-ax, by-ay, bz-az};
  REAL edge1[3] = {cx-ax, cy-ay, cz-az};
  REAL v0[3] = {ax-px, ay-py, az-pz};

  REAL a = DOT(edge0,edge0);
  REAL b = DOT(edge0,edge1);
  REAL c = DOT(edge1,edge1);
  REAL d = DOT(edge0,v0);
  REAL e = DOT(edge1,v0);

  REAL det = a*c - b*b;
  REAL s = b*e - c*d;
  REAL t = b*d - a*e;

  if (s + t < det)
  {
    if (s < 0.f)
    {
      if (t < 0.f)
      {
	if (d < 0.f)
	{
	  s = clamp(-d/a, 0.f, 1.f);
	  t = 0.f;
	}
	else
	{
	  s = 0.f;
	  t = clamp(-e/c, 0.f, 1.f);
	}
      }
      else
      {
	s = 0.f;
	t = clamp(-e/c, 0.f, 1.f);
      }
    }
    else if (t < 0.f)
    {
      s = clamp(-d/a, 0.f, 1.f);
      t = 0.f;
    }
    else
    {
      REAL invDet = 1.f / det;
      s *= invDet;
      t *= invDet;
    }
  }
  else
  {
    if (s < 0.f)
    {
      REAL tmp0 = b+d;
      REAL tmp1 = c+e;
      if (tmp1 > tmp0)
      {
	REAL numer = tmp1 - tmp0;
	REAL denom = a-2*b+c;
	s = clamp(numer/denom, 0.f, 1.f);
	t = 1-s;
      }
      else
      {
	t = clamp(-e/c, 0.f, 1.f);
	s = 0.f;
      }
    }
    else if (t < 0.f)
    {
      if (a+d > b+e)
      {
	REAL numer = c+e-b-d;
	REAL denom = a-2*b+c;
	s = clamp(numer/denom, 0.f, 1.f);
	t = 1-s;
      }
      else
      {
	s = clamp(-e/c, 0.f, 1.f);
	t = 0.f;
      }
    }
    else
    {
      REAL numer = c+e-b-d;
      REAL denom = a-2*b+c;
      s = clamp(numer/denom, 0.f, 1.f);
      t = 1.f - s;
    }
  }

  point[0][j] = ax + s*edge0[0] + t*edge1[0];
  point[1][j] = ay + s*edge0[1] + t*edge1[1];
  point[2][j] = az + s*edge0[2] + t*edge1[2];
  REAL q[3] = {px-point[0][j], py-point[1][j], pz-point[2][j]};
  REAL len = LEN (q);
  REAL ilen = 1.0/len;
  normal[0][j] = ilen*q[0];
  normal[1][j] = ilen*q[1];
  normal[2][j] = ilen*q[2];
  return r - len;
}

/* update existing contact points */
task void update_existing (uniform int span, uniform int parnum, uniform master_conpnt master[],
  uniform REAL * uniform center[6], uniform REAL * uniform radii[3],
  uniform REAL * uniform orient[18], uniform REAL * uniform tri[3][3])
{
  uniform int start = taskIndex*span;
  uniform int end = taskIndex == taskCount-1 ? parnum: start+span;

  for (uniform int l = start; l < end; l ++)
  {
    for (uniform master_conpnt * uniform con = &master[l]; con; con = con->next)
    {
      foreach (k = 0 ... con->size)
      {
	int i = con->master[k];
	int j = con->slave[1][k];

	if (radii[1][i] < 0.0) /* sphere- */
	{
	  if (j < 0) /* sphere-triangle */
	  {
	    int u = -j-1;
	    con->depth[k] = triangle_sphere (tri[0][0][u], tri[0][1][u], tri[0][2][u],
	                                     tri[1][0][u], tri[1][1][u], tri[1][2][u],
					     tri[2][0][u], tri[2][1][u], tri[2][2][u],
					     center[0][i], center[1][i], center[2][i],
					     radii[0][i], con->point, con->normal, k);
	  }
	  else if (radii[1][j] < 0.0) /* sphere-sphere */
	  {
	    REAL p[3], q[3], c[3], len, ilen;

	    p[0] = center[0][i];
	    p[1] = center[1][i];
	    p[2] = center[2][i];
	    c[0] = center[0][j];
	    c[1] = center[1][j];
	    c[2] = center[2][j];
	    q[0] = p[0]-c[0];
	    q[1] = p[1]-c[1];
	    q[2] = p[2]-c[2];
	    len = LEN(q);
	    ilen = 1.0/len;
	    con->point[0][k] = 0.5*(p[0]+c[0]);
	    con->point[1][k] = 0.5*(p[1]+c[1]); 
	    con->point[2][k] = 0.5*(p[2]+c[2]); 
	    con->normal[0][k] = ilen*q[0];
	    con->normal[1][k] = ilen*q[1];
	    con->normal[2][k] = ilen*q[2];
	    con->depth[k] = radii[0][i]+radii[0][j] - len;
	  }
	  else /* sphere-ellipsoid */
	  {
	    assert (0); /* TODO */
	  }
	}
	else /* ellipsoid- */
	{
	  if (j < 0) /* ellipsoid-triangle */
	  {
	    assert (0); /* TODO */
	  }
	  else if (radii[1][j] < 0.0) /* ellipsoid-sphere */
	  {
	    assert (0); /* TODO */
	  }
	  else /* ellipsoid-ellipsoid */
	  {
	    assert (0); /* TODO */
	  }
	}
      }
    }
  }
}

/* allocate global array of master contact points */
export uniform master_conpnt * uniform master_alloc (uniform master_conpnt * uniform old, uniform int nold, uniform int size)
{
  uniform master_conpnt * uniform con = uniform new uniform master_conpnt [size];

  if (nold)
  {
    memcpy (con, old, nold * sizeof (uniform master_conpnt));

    delete old;
  }

  for (uniform int i = nold; i < size; i ++)
  {
    con[i].size = 0;
    con[i].next = NULL;
    con[i].lock = -1;
  }
 
  return con;
}

/* free global array of master contact points */
export void master_free (uniform master_conpnt * uniform con, uniform int size)
{
  for (uniform int i = 0; i < size; i ++)
  {
    uniform master_conpnt * uniform ptr = con[i].next;
    while (ptr)
    {
      uniform master_conpnt * uniform next = ptr->next;
      delete ptr;
      ptr = next;
    }
  }

  delete con;
}

/* allocate global array of slave contact points */
export uniform slave_conpnt * uniform slave_alloc (uniform slave_conpnt * uniform old, uniform int nold, uniform int size)
{
  uniform slave_conpnt * uniform con = uniform new uniform slave_conpnt [size];

  if (nold)
  {
    memcpy (con, old, nold * sizeof (uniform slave_conpnt));

    delete old;
  }

  for (uniform int i = 0; i < size; i ++)
  {
    con[i].size = 0;
    con[i].next = NULL;
    con[i].lock = -1;
  }
 
  return con;
}

/* free global array of slave contact points */
export void slave_free (uniform slave_conpnt * uniform con, uniform int size)
{
  for (uniform int i = 0; i < size; i ++)
  {
    uniform slave_conpnt * uniform ptr = con[i].next;
    while (ptr)
    {
      uniform slave_conpnt * uniform next = ptr->next;
      delete ptr;
      ptr = next;
    }
  }

  delete con;
}

/* perform contact detection */
export void condet (uniform int threads, uniform partitioning tree[], uniform master_conpnt master[],
  uniform int parnum, uniform int ellnum, uniform int ellcol[], uniform int part[], uniform REAL * uniform center[6],
  uniform REAL * uniform radii[3], uniform REAL * uniform orient[18], uniform int trinum,
  uniform int tricol[], uniform int triobs[], uniform REAL * uniform tri[3][3])
{
  if (tree == NULL) return;

  launch [threads] update_existing (parnum/threads, parnum, master, center, radii, orient, tri);
  sync;

  launch [threads] test_ellipsoids (ellnum/threads, tree, ellnum, ellcol, part, center, radii, orient, master);

  launch [threads] test_triangles (trinum/threads, tree, trinum, tricol, triobs, tri, master);
}

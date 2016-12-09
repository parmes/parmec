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

#include <Python.h>
#include <structmember.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include "macros.h"
#include "parmec.h"
#include "mem.h"
#include "map.h"

namespace parmec
{
int output_frame = 0; /* output files frame */
}

using namespace parmec;
using namespace std;

/* output vtk dataset of rigid body data */
static void output_rb_dataset (int num, int *set, int ent, ofstream &out)
{
  int i, j;

  out << "DATASET UNSTRUCTURED_GRID\n";

  out << "POINTS " << num << " float\n";
  for (i = 0; i < num; i ++)
  {
    j = set[i];
    out << position [0][j] << " " << position[1][j] << " " << position[2][j] << "\n";
  }

  if (ent & (OUT_NUMBER|OUT_DISPL|OUT_ORIENT|OUT_LINVEL|OUT_ANGVEL|OUT_FORCE|OUT_TORQUE))
  {
    out << "POINT_DATA " << num << "\n";
  }

  if (ent & OUT_NUMBER)
  {
    out << "SCALARS number int\n";
    out << "LOOKUP_TABLE default\n";
    for (i = 0; i < num; i ++)
    {
      out << set[i] << "\n";
    }
  }

  if (ent & OUT_DISPL)
  {
    out << "VECTORS displ float\n";

    for (i = 0; i < num; i ++)
    {
      j = set[i];

      REAL d[3] = {position[0][j]-position[3][j],
                   position[1][j]-position[4][j],
		   position[2][j]-position[5][j]};

      out << d[0] << " " << d[1] << " " << d[2] << "\n";
    }
  }

  if (ent & OUT_ORIENT)
  {
    out << "VECTORS orient1 float\n";
    for (i = 0; i < num; i ++)
    {
      j = set[i];
      out << rotation[0][j] << " " << rotation[1][j] << " " << rotation[2][j] << "\n";
    }

    out << "VECTORS orient2 float\n";
    for (i = 0; i < num; i ++)
    {
      j = set[i];
      out << rotation[3][j] << " " << rotation[4][j] << " " << rotation[5][j] << "\n";
    }

    out << "VECTORS orient3 float\n";
    for (i = 0; i < num; i ++)
    {
      j = set[i];
      out << rotation[6][j] << " " << rotation[7][j] << " " << rotation[8][j] << "\n";
    }
  }

  if (ent & OUT_LINVEL)
  {
    out << "VECTORS linvel float\n";

    for (i = 0; i < num; i ++)
    {
      j = set[i];

      out << linear[0][j] << " " << linear[1][j] << " " << linear[2][j] << "\n";
    }
  }

  if (ent & OUT_ANGVEL)
  {
    out << "VECTORS angvel float\n";

    for (i = 0; i < num; i ++)
    {
      j = set[i];

      out << angular[0][j] << " " << angular[1][j] << " " << angular[2][j] << "\n";
    }
  }

  if (ent & OUT_FORCE)
  {
    out << "VECTORS force float\n";

    for (i = 0; i < num; i ++)
    {
      j = set[i];

      out << force[0][j] << " " << force[1][j] << " " << force[2][j] << "\n";
    }
  }

  if (ent & OUT_TORQUE)
  {
    out << "VECTORS torque float\n";

    for (i = 0; i < num; i ++)
    {
      j = set[i];

      out << torque[0][j] << " " << torque[1][j] << " " << torque[2][j] << "\n";
    }
  }
}

/* output vtk dataset of triangles */
static void output_triangle_dataset (int num, int *set, int ent, ofstream &out)
{
  int i, j;

  out << "DATASET UNSTRUCTURED_GRID\n";

  out << "POINTS " << 3*num << " float\n";
  for (i = 0; i < num; i ++)
  {
    j = set[i];
    out << tri [0][0][j] << " " << tri[0][1][j] << " " << tri[0][2][j] << " "
	<< tri [1][0][j] << " " << tri[1][1][j] << " " << tri[1][2][j] << " "
	<< tri [2][0][j] << " " << tri[2][1][j] << " " << tri[2][2][j] << "\n";
  }

  out << "CELLS " << num << " " << 4*num << "\n";
  for (i = 0; i < num; i ++)
  {
    out << 3 << " " << 3*i << " " << 3*i+1 << " " << 3*i+2 << "\n";
  }

  out << "CELL_TYPES " << num << "\n";
  for (i = 0; i < num; i ++)
  {
    out << 5 << "\n";
  }

  if (ent & (OUT_DISPL|OUT_LINVEL))
  {
    out << "POINT_DATA " << 3*num << "\n";
  }

  if (ent & OUT_DISPL)
  {
    out << "VECTORS displ float\n";

    for (i = 0; i < num; i ++)
    {
      j = set[i];

      REAL p1[3] = {tri [0][0][j], tri[0][1][j], tri[0][2][j]};
      REAL p2[3] = {tri [1][0][j], tri[1][1][j], tri[1][2][j]};
      REAL p3[3] = {tri [2][0][j], tri[2][1][j], tri[2][2][j]};

      j = triobs[set[i]];

      if (j >= 0) /* particle */
      {
	REAL L[9], X[3], x[3], P[3], d[3];

	L[0] = rotation[0][j];
	L[1] = rotation[1][j];
	L[2] = rotation[2][j];
	L[3] = rotation[3][j];
	L[4] = rotation[4][j];
	L[5] = rotation[5][j];
	L[6] = rotation[6][j];
	L[7] = rotation[7][j];
	L[8] = rotation[8][j];

	x[0] = position[0][j];
	x[1] = position[1][j];
	x[2] = position[2][j];

	X[0] = position[3][j];
	X[1] = position[4][j];
	X[2] = position[5][j];

	SUB (p1, x, d);
	TVADDMUL (X, L, d, P);
	SUB (p1, P, d);
        out << d[0] << " " << d[1] << " " << d[2] << "\n";

	SUB (p2, x, d);
	TVADDMUL (X, L, d, P);
	SUB (p2, P, d);
        out << d[0] << " " << d[1] << " " << d[2] << "\n";

	SUB (p3, x, d);
	TVADDMUL (X, L, d, P);
	SUB (p3, P, d);
        out << d[0] << " " << d[1] << " " << d[2] << "\n";
      }
      else if (j == -1) /* static obstacle */
      {
	out << "0.0 0.0 0.0\n";
	out << "0.0 0.0 0.0\n";
	out << "0.0 0.0 0.0\n";
      }
      else /* moving obstacle */
      {
	j = -j-2;

	REAL x[3] = {obspnt[3*j], obspnt[3*j+1], obspnt[3*j+2]}, d[3];

	SUB (p1, x, d);
        out << d[0] << " " << d[1] << " " << d[2] << "\n";
	SUB (p2, x, d);
        out << d[0] << " " << d[1] << " " << d[2] << "\n";
	SUB (p3, x, d);
        out << d[0] << " " << d[1] << " " << d[2] << "\n";
      }
    }
  }

  if (ent & OUT_LINVEL)
  {
    out << "VECTORS linvel float\n";

    for (i = 0; i < num; i ++)
    {
      j = set[i];

      REAL p1[3] = {tri [0][0][j], tri[0][1][j], tri[0][2][j]};
      REAL p2[3] = {tri [1][0][j], tri[1][1][j], tri[1][2][j]};
      REAL p3[3] = {tri [2][0][j], tri[2][1][j], tri[2][2][j]};

      j = triobs[set[i]];

      if (j >= 0) /* particle */
      {
	REAL x[3], a[3], o[3], v[3], w[3];

	x[0] = position[0][j];
	x[1] = position[1][j];
	x[2] = position[2][j];

	o[0] = angular[3][j];
	o[1] = angular[4][j];
	o[2] = angular[5][j];

	v[0] = linear[0][j];
	v[1] = linear[1][j];
	v[2] = linear[2][j];

        COPY (v, w);
	SUB (p1, x, a);
	PRODUCTADD (o, a, w);
        out << w[0] << " " << w[1] << " " << w[2] << "\n";

        COPY (v, w);
	SUB (p2, x, a);
	PRODUCTADD (o, a, w);
        out << w[0] << " " << w[1] << " " << w[2] << "\n";

        COPY (v, w);
	SUB (p3, x, a);
	PRODUCTADD (o, a, w);
        out << w[0] << " " << w[1] << " " << w[2] << "\n";
      }
      else if (j == -1) /* static obstacle */
      {
	out << "0.0 0.0 0.0\n";
	out << "0.0 0.0 0.0\n";
	out << "0.0 0.0 0.0\n";
      }
      else /* moving obstacle */
      {
	j = -j-2;

	REAL x[3] = {obspnt[3*j], obspnt[3*j+1], obspnt[3*j+2]};
	REAL o[3] = {obsang[3*j], obsang[3*j+1], obsang[3*j+2]};
	REAL v[3] = {obslin[3*j], obslin[3*j+1], obslin[3*j+2]};
	REAL a[3], w[3];

        COPY (v, w);
	SUB (p1, x, a);
	PRODUCTADD (o, a, w);
        out << w[0] << " " << w[1] << " " << w[2] << "\n";

        COPY (v, w);
	SUB (p2, x, a);
	PRODUCTADD (o, a, w);
        out << w[0] << " " << w[1] << " " << w[2] << "\n";

        COPY (v, w);
	SUB (p3, x, a);
	PRODUCTADD (o, a, w);
        out << w[0] << " " << w[1] << " " << w[2] << "\n";
      }
    }
  }

  if (ent & (OUT_NUMBER|OUT_COLOR|OUT_ANGVEL|OUT_FORCE|OUT_TORQUE))
  {
    out << "CELL_DATA " << num << "\n";
  }

  if (ent & OUT_NUMBER)
  {
    out << "SCALARS number int\n";
    out << "LOOKUP_TABLE default\n";
    for (i = 0; i < num; i ++)
    {
      out << triobs[set[i]] << "\n";
    }
  }


  if (ent & OUT_COLOR)
  {
    out << "SCALARS colors int\n";
    out << "LOOKUP_TABLE default\n";
    for (i = 0; i < num; i ++)
    {
      out << tricol[set[i]] << "\n";
    }
  }

  if (ent & OUT_ANGVEL)
  {
    out << "VECTORS angvel float\n";
    for (i = 0; i < num; i ++)
    {
      j = triobs[set[i]];

      if (j >= 0) /* particle */
      {
        out << angular[3][j] << " " << angular[4][j] << " " << angular[5][j] << "\n";
      }
      else if (j == -1) /* static obstacle */
      {
	out << "0.0 0.0 0.0\n";
      }
      else /* moving obstacle */
      {
	j = -j-2;

        out << obsang[3*j] << " " << obsang[3*j+1] << " " << obsang[3*j+2] << "\n";
      }
    }
  }

  if (ent & OUT_FORCE)
  {
    out << "VECTORS force float\n";
    for (i = 0; i < num; i ++)
    {
      j = triobs[set[i]];

      if (j >= 0)
      {
        out << force[0][j] << " " << force[1][j] << " " << force[2][j] << "\n";
      }
      else
      {
	out << "0.0 0.0 0.0\n";
      }
    }
  }

  if (ent & OUT_TORQUE)
  {
    out << "VECTORS torque float\n";
    for (i = 0; i < num; i ++)
    {
      j = triobs[set[i]];

      if (j >= 0)
      {
        out << torque[0][j] << " " << torque[1][j] << " " << torque[2][j] << "\n";
      }
      else
      {
	out << "0.0 0.0 0.0\n";
      }
    }
  }
}

/* find triangles belonging to a particle set [part0, part1) */
static int find_triangle_set (int *part0, int *part1, int *triangles)
{
  int num = 0;

  for (int *p = part0; p < part1; p ++)
  {
    int s = 0, e = trinum, i, j;

    while (s < e)
    {
      i = (s+e)/2;

      for (j = triobs[i]; i < e && j < 0; i++) j = triobs[i]; /* skip obstacles */

      if (j == *p) break;
      else if (j < *p) s = i;
      else e = i;
    }

    if (j == *p)
    {
      while (i > 0 && triobs[i-1] == j) i --; /* find start of range */

      while (i < trinum && triobs[i] == j)
      {
	triangles[num] = i;
	num ++;
	i ++;
      }
    }
  }

  return num;
}

/* output vtk dataset of spring data */
static void output_spring_dataset (int num, int *set, int ent, ofstream &out)
{
  int i, j;

  out << "DATASET UNSTRUCTURED_GRID\n";

  out << "POINTS " << num << " float\n";
  for (i = 0; i < num; i ++)
  {
    j = set[i];
    out << 0.5*(sprpnt[0][0][j]+sprpnt[1][0][j]) << " "
        << 0.5*(sprpnt[0][1][j]+sprpnt[1][1][j]) << " "
	<< 0.5*(sprpnt[0][2][j]+sprpnt[1][2][j]) << "\n";
  }

  if (ent & (OUT_NUMBER|OUT_DISPL|OUT_F|OUT_SF|OUT_PAIR))
  {
    out << "POINT_DATA " << num << "\n";
  }

  if (ent & OUT_NUMBER)
  {
    out << "SCALARS number int\n";
    out << "LOOKUP_TABLE default\n";
    for (i = 0; i < num; i ++)
    {
      j = set[i];
      out << sprid[j] << "\n";
    }
  }

  if (ent & OUT_DISPL)
  {
    out << "SCALARS displ float\n";
    out << "LOOKUP_TABLE default\n";
    for (i = 0; i < num; i ++)
    {
      j = set[i];
      out << stroke[0][j] << "\n";
    }
  }

  if (ent & OUT_ORIENT)
  {
    out << "VECTORS orient float\n";
    for (i = 0; i < num; i ++)
    {
      j = set[i];
      out << (sprpnt[1][0][j]-sprpnt[0][0][j]) << " "
	  << (sprpnt[1][1][j]-sprpnt[0][1][j]) << " "
	  << (sprpnt[1][2][j]-sprpnt[0][2][j]) << "\n";
    }
  }

  if (ent & OUT_F)
  {
    out << "SCALARS F float\n";
    out << "LOOKUP_TABLE default\n";
    for (i = 0; i < num; i ++)
    {
      j = set[i];
      out << sprfrc[0][j] << "\n";
    }
  }

  if (ent & OUT_SF)
  {
    out << "SCALARS SF float\n";
    out << "LOOKUP_TABLE default\n";
    for (i = 0; i < num; i ++)
    {
      j = set[i];
      out << sprfrc[1][j] << "\n";
    }
  }
}

/* output files */
void output_files ()
{
  ostringstream oss;
  ofstream out;
  int i, j;

  if (ellnum)
  {
    oss.str("");
    oss.clear();
    oss << outpath << ".dump";
    if (output_frame) out.open (oss.str().c_str(), ios::app);
    else out.open (oss.str().c_str());

    out << "ITEM: TIMESTEP\n";
    out << curtime << "\n";
    out << "ITEM: NUMBER OF ATOMS\n";
    out << ellnum << "\n";
    out << "ITEM: BOX BOUNDS\n";
    out << "0 1\n";
    out << "0 1\n";
    out << "0 1\n";
    out << "ITEM: ATOMS id x y z radius\n";
    for (i = 0; i < ellnum; i ++)
    {
      out << i+1 << " " << center[0][i] << " " << center[1][i] << " " << center[2][i] << " " <<  radii[0][i] << "\n";
    }

    out.close();

    /* TODO --> include output entities and subsets */
  }

  if (parnum)
  {
    int num, *set; /* number of and the of unselected particles */

    ERRMEM (set = new int[parnum]);

    for (j = -1; j < outnum; j ++)
    {
      if (j < 0 && (outrest[1] & OUT_MODE_RB)) /* output unselected particles */
      {
	oss.str("");
	oss.clear();
	oss << outpath << j+1 << "rb.vtk." << output_frame;
	out.open (oss.str().c_str());

	out << "# vtk DataFile Version 2.0\n";
	out << "PARMEC rigid bodies output at time " << curtime << "\n";
	out << "ASCII\n";

	for (num = i = 0; i < parnum; i ++)
	{
	  if (flags[i] & OUTREST)
	  {
	    set[num ++] = i;
	  }
	}

	output_rb_dataset (num, set, outrest[0], out);

        out.close();
      }
      else if (outmode[j] & OUT_MODE_RB) /* output selected particles */
      {
	oss.str("");
	oss.clear();
	oss << outpath << j+1 << "rb.vtk." << output_frame;
	out.open (oss.str().c_str());

	out << "# vtk DataFile Version 2.0\n";
	out << "PARMEC rigid bodies output at time " << curtime << "\n";
	out << "ASCII\n";

	output_rb_dataset (outidx[j+1]-outidx[j], &outpart[outidx[j]], outent[j], out);

        out.close();
      }
    }

    delete set;
  }

  if (trinum)
  {
    int num, *set; /* number of and the set of triangles */

    ERRMEM (set = new int[trinum]);

    for (j = -1; j < outnum; j ++)
    {
      if (j < 0 && (outrest[1] & OUT_MODE_MESH)) /* output unselected triangles */
      {
	oss.str("");
	oss.clear();
	oss << outpath << j+1 << ".vtk." << output_frame;
	out.open (oss.str().c_str());

	out << "# vtk DataFile Version 2.0\n";
	out << "PARMEC triangles output at time " << curtime << "\n";
	out << "ASCII\n";

	for (num = i = 0; i < trinum; i ++)
	{
	  if (triobs[i] >= 0 && (flags[triobs[i]] & OUTREST)) /* triangles of unselected particles */
	  {
	    set[num ++] = i;
	  }
	  else if (triobs[i] < 0) /* triangles of obstacles */
	  {
	    set[num ++] = i;
	  }
	}

	output_triangle_dataset (num, set, outrest[0], out);

	out.close();
      }
      else if (outmode[j] & OUT_MODE_MESH) /* output selected triangles */
      {
	oss.str("");
	oss.clear();
	oss << outpath << j+1 << ".vtk." << output_frame;
	out.open (oss.str().c_str());

	out << "# vtk DataFile Version 2.0\n";
	out << "PARMEC triangles output at time " << curtime << "\n";
	out << "ASCII\n";

	num = find_triangle_set (&outpart[outidx[j]], &outpart[outidx[j+1]], set);

	if (num) output_triangle_dataset (num, set, outent[j], out);

        out.close();
      }
    }

    delete set;
  }

  /* TODO --> *cd.vtk.* contact data output */

  if (sprnum)
  {
    int num, *set; /* number of and the set of springs */
    MAP **map, *item; /* map of springs attached to particles, and iterator */
    MEM mem;

    MEM_Init (&mem, sizeof (MAP), 1024);
    ERRMEM (map = static_cast<MAP**>(MEM_CALLOC(parnum * sizeof(MAP*))));
    ERRMEM (set = new int[sprnum]);

    for (j = -1; j < outnum; j ++)
    {
      if (j < 0 && (outrest[1] & OUT_MODE_SD)) /* output springs attached to unselected particles */
      {
	oss.str("");
	oss.clear();
	oss << outpath << j+1 << "sd.vtk." << output_frame;
	out.open (oss.str().c_str());

	out << "# vtk DataFile Version 2.0\n";
	out << "PARMEC springs output at time " << curtime << "\n";
	out << "ASCII\n";

	for (num = i = 0; i < sprnum; i ++)
	{
	  if (flags[sprpart[0][i]] && OUTREST || /* first or second particle is unselected */
	  (sprpart[1][i] >= 0 && (flags[sprpart[1][i]] & OUTREST)))
	  {
	    set[num ++] = i;
	  }

	  MAP_Insert (&mem, &map[sprpart[0][i]], (void*)(long)i, NULL, NULL); /* map springs to particles */
	  if (sprpart[1][i] >= 0) MAP_Insert (&mem, &map[sprpart[1][i]], (void*)(long)i, NULL, NULL);
	}

	output_spring_dataset (num, set, outrest[0], out);

        out.close();
      }
      else if (outmode[j] & OUT_MODE_SD) /* output springs attached to selected particles */
      {
	oss.str("");
	oss.clear();
	oss << outpath << j+1 << "sd.vtk." << output_frame;
	out.open (oss.str().c_str());

	out << "# vtk DataFile Version 2.0\n";
	out << "PARMEC springs output at time " << curtime << "\n";
	out << "ASCII\n";

	for (num = 0, i = outidx[j]; i < outidx[j+1]; i ++)
	{
	  for (item = MAP_First (map[outpart[i]]); item; item = MAP_Next(item))
	  {
	    set[num ++] = (int)(long)item->key;
	  }
	}

	if (num) output_spring_dataset (num, set, outent[j], out);

        out.close();
      }
    }

    delete set;
    free (map);
    MEM_Release (&mem);
  }

  output_frame ++; 
}

/* output time history */
void output_history ()
{
  int i, j, k;

  for (i = 0; i < hisnum; i ++) /* append time histories */
  {
    if (hisent[i] == HIS_TIME)
    {
      PyList_Append ((PyObject*)history[i], PyFloat_FromDouble(curtime));
    }
    else switch (hiskind[i]&(HIS_LIST|HIS_SPHERE|HIS_BOX))
    {
    case HIS_LIST:
    {
      if (hiskind[i] & HIS_POINT) /* one particle point based */
      {
	k = hislst[hisidx[i]];

	REAL x[3] = {position[0][k], position[1][k], position[2][k]};
	REAL X[3] = {position[3][k], position[4][k], position[5][k]};
	REAL v[3] = {linear[0][k], linear[1][k], linear[2][k]};
	REAL o[3] = {angular[3][k], angular[4][k], angular[5][k]};
	REAL L[9] = {rotation[0][k], rotation[1][k], rotation[2][k],
	             rotation[3][k], rotation[4][k], rotation[5][k],
		     rotation[6][k], rotation[7][k], rotation[8][k]};
        REAL P[3] = {source[0][i], source[1][i], source[2][i]};
	REAL Q[3], p[3], a[3], value;

	SUB (P, X, Q);
	NVADDMUL (x, L, Q, p);
	SUB (p, x, a);

	switch (hisent[i])
	{
	case HIS_PX:
	  value = p[0];
	break;
	case HIS_PY:
	  value = p[1];
	break;
	case HIS_PZ:
	  value = p[2];
	break;
	case HIS_PL:
	  value = LEN(p);
	break;
	case HIS_DX:
	  value = p[0]-P[0];
	break;
	case HIS_DY:
	  value = p[1]-P[1];
	break;
	case HIS_DZ:
	  value = p[2]-P[2];
	break;
	case HIS_DL:
	{
	  REAL q[3] = {p[0]-P[0], p[1]-P[1], p[2]-P[2]};
	  value = LEN(q);
	}
	break;
	case HIS_VX:
	  value = v[0] + a[1]*o[2] - a[2]*o[1];
	break;
	case HIS_VY:
	  value = v[1] + a[2]*o[0] - a[0]*o[2];
	break;
	case HIS_VZ:
	  value = v[2] + a[0]*o[1] - a[1]*o[0];
	break;
	case HIS_VL:
	{
	  REAL q[3] = {v[0] + a[1]*o[2] - a[2]*o[1], v[1] + a[2]*o[0] - a[0]*o[2], v[2] + a[0]*o[1] - a[1]*o[0]};
	  value = LEN(q);
	}
	break;
	case HIS_OX:
	  value = o[0];
	break;
	case HIS_OY:
	  value = o[1];
	break;
	case HIS_OZ:
	  value = o[2];
	break;
	case HIS_OL:
	{
	  value = LEN(o);
	}
	break;
	case HIS_FX:
	  value = force[0][k];
	break;
	case HIS_FY:
	  value = force[1][k];
	break;
	case HIS_FZ:
	  value = force[2][k];
	break;
	case HIS_FL:
	{
	  REAL q[3] = {force[0][k], force[1][k], force[2][k]};
	  value = LEN(q);
	}
	break;
	case HIS_TX:
	  value = torque[0][k];
	break;
	case HIS_TY:
	  value = torque[1][k];
	break;
	case HIS_TZ:
	  value = torque[2][k];
	break;
	case HIS_TL:
	{
	  REAL q[3] = {torque[0][k], torque[1][k], torque[2][k]};
	  value = LEN(q);
	}
	break;
	}

        PyList_Append ((PyObject*)history[i], PyFloat_FromDouble(value));
      }
      else /* particle or spring list based */
      {
	REAL value = 0.0;

	for (j = hisidx[i]; j < hisidx[i+1]; j ++)
	{
	  k = hislst[j];

	  switch (hisent[i])
	  {
	  case HIS_PX:
	    value += position[0][k];
	  break;
	  case HIS_PY:
	    value += position[1][k];
	  break;
	  case HIS_PZ:
	    value += position[2][k];
	  break;
	  case HIS_PL:
	  {
	    REAL q[3] = {position[0][k], position[1][k], position[2][k]};
	    value += LEN(q);
	  }
	  break;
	  case HIS_DX:
	    value += position[0][k]-position[3][k];
	  break;
	  case HIS_DY:
	    value += position[1][k]-position[4][k];
	  break;
	  case HIS_DZ:
	    value += position[2][k]-position[5][k];
	  break;
	  case HIS_DL:
	  {
	    REAL q[3] = {position[0][k]-position[3][k], position[1][k]-position[4][k], position[2][k]-position[5][k]};
	    value += LEN(q);
	  }
	  break;
	  case HIS_VX:
	    value += linear[0][k];
	  break;
	  case HIS_VY:
	    value += linear[1][k];
	  break;
	  case HIS_VZ:
	    value += linear[2][k];
	  break;
	  case HIS_VL:
	  {
	    REAL q[3] = {linear[0][k], linear[1][k], linear[2][k]};
	    value += LEN(q);
	  }
	  break;
	  case HIS_OX:
	    value += angular[3][k];
	  break;
	  case HIS_OY:
	    value += angular[4][k];
	  break;
	  case HIS_OZ:
	    value += angular[5][k];
	  break;
	  case HIS_OL:
	  {
	    REAL q[3] = {angular[3][k], angular[4][k], angular[5][k]};
	    value += LEN(q);
	  }
	  break;
          case HIS_FX:
	    value += force[0][k];
	  break;
	  case HIS_FY:
	    value += force[1][k];
	  break;
	  case HIS_FZ:
	    value += force[2][k];
	  break;
	  case HIS_FL:
	  {
	    REAL q[3] = {force[0][k], force[1][k], force[2][k]};
	    value += LEN(q);
	  }
	  break;
          case HIS_TX:
	    value += torque[0][k];
	  break;
	  case HIS_TY:
	    value += torque[1][k];
	  break;
	  case HIS_TZ:
	    value += torque[2][k];
	  break;
	  case HIS_TL:
	  {
	    REAL q[3] = {torque[0][k], torque[1][k], torque[2][k]};
	    value += LEN(q);
	  }
	  break;
	  case HIS_STROKE:
	  {
	    int l = sprmap[k];
	    value += stroke[0][l];
	  }
	  break;
	  case HIS_STF:
	  {
	    int l = sprmap[k];
	    value += sprfrc[0][l];
	  }
	  break;
	  case HIS_SF:
	  {
	    int l = sprmap[k];
	    value += sprfrc[1][l];
	  }
	  break;
	  }
	}

	j = hisidx[i+1]-hisidx[i];

        PyList_Append ((PyObject*)history[i], PyFloat_FromDouble(value/(REAL)j));
      }
    }
    break;
    case HIS_SPHERE:
    {
      ASSERT (0, "Sphere based time history is not yet implemented");
    }
    break;
    case HIS_BOX:
    {
      ASSERT (0, "Box based time history is not yet implemented");
    }
    break;
    }
  }
}

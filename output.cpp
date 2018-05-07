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
#include <hdf5.h>
#include <hdf5_hl.h>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include "macros.h"
#include "parmec.h"
#include "mem.h"
#include "map.h"

#if MED
extern "C" {
#include "med.h"
}
static med_idt med_fid = (med_idt)-1;
#endif

using namespace parmec;
using namespace std;

/* output dump dataset of spheres and ellipsoids */
static void output_dump_dataset (int num, int *set, int ent, ofstream &out)
{
  int i, j, k;

  out << "ITEM: TIMESTEP\n";
  out << curtime << "\n";
  out << "ITEM: NUMBER OF ATOMS\n";
  out << num << "\n";
  out << "ITEM: BOX BOUNDS\n";
  out << "0 1\n";
  out << "0 1\n";
  out << "0 1\n";
  out << "ITEM: ATOMS id x y z radius";
  if (ent & OUT_LINVEL) out << " vx vy vz";
  if (ent & OUT_ANGVEL) out << " omegax omegay omegaz";
  if (ent & OUT_FORCE) out << " fx fy fz";
  if (ent & OUT_TORQUE) out << " tqx tqy tqz";
  out << "\n";

  for (i = 0; i < num; i ++)
  {
    j = set[i];

    k = part[j];

    out << j+1 << " " << center[0][j] << " " << center[1][j] << " " << center[2][j] << " " <<  radii[0][j];

    if (ent & OUT_LINVEL) out << " " << linear[0][k] << " " << linear[1][k] << " " << linear[2][k];

    if (ent & OUT_ANGVEL) out << " " << angular[0][k] << " " << angular[1][k] << " " << angular[2][k];

    if (ent & OUT_FORCE) out << " " << force[0][k] << " " << force[1][k] << " " << force[2][k];

    if (ent & OUT_TORQUE) out << " " << torque[0][k] << " " << torque[1][k] << " " << torque[2][k];

    out << "\n";
  }

  /* NUMBER is already handled by the 'id' */

  /* TODO: handle the below by means of d_name fields */

#if 0
  if (ent & OUT_DISPL)
  {
    out << "VECTORS DISPL float\n";

    for (i = 0; i < num; i ++)
    {
      j = set[i];

      REAL d[3] = {position[0][j]-position[3][j],
                   position[1][j]-position[4][j],
		   position[2][j]-position[5][j]};

      out << d[0] << " " << d[1] << " " << d[2] << "\n";
    }
  }

  if (ent & OUT_ORIENT1)
  {
    out << "VECTORS ORIENT1 float\n";
    for (i = 0; i < num; i ++)
    {
      j = set[i];
      out << rotation[0][j] << " " << rotation[1][j] << " " << rotation[2][j] << "\n";
    }
  }

  if (ent & OUT_ORIENT2)
  {
    out << "VECTORS ORIENT2 float\n";
    for (i = 0; i < num; i ++)
    {
      j = set[i];
      out << rotation[3][j] << " " << rotation[4][j] << " " << rotation[5][j] << "\n";
    }
  }

  if (ent & OUT_ORIENT3)
  {
    out << "VECTORS ORIENT3 float\n";
    for (i = 0; i < num; i ++)
    {
      j = set[i];
      out << rotation[6][j] << " " << rotation[7][j] << " " << rotation[8][j] << "\n";
    }
  }
#endif
}

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

  if (ent & (OUT_NUMBER|OUT_DISPL|OUT_ORIENT1|OUT_ORIENT2|OUT_ORIENT3|OUT_LINVEL|OUT_ANGVEL|OUT_FORCE|OUT_TORQUE))
  {
    out << "POINT_DATA " << num << "\n";
  }

  if (ent & OUT_NUMBER)
  {
    out << "SCALARS NUMBER int\n";
    out << "LOOKUP_TABLE default\n";
    for (i = 0; i < num; i ++)
    {
      out << set[i] << "\n";
    }
  }

  if (ent & OUT_DISPL)
  {
    out << "VECTORS DISPL float\n";

    for (i = 0; i < num; i ++)
    {
      j = set[i];

      REAL d[3] = {position[0][j]-position[3][j],
                   position[1][j]-position[4][j],
		   position[2][j]-position[5][j]};

      out << d[0] << " " << d[1] << " " << d[2] << "\n";
    }
  }

  if (ent & OUT_ORIENT1)
  {
    out << "VECTORS ORIENT1 float\n";
    for (i = 0; i < num; i ++)
    {
      j = set[i];
      out << rotation[0][j] << " " << rotation[1][j] << " " << rotation[2][j] << "\n";
    }
  }

  if (ent & OUT_ORIENT2)
  {
    out << "VECTORS ORIENT2 float\n";
    for (i = 0; i < num; i ++)
    {
      j = set[i];
      out << rotation[3][j] << " " << rotation[4][j] << " " << rotation[5][j] << "\n";
    }
  }

  if (ent & OUT_ORIENT3)
  {
    out << "VECTORS ORIENT3 float\n";
    for (i = 0; i < num; i ++)
    {
      j = set[i];
      out << rotation[6][j] << " " << rotation[7][j] << " " << rotation[8][j] << "\n";
    }
  }

  if (ent & OUT_LINVEL)
  {
    out << "VECTORS LINVEL float\n";

    for (i = 0; i < num; i ++)
    {
      j = set[i];

      out << linear[0][j] << " " << linear[1][j] << " " << linear[2][j] << "\n";
    }
  }

  if (ent & OUT_ANGVEL)
  {
    out << "VECTORS ANGVEL float\n";

    for (i = 0; i < num; i ++)
    {
      j = set[i];

      out << angular[0][j] << " " << angular[1][j] << " " << angular[2][j] << "\n";
    }
  }

  if (ent & OUT_FORCE)
  {
    out << "VECTORS FORCE float\n";

    for (i = 0; i < num; i ++)
    {
      j = set[i];

      out << force[0][j] << " " << force[1][j] << " " << force[2][j] << "\n";
    }
  }

  if (ent & OUT_TORQUE)
  {
    out << "VECTORS TORQUE float\n";

    for (i = 0; i < num; i ++)
    {
      j = set[i];

      out << torque[0][j] << " " << torque[1][j] << " " << torque[2][j] << "\n";
    }
  }
}

/* output hdf5 dataset of rigid body data */
static void h5_rb_dataset (int num, int *set, int ent, hid_t h5_step)
{
  double *data, *pdata;
  int i, j, *numb;

  ERRMEM (data = new double [9*num]);
  for (i = 0; i < num; i ++)
  {
    j = set[i];
    data[3*i+0] = position[0][j];
    data[3*i+1] = position[1][j];
    data[3*i+2] = position[2][j];
  }
  hsize_t dims[2] = {num, 3};
  ASSERT (H5LTmake_dataset_double (h5_step, "GEOM", 2, dims, data) >= 0, "HDF5 file write error");

  if (ent & OUT_NUMBER)
  {
    ERRMEM (numb = new int[num]);

    for (i = 0; i < num; i ++)
    {
      numb[i] = set[i];
    }

    hsize_t length = num;
    ASSERT (H5LTmake_dataset_int (h5_step, "NUMBER", 1, &length, numb) >= 0, "HDF5 file write error");

    delete [] numb;
  }

  if (ent & OUT_DISPL)
  {
    pdata = data;
    for (i = 0; i < num; i ++)
    {
      j = set[i];

      REAL d[3] = {position[0][j]-position[3][j],
                   position[1][j]-position[4][j],
		   position[2][j]-position[5][j]};

      pdata[0] = d[0]; pdata[1] = d[1]; pdata[2] = d[2]; pdata +=3;
    }

    ASSERT (H5LTmake_dataset_double (h5_step, "DISPL", 2, dims, data) >= 0, "HDF5 file write error");
  }

  if (ent & OUT_LINVEL)
  {
    pdata = data;
    for (i = 0; i < num; i ++)
    {
      j = set[i];

      pdata[0] = linear[0][j]; pdata[1] = linear[1][j]; pdata[2] = linear[2][j]; pdata += 3;
    }

    ASSERT (H5LTmake_dataset_double (h5_step, "LINVEL", 2, dims, data) >= 0, "HDF5 file write error");
  }

  if (ent & OUT_ANGVEL)
  {
    pdata = data;
    for (i = 0; i < num; i ++)
    {
      j = set[i];

      pdata[0] = angular[0][j]; pdata[1] = angular[1][j]; pdata[2] = angular[2][j]; pdata += 3;
    }

    ASSERT (H5LTmake_dataset_double (h5_step, "ANGVEL", 2, dims, data) >= 0, "HDF5 file write error");
  }

  if (ent & OUT_FORCE)
  {
    pdata = data;
    for (i = 0; i < num; i ++)
    {
      j = set[i];

      pdata[0] = force[0][j]; pdata[1] = force[1][j]; pdata[2] = force[2][j]; pdata += 3;
    }

    ASSERT (H5LTmake_dataset_double (h5_step, "FORCE", 2, dims, data) >= 0, "HDF5 file write error");
  }

  if (ent & OUT_TORQUE)
  {
    pdata = data;
    for (i = 0; i < num; i ++)
    {
      j = set[i];

      pdata[0] = torque[0][j]; pdata[1] = torque[1][j]; pdata[2] = torque[2][j]; pdata += 3;
    }

    ASSERT (H5LTmake_dataset_double (h5_step, "TORQUE", 2, dims, data) >= 0, "HDF5 file write error");
  }

  if (ent & OUT_ORIENT)
  {
    pdata = data;
    for (i = 0; i < num; i ++)
    {
      j = set[i];
      pdata[0] = rotation[0][j]; pdata[1] = rotation[1][j]; pdata[2] = rotation[2][j]; pdata += 3;
      pdata[0] = rotation[3][j]; pdata[1] = rotation[4][j]; pdata[2] = rotation[5][j]; pdata += 3;
      pdata[0] = rotation[6][j]; pdata[1] = rotation[7][j]; pdata[2] = rotation[8][j]; pdata += 3;
    }

    hsize_t dims[2] = {num, 9};
    ASSERT (H5LTmake_dataset_double (h5_step, "ORIENT", 2, dims, data) >= 0, "HDF5 file write error");
  }

  if (ent & OUT_ORIENT1)
  {
    pdata = data;
    for (i = 0; i < num; i ++)
    {
      j = set[i];
      pdata[0] = rotation[0][j]; pdata[1] = rotation[1][j]; pdata[2] = rotation[2][j]; pdata += 3;
    }

    ASSERT (H5LTmake_dataset_double (h5_step, "ORIENT1", 2, dims, data) >= 0, "HDF5 file write error");
  }

  if (ent & OUT_ORIENT2)
  {
    pdata = data;
    for (i = 0; i < num; i ++)
    {
      j = set[i];
      pdata[0] = rotation[3][j]; pdata[1] = rotation[4][j]; pdata[2] = rotation[5][j]; pdata += 3;
    }

    ASSERT (H5LTmake_dataset_double (h5_step, "ORIENT2", 2, dims, data) >= 0, "HDF5 file write error");
  }

  if (ent & OUT_ORIENT3)
  {
    pdata = data;
    for (i = 0; i < num; i ++)
    {
      j = set[i];
      pdata[0] = rotation[6][j]; pdata[1] = rotation[7][j]; pdata[2] = rotation[8][j]; pdata += 3;
    }

    ASSERT (H5LTmake_dataset_double (h5_step, "ORIENT3", 2, dims, data) >= 0, "HDF5 file write error");
  }

  delete [] data;
}

/* output vtk dataset of triangles */
static void vtk_triangle_dataset (int num, int *set, int ent, ofstream &out)
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

/* output hdf5 dataset of triangles */
static void h5_triangle_dataset (int num, int *set, int ent, hid_t h5_step)
{
  double *data, *pdata;
  int i, j, *topo;

  ERRMEM (data = new double [9*num]);
  for (i = 0; i < num; i ++)
  {
    j = set[i];
    data[9*i+0] = tri[0][0][j];
    data[9*i+1] = tri[0][1][j];
    data[9*i+2] = tri[0][2][j];
    data[9*i+3] = tri[1][0][j];
    data[9*i+4] = tri[1][1][j];
    data[9*i+5] = tri[1][2][j];
    data[9*i+6] = tri[2][0][j];
    data[9*i+7] = tri[2][1][j]; 
    data[9*i+8] = tri[2][2][j];
  }
  hsize_t dims[2] = {3*num, 3};
  ASSERT (H5LTmake_dataset_double (h5_step, "GEOM", 2, dims, data) >= 0, "HDF5 file write error");

  ERRMEM (topo = new int[4*num]);
  for (i = 0; i < num; i ++)
  {
    topo[4*i+0] = 4;
    topo[4*i+1] = 3*i;
    topo[4*i+2] = 3*i+1;
    topo[4*i+3] = 3*i+2;
  }
  hsize_t length = 4*num;
  ASSERT (H5LTmake_dataset_int (h5_step, "TOPO", 1, &length, topo) >= 0, "HDF5 file write error");

  if (ent & OUT_DISPL)
  {
    pdata = data;
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
	pdata[0] = d[0]; pdata[1] = d[1]; pdata[2] = d[2]; pdata +=3;

	SUB (p2, x, d);
	TVADDMUL (X, L, d, P);
	SUB (p2, P, d);
	pdata[0] = d[0]; pdata[1] = d[1]; pdata[2] = d[2]; pdata +=3;

	SUB (p3, x, d);
	TVADDMUL (X, L, d, P);
	SUB (p3, P, d);
	pdata[0] = d[0]; pdata[1] = d[1]; pdata[2] = d[2]; pdata +=3;
      }
      else if (j == -1) /* static obstacle */
      {
	pdata[0] = 0.0; pdata[1] = 0.0; pdata[2] = 0.0; pdata +=3;
	pdata[0] = 0.0; pdata[1] = 0.0; pdata[2] = 0.0; pdata +=3;
	pdata[0] = 0.0; pdata[1] = 0.0; pdata[2] = 0.0; pdata +=3;
      }
      else /* moving obstacle */
      {
	j = -j-2;

	REAL x[3] = {obspnt[3*j], obspnt[3*j+1], obspnt[3*j+2]}, d[3];

	SUB (p1, x, d);
	pdata[0] = d[0]; pdata[1] = d[1]; pdata[2] = d[2]; pdata +=3;
	SUB (p2, x, d);
	pdata[0] = d[0]; pdata[1] = d[1]; pdata[2] = d[2]; pdata +=3;
	SUB (p3, x, d);
	pdata[0] = d[0]; pdata[1] = d[1]; pdata[2] = d[2]; pdata +=3;
      }
    }

    ASSERT (H5LTmake_dataset_double (h5_step, "DISPL", 2, dims, data) >= 0, "HDF5 file write error");
  }

  if (ent & OUT_LINVEL)
  {
    pdata = data;
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
	pdata[0] = w[0]; pdata[1] = w[1]; pdata[2] = w[2]; pdata +=3;

        COPY (v, w);
	SUB (p2, x, a);
	PRODUCTADD (o, a, w);
	pdata[0] = w[0]; pdata[1] = w[1]; pdata[2] = w[2]; pdata +=3;

        COPY (v, w);
	SUB (p3, x, a);
	PRODUCTADD (o, a, w);
	pdata[0] = w[0]; pdata[1] = w[1]; pdata[2] = w[2]; pdata +=3;
      }
      else if (j == -1) /* static obstacle */
      {
	pdata[0] = 0.0; pdata[1] = 0.0; pdata[2] = 0.0; pdata +=3;
	pdata[0] = 0.0; pdata[1] = 0.0; pdata[2] = 0.0; pdata +=3;
	pdata[0] = 0.0; pdata[1] = 0.0; pdata[2] = 0.0; pdata +=3;
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
	pdata[0] = w[0]; pdata[1] = w[1]; pdata[2] = w[2]; pdata +=3;

        COPY (v, w);
	SUB (p2, x, a);
	PRODUCTADD (o, a, w);
	pdata[0] = w[0]; pdata[1] = w[1]; pdata[2] = w[2]; pdata +=3;

        COPY (v, w);
	SUB (p3, x, a);
	PRODUCTADD (o, a, w);
	pdata[0] = w[0]; pdata[1] = w[1]; pdata[2] = w[2]; pdata +=3;
      }
    }

    ASSERT (H5LTmake_dataset_double (h5_step, "LINVEL", 2, dims, data) >= 0, "HDF5 file write error");
  }

  length = num;

  if (ent & OUT_NUMBER)
  {
    for (i = 0; i < num; i ++)
    {
      topo[i] = triobs[set[i]];
    }

    ASSERT (H5LTmake_dataset_int (h5_step, "NUMBER", 1, &length, topo) >= 0, "HDF5 file write error");
  }

  if (ent & OUT_COLOR)
  {
    for (i = 0; i < num; i ++)
    {
      topo[i] = tricol[set[i]];
    }

    ASSERT (H5LTmake_dataset_int (h5_step, "COLOR", 1, &length, topo) >= 0, "HDF5 file write error");
  }

  if (ent & OUT_ANGVEL)
  {
    pdata = data;
    for (i = 0; i < num; i ++)
    {
      j = triobs[set[i]];

      if (j >= 0) /* particle */
      {
	pdata[0] = angular[3][j]; pdata[1] = angular[4][j]; pdata[2] = angular[5][j]; pdata += 3;
      }
      else if (j == -1) /* static obstacle */
      {
	pdata[0] = 0.0; pdata[1] = 0.0; pdata[2] = 0.0; pdata += 3;
      }
      else /* moving obstacle */
      {
	j = -j-2;

        pdata[0] = obsang[3*j]; pdata[1] = obsang[3*j+1]; pdata[2] =  obsang[3*j+2]; pdata += 3;
      }
    }

    hsize_t dims[2] = {num, 3};
    ASSERT (H5LTmake_dataset_double (h5_step, "ANGVEL", 2, dims, data) >= 0, "HDF5 file write error");
  }

  if (ent & OUT_FORCE)
  {
    pdata = data;
    for (i = 0; i < num; i ++)
    {
      j = triobs[set[i]];

      if (j >= 0)
      {
        pdata[0] = force[0][j]; pdata[1] = force[1][j]; pdata[2] = force[2][j]; pdata += 3;
      }
      else
      {
	pdata[0] = 0.0; pdata[1] = 0.0; pdata[2] = 0.0; pdata += 3;
      }
    }

    hsize_t dims[2] = {num, 3};
    ASSERT (H5LTmake_dataset_double (h5_step, "FORCE", 2, dims, data) >= 0, "HDF5 file write error");
  }

  if (ent & OUT_TORQUE)
  {
    pdata = data;
    for (i = 0; i < num; i ++)
    {
      j = triobs[set[i]];

      if (j >= 0)
      {
        pdata[0] = torque[0][j]; pdata[1] = torque[1][j]; pdata[2] = torque[2][j]; pdata += 3;
      }
      else
      {
	pdata[0] = 0.0; pdata[1] = 0.0; pdata[2] = 0.0; pdata += 3;
      }
    }

    hsize_t dims[2] = {num, 3};
    ASSERT (H5LTmake_dataset_double (h5_step, "TORQUE", 2, dims, data) >= 0, "HDF5 file write error");
  }

  delete [] data;
  delete [] topo;
}

#if MED
/* rotation matrix to quaternion conversion:
 * http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/ 
 * to be used along with the med-3.2.0/tests/usecases/c/UsesCase_MEDmesh_9.c example 
 * to support the MEDmeshNodeCoordinateTrsfWr function call */
static void rot2qua (REAL *r, REAL *q)
{
  q[0] = 0.5*sqrt(1. + r[0] + r[4] + r[8]);
  q[1] = (r[5] - r[7])/(4.*q[0]);
  q[2] = (r[6] - r[2])/(4.*q[0]);
  q[3] = (r[1] - r[3])/(4.*q[0]);
}

/* output MED dataset of triangles; based on
 * https://gitlab.onelab.info/gmsh/gmsh/blob/master/Geo/GModelIO_MED.cpp
 * med-3.2.0/tests/usecases/c/UsesCase_MEDmesh_{1,6,9}.c 
 * med-3.2.0/tests/usecases/c/UsesCase_MEDfield_{4}.c */
static void med_triangle_dataset (int num, int *set, int ent, med_idt fid)
{
  int i, j;

  /* mesh name */
  char meshName[MED_NAME_SIZE + 1] = "PARMEC mesh";

  if (output_frame == 0)
  {
    /* create 3d unstructured mesh */
    char dtUnit[MED_SNAME_SIZE + 1] = "";
    char axisName[3 * MED_SNAME_SIZE + 1] = "";
    char axisUnit[3 * MED_SNAME_SIZE + 1] = "";
    ASSERT(MEDmeshCr(fid, meshName, 3, 3, MED_UNSTRUCTURED_MESH, "PARMEC mesh in MED format", dtUnit,
	   MED_SORT_DTIT, MED_CARTESIAN, axisName, axisUnit) >= 0, "Could not create MED mesh");

    /* create family 0 : by default, all mesh entities family number is 0 */
    ASSERT (MEDfamilyCr(fid, meshName, MED_NO_NAME, 0, 0, MED_NO_GROUP) >= 0, "Could not create MED family");
  }

  /* write nodes */
  {
    std::vector<med_float> coord;
    for (i = 0; i < num; i ++)
    {
      j = set[i];
      coord.push_back(tri[0][0][j]);
      coord.push_back(tri[0][1][j]);
      coord.push_back(tri[0][2][j]);
      coord.push_back(tri[1][0][j]);
      coord.push_back(tri[1][1][j]);
      coord.push_back(tri[1][2][j]);
      coord.push_back(tri[2][0][j]);
      coord.push_back(tri[2][1][j]);
      coord.push_back(tri[2][2][j]);
    }

    if (output_frame == 0)
    {
      ASSERT (MEDmeshNodeCoordinateWr(fid, meshName, MED_NO_DT, MED_NO_IT, 0.,
	      MED_FULL_INTERLACE, coord.size()/3, &coord[0]) >= 0, "Could not write MED nodes");
    }
    else
    {
      ASSERT (MEDmeshNodeCoordinateWr(fid, meshName, output_frame, 1, curtime-curtime_output,
	      MED_FULL_INTERLACE, coord.size()/3, &coord[0]) >= 0, "Could not write MED nodes");
    }
  }

  /* write elements */
  if (output_frame == 0)
  {
    std::vector<med_int> conn;
    for (i = 0; i < num; i ++)
    {
      conn.push_back(3*i+1);
      conn.push_back(3*i+2);
      conn.push_back(3*i+3);
    }

    ASSERT (MEDmeshElementConnectivityWr(fid, meshName, MED_NO_DT, MED_NO_IT, 0., MED_CELL, MED_TRIA3,
	    MED_NODAL, MED_FULL_INTERLACE, conn.size()/3, &conn[0]) >= 0, "Could not write MED elements");
  }

  if (ent & OUT_DISPL)
  {
    const char fieldName[MED_NAME_SIZE+1] = "displ";
    const med_int ncomponent = 3;
    const char componentName[3*MED_SNAME_SIZE+1] = "dx              dy              dz";
    const char componentUnit[3*MED_SNAME_SIZE+1] = "1.0             1.0             1.0";
    const char dtUnit[MED_SNAME_SIZE+1] = "s";

    if (output_frame == 0)
    {
      ASSERT (MEDfieldCr(fid, fieldName, MED_FLOAT64, ncomponent, componentName,
        componentUnit, dtUnit, meshName) >= 0, "Could not create MED field");
    }

    std::vector<med_float> values;
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
	values.push_back (d[0]);
	values.push_back (d[1]);
	values.push_back (d[2]);

	SUB (p2, x, d);
	TVADDMUL (X, L, d, P);
	SUB (p2, P, d);
	values.push_back (d[0]);
	values.push_back (d[1]);
	values.push_back (d[2]);

	SUB (p3, x, d);
	TVADDMUL (X, L, d, P);
	SUB (p3, P, d);
	values.push_back (d[0]);
	values.push_back (d[1]);
	values.push_back (d[2]);
      }
      else if (j == -1) /* static obstacle */
      {
	values.push_back (0.);
	values.push_back (0.);
	values.push_back (0.);
	values.push_back (0.);
	values.push_back (0.);
	values.push_back (0.);
	values.push_back (0.);
	values.push_back (0.);
	values.push_back (0.);
      }
      else /* moving obstacle */
      {
	j = -j-2;

	REAL x[3] = {obspnt[3*j], obspnt[3*j+1], obspnt[3*j+2]}, d[3];

	SUB (p1, x, d);
	values.push_back (d[0]);
	values.push_back (d[1]);
	values.push_back (d[2]);
	SUB (p2, x, d);
	values.push_back (d[0]);
	values.push_back (d[1]);
	values.push_back (d[2]);
	SUB (p3, x, d);
	values.push_back (d[0]);
	values.push_back (d[1]);
	values.push_back (d[2]);
      }
    }

    if (output_frame == 0)
    {
      ASSERT (MEDfieldValueWr(fid, fieldName, MED_NO_DT, MED_NO_IT, 0., MED_NODE, MED_NONE, MED_FULL_INTERLACE, 
        MED_ALL_CONSTITUENT, values.size()/3, (unsigned char*) &values[0]) >= 0, "Could not write MED field values");
    }
    else
    {
      ASSERT (MEDfieldValueWr(fid, fieldName, output_frame, 1, curtime-curtime_output, MED_NODE, MED_NONE, MED_FULL_INTERLACE, 
        MED_ALL_CONSTITUENT, values.size()/3, (unsigned char*) &values[0]) >= 0, "Could not write MED field values");
    }
  }

  if (ent & OUT_LINVEL)
  {
    const char fieldName[MED_NAME_SIZE+1] = "linvel";
    const med_int ncomponent = 3;
    const char componentName[3*MED_SNAME_SIZE+1] = "vx              vy              vz";
    const char componentUnit[3*MED_SNAME_SIZE+1] = "1.0             1.0             1.0";
    const char dtUnit[MED_SNAME_SIZE+1] = "s";

    if (output_frame == 0)
    {
      ASSERT (MEDfieldCr(fid, fieldName, MED_FLOAT64, ncomponent, componentName,
        componentUnit, dtUnit, meshName) >= 0, "Could not create MED field");
    }

    std::vector<med_float> values;
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
	values.push_back(w[0]);
	values.push_back(w[1]);
	values.push_back(w[2]);

        COPY (v, w);
	SUB (p2, x, a);
	PRODUCTADD (o, a, w);
	values.push_back(w[0]);
	values.push_back(w[1]);
	values.push_back(w[2]);

        COPY (v, w);
	SUB (p3, x, a);
	PRODUCTADD (o, a, w);
	values.push_back(w[0]);
	values.push_back(w[1]);
	values.push_back(w[2]);
      }
      else if (j == -1) /* static obstacle */
      {
	values.push_back(0.);
	values.push_back(0.);
	values.push_back(0.);
	values.push_back(0.);
	values.push_back(0.);
	values.push_back(0.);
	values.push_back(0.);
	values.push_back(0.);
	values.push_back(0.);
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
	values.push_back(w[0]);
	values.push_back(w[1]);
	values.push_back(w[2]);

        COPY (v, w);
	SUB (p2, x, a);
	PRODUCTADD (o, a, w);
	values.push_back(w[0]);
	values.push_back(w[1]);
	values.push_back(w[2]);

        COPY (v, w);
	SUB (p3, x, a);
	PRODUCTADD (o, a, w);
	values.push_back(w[0]);
	values.push_back(w[1]);
	values.push_back(w[2]);
      }
    }

    if (output_frame == 0)
    {
      ASSERT (MEDfieldValueWr(fid, fieldName, MED_NO_DT, MED_NO_IT, 0., MED_NODE, MED_NONE, MED_FULL_INTERLACE, 
        MED_ALL_CONSTITUENT, values.size()/3, (unsigned char*) &values[0]) >= 0, "Could not write MED field values");
    }
    else
    {
      ASSERT (MEDfieldValueWr(fid, fieldName, output_frame, 1, curtime-curtime_output, MED_NODE, MED_NONE, MED_FULL_INTERLACE, 
        MED_ALL_CONSTITUENT, values.size()/3, (unsigned char*) &values[0]) >= 0, "Could not write MED field values");
    }
  }

  if (ent & OUT_NUMBER)
  {
    const char fieldName[MED_NAME_SIZE+1] = "number";
    const med_int ncomponent = 1;
    const char componentName[MED_SNAME_SIZE+1] = "number";
    const char componentUnit[MED_SNAME_SIZE+1] = "1";
    const char dtUnit[MED_SNAME_SIZE+1] = "s";

    if (output_frame == 0)
    {
      ASSERT (MEDfieldCr(fid, fieldName, MED_FLOAT64, ncomponent, componentName,
        componentUnit, dtUnit, meshName) >= 0, "Could not create MED field");
    }

    std::vector<med_float> values;
    for (i = 0; i < num; i ++)
    {
      values.push_back(triobs[set[i]]);
    }

    if (output_frame == 0)
    {
      ASSERT (MEDfieldValueWr(fid, fieldName, MED_NO_DT, MED_NO_IT, 0., MED_CELL, MED_TRIA3, MED_FULL_INTERLACE, 
        MED_ALL_CONSTITUENT, values.size(), (unsigned char*) &values[0]) >= 0, "Could not write MED field values");
    }
    else
    {
      ASSERT (MEDfieldValueWr(fid, fieldName, output_frame, 1, curtime-curtime_output, MED_CELL, MED_TRIA3, MED_FULL_INTERLACE, 
        MED_ALL_CONSTITUENT, values.size(), (unsigned char*) &values[0]) >= 0, "Could not write MED field values");
    }
  }

  if (ent & OUT_COLOR)
  {
    const char fieldName[MED_NAME_SIZE+1] = "color";
    const med_int ncomponent = 1;
    const char componentName[MED_SNAME_SIZE+1] = "color";
    const char componentUnit[MED_SNAME_SIZE+1] = "1";
    const char dtUnit[MED_SNAME_SIZE+1] = "s";

    if (output_frame == 0)
    {
      ASSERT (MEDfieldCr(fid, fieldName, MED_FLOAT64, ncomponent, componentName,
        componentUnit, dtUnit, meshName) >= 0, "Could not create MED field");
    }

    std::vector<med_float> values;
    for (i = 0; i < num; i ++)
    {
      values.push_back(tricol[set[i]]);
    }

    if (output_frame == 0)
    {
      ASSERT (MEDfieldValueWr(fid, fieldName, MED_NO_DT, MED_NO_IT, 0., MED_CELL, MED_TRIA3, MED_FULL_INTERLACE, 
        MED_ALL_CONSTITUENT, values.size(), (unsigned char*) &values[0]) >= 0, "Could not write MED field values");
    }
    else
    {
      ASSERT (MEDfieldValueWr(fid, fieldName, output_frame, 1, curtime-curtime_output, MED_CELL, MED_TRIA3, MED_FULL_INTERLACE, 
        MED_ALL_CONSTITUENT, values.size(), (unsigned char*) &values[0]) >= 0, "Could not write MED field values");
    }
  }

  if (ent & OUT_ANGVEL)
  {
    const char fieldName[MED_NAME_SIZE+1] = "angvel";
    const med_int ncomponent = 3;
    const char componentName[3*MED_SNAME_SIZE+1] = "ox              oy              oz";
    const char componentUnit[3*MED_SNAME_SIZE+1] = "1.0             1.0             1.0";
    const char dtUnit[MED_SNAME_SIZE+1] = "s";

    if (output_frame == 0)
    {
      ASSERT (MEDfieldCr(fid, fieldName, MED_FLOAT64, ncomponent, componentName,
        componentUnit, dtUnit, meshName) >= 0, "Could not create MED field");
    }

    std::vector<med_float> values;
    for (i = 0; i < num; i ++)
    {
      j = triobs[set[i]];

      if (j >= 0) /* particle */
      {
        values.push_back(angular[3][j]);
	values.push_back(angular[4][j]);
	values.push_back(angular[5][j]);
      }
      else if (j == -1) /* static obstacle */
      {
        values.push_back(0.);
	values.push_back(0.);
	values.push_back(0.);
      }
      else /* moving obstacle */
      {
	j = -j-2;

        values.push_back(obsang[3*j]); 
	values.push_back(obsang[3*j+1]);
	values.push_back(obsang[3*j+2]);
      }
    }

    if (output_frame == 0)
    {
      ASSERT (MEDfieldValueWr(fid, fieldName, MED_NO_DT, MED_NO_IT, 0., MED_CELL, MED_TRIA3, MED_FULL_INTERLACE, 
        MED_ALL_CONSTITUENT, values.size()/3, (unsigned char*) &values[0]) >= 0, "Could not write MED field values");
    }
    else
    {
      ASSERT (MEDfieldValueWr(fid, fieldName, output_frame, 1, curtime-curtime_output, MED_CELL, MED_TRIA3, MED_FULL_INTERLACE, 
        MED_ALL_CONSTITUENT, values.size()/3, (unsigned char*) &values[0]) >= 0, "Could not write MED field values");
    }
  }

  if (ent & OUT_FORCE)
  {
    const char fieldName[MED_NAME_SIZE+1] = "force";
    const med_int ncomponent = 3;
    const char componentName[3*MED_SNAME_SIZE+1] = "fx              fy              fz";
    const char componentUnit[3*MED_SNAME_SIZE+1] = "1.0             1.0             1.0";
    const char dtUnit[MED_SNAME_SIZE+1] = "s";

    if (output_frame == 0)
    {
      ASSERT (MEDfieldCr(fid, fieldName, MED_FLOAT64, ncomponent, componentName,
        componentUnit, dtUnit, meshName) >= 0, "Could not create MED field");
    }

    std::vector<med_float> values;
    for (i = 0; i < num; i ++)
    {
      j = triobs[set[i]];

      if (j >= 0) /* particle */
      {
        values.push_back(force[0][j]); 
	values.push_back(force[1][j]); 
	values.push_back(force[2][j]);
      }
      else /* obstacle */
      {
        values.push_back(0.);
	values.push_back(0.);
	values.push_back(0.);
      }
    }

    if (output_frame == 0)
    {
      ASSERT (MEDfieldValueWr(fid, fieldName, MED_NO_DT, MED_NO_IT, 0., MED_CELL, MED_TRIA3, MED_FULL_INTERLACE, 
        MED_ALL_CONSTITUENT, values.size()/3, (unsigned char*) &values[0]) >= 0, "Could not write MED field values");
    }
    else
    {
      ASSERT (MEDfieldValueWr(fid, fieldName, output_frame, 1, curtime-curtime_output, MED_CELL, MED_TRIA3, MED_FULL_INTERLACE, 
        MED_ALL_CONSTITUENT, values.size()/3, (unsigned char*) &values[0]) >= 0, "Could not write MED field values");
    }
  }

  if (ent & OUT_TORQUE)
  {
    const char fieldName[MED_NAME_SIZE+1] = "torque";
    const med_int ncomponent = 3;
    const char componentName[3*MED_SNAME_SIZE+1] = "tx              ty              tz";
    const char componentUnit[3*MED_SNAME_SIZE+1] = "1.0             1.0             1.0";
    const char dtUnit[MED_SNAME_SIZE+1] = "s";

    if (output_frame == 0)
    {
      ASSERT (MEDfieldCr(fid, fieldName, MED_FLOAT64, ncomponent, componentName,
        componentUnit, dtUnit, meshName) >= 0, "Could not create MED field");
    }

    std::vector<med_float> values;
    for (i = 0; i < num; i ++)
    {
      j = triobs[set[i]];

      if (j >= 0) /* particle */
      {
        values.push_back(torque[0][j]); 
	values.push_back(torque[1][j]); 
	values.push_back(torque[2][j]);
      }
      else /* obstacle */
      {
        values.push_back(0.);
	values.push_back(0.);
	values.push_back(0.);
      }
    }

    if (output_frame == 0)
    {
      ASSERT (MEDfieldValueWr(fid, fieldName, MED_NO_DT, MED_NO_IT, 0., MED_CELL, MED_TRIA3, MED_FULL_INTERLACE, 
        MED_ALL_CONSTITUENT, values.size()/3, (unsigned char*) &values[0]) >= 0, "Could not write MED field values");
    }
    else
    {
      ASSERT (MEDfieldValueWr(fid, fieldName, output_frame, 1, curtime-curtime_output, MED_CELL, MED_TRIA3, MED_FULL_INTERLACE, 
        MED_ALL_CONSTITUENT, values.size()/3, (unsigned char*) &values[0]) >= 0, "Could not write MED field values");
    }
  }
}
#endif

/* find ellipsoids belonging to a particle set [part0, part1) */
static int find_ellipsoid_set (int *part0, int *part1, int *ellipsoids)
{
  int num = 0;

  for (int *p = part0; p < part1; p ++)
  {
    int s = 0, e = ellnum, i, j;

    while (s < e)
    {
      i = (s+e)/2;

      j = part[i];

      if (j == *p) break;
      else if (j < *p) s = i;
      else e = i;
    }

    if (j == *p)
    {
      while (i > 0 && part[i-1] == j) i --; /* find start of range */

      while (i < ellnum && part[i] == j)
      {
	ellipsoids[num] = i;
	num ++;
	i ++;
      }
    }
  }

  return num;
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

  if (ent & (OUT_NUMBER|OUT_DISPL|OUT_F|OUT_SF|OUT_PAIR|OUT_SS))
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

  if (ent & OUT_LENGTH)
  {
    out << "SCALARS length float\n";
    out << "LOOKUP_TABLE default\n";
    for (i = 0; i < num; i ++)
    {
      j = set[i];
      REAL q[3] = {sprpnt[1][0][j]-sprpnt[0][0][j],
		   sprpnt[1][1][j]-sprpnt[0][1][j],
		   sprpnt[1][2][j]-sprpnt[0][2][j]};
      out << LEN(q) << "\n";
    }
  }

  if (ent & OUT_ORIENT)
  {
    out << "VECTORS orient float\n";
    for (i = 0; i < num; i ++)
    {
      j = set[i];
      out << sprdir[0][j] << " "
	  << sprdir[1][j] << " "
	  << sprdir[2][j] << "\n";
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

  if (ent & OUT_SS)
  {
    out << "SCALARS SS float\n";
    out << "LOOKUP_TABLE default\n";
    for (i = 0; i < num; i ++)
    {
      j = set[i];
      out << (float)unspring[j] << "\n";
    }
  }
}

/* output hdf5 dataset of spring data */
static void h5_spring_dataset (int num, int *set, int ent, hid_t h5_step)
{
  double *data, *pdata;
  int i, j, *numb;

  ERRMEM (data = new double [3*num]);
  for (i = 0; i < num; i ++)
  {
    j = set[i];
    data[3*i+0] = 0.5*(sprpnt[0][0][j]+sprpnt[1][0][j]);
    data[3*i+1] = 0.5*(sprpnt[0][1][j]+sprpnt[1][1][j]);
    data[3*i+2] = 0.5*(sprpnt[0][2][j]+sprpnt[1][2][j]);
  }
  hsize_t dims[2] = {num, 3};
  ASSERT (H5LTmake_dataset_double (h5_step, "GEOM", 2, dims, data) >= 0, "HDF5 file write error");

  if (ent & OUT_NUMBER)
  {
    ERRMEM (numb = new int[num]);

    for (i = 0; i < num; i ++)
    {
      j = set[i];
      numb[i] = sprid[j];
    }

    hsize_t length = num;
    ASSERT (H5LTmake_dataset_int (h5_step, "NUMBER", 1, &length, numb) >= 0, "HDF5 file write error");

    delete [] numb;
  }

  if (ent & OUT_DISPL)
  {
    for (i = 0; i < num; i ++)
    {
      j = set[i];
      data[i] = stroke[0][j];
    }

    hsize_t length = num;
    ASSERT (H5LTmake_dataset_double (h5_step, "DISPL", 1, &length, data) >= 0, "HDF5 file write error");
  }

  if (ent & OUT_LENGTH)
  {
    for (i = 0; i < num; i ++)
    {
      j = set[i];
      REAL q[3] = {sprpnt[1][0][j]-sprpnt[0][0][j],
		   sprpnt[1][1][j]-sprpnt[0][1][j],
		   sprpnt[1][2][j]-sprpnt[0][2][j]};
      data[i] = LEN(q);
    }

    hsize_t length = num;
    ASSERT (H5LTmake_dataset_double (h5_step, "LENGTH", 1, &length, data) >= 0, "HDF5 file write error");
  }

  if (ent & OUT_ORIENT)
  {
    for (pdata = data, i = 0; i < num; i ++, pdata += 3)
    {
      j = set[i];
      pdata[0] = sprdir[0][j];
      pdata[1] = sprdir[1][j];
      pdata[2] = sprdir[2][j];
    }

    ASSERT (H5LTmake_dataset_double (h5_step, "ORIENT", 2, dims, data) >= 0, "HDF5 file write error");
  }

  if (ent & OUT_F)
  {
    for (i = 0; i < num; i ++)
    {
      j = set[i];
      data[i] = sprfrc[0][j];
    }

    hsize_t length = num;
    ASSERT (H5LTmake_dataset_double (h5_step, "F", 1, &length, data) >= 0, "HDF5 file write error");
  }

  if (ent & OUT_SF)
  {
    for (i = 0; i < num; i ++)
    {
      j = set[i];
      data[i] = sprfrc[1][j];
    }

    hsize_t length = num;
    ASSERT (H5LTmake_dataset_double (h5_step, "SF", 1, &length, data) >= 0, "HDF5 file write error");
  }

  if (ent & OUT_SS)
  {
    for (i = 0; i < num; i ++)
    {
      j = set[i];
      data[i] = unspring[j];
    }

    hsize_t length = num;
    ASSERT (H5LTmake_dataset_double (h5_step, "SS", 1, &length, data) >= 0, "HDF5 file write error");
  }

  delete [] data;
}

/* append an XMF file */
static void append_xmf_file (const char *xmf_path, int mode, int elements, int nodes, int topo_size, const char *label, const char *h5file, int ent)
{
  FILE *xmf_file;

  if (curtime == 0.0)
  {
    ASSERT(xmf_file = fopen (xmf_path, "w"), "XMF markup file open failed");

    fprintf (xmf_file, "<Xdmf>\n");
    fprintf (xmf_file, "<Domain>\n");
    fprintf (xmf_file, "<Grid GridType=\"Collection\" CollectionType=\"Temporal\">\n");
  }
  else
  {
    ASSERT (xmf_file = fopen (xmf_path, "r+"), "XMF markup file open failed");
    fseek (xmf_file, -1, SEEK_END);
    int ln = 0;
    for(;;) /* move back three lines */
    {
      if (fgetc(xmf_file) == '\n') ln ++;
      if (ln == 4) break;
      else fseek (xmf_file, -2, SEEK_CUR);
    }
    int pos = ftell(xmf_file);
    char *mem = (char*)malloc(pos);
    ERRMEM (mem);
    fseek (xmf_file, 0, SEEK_SET);
    ASSERT (fread (mem, sizeof(char), pos, xmf_file) == pos, "XMF markup file read failed"); /* read until the last three lines */
    fclose (xmf_file);
    ASSERT (xmf_file = fopen (xmf_path, "w"), "XMF markup file open failed");
    fwrite (mem, sizeof(char), pos, xmf_file); /* effectively truncate the last three lines */
    free (mem);
  }

  ASSERT (xmf_file, "XDMF file open error");

  fprintf (xmf_file, "<Grid Name=\"%s\" Type=\"Uniform\">\n", label);
  fprintf (xmf_file, "<Time Type=\"Single\" Value=\"%f\" />\n", curtime);

  switch (mode)
  {
  case OUT_MODE_MESH:
    fprintf (xmf_file, "<Topology Type=\"Mixed\" NumberOfElements=\"%d\">\n", elements);
    fprintf (xmf_file, "<DataStructure Dimensions=\"%d\" NumberType=\"Int\" Format=\"HDF\">\n", topo_size);
    fprintf (xmf_file, "%s:/%d/TOPO\n", h5file, output_frame);
    fprintf (xmf_file, "</DataStructure>\n");
    fprintf (xmf_file, "</Topology>\n");
    break;
  case OUT_MODE_RB:
  case OUT_MODE_SD:
    fprintf (xmf_file, "<Topology Type=\"Polyvertex\" NumberOfElements=\"%d\" NodesPerElement=\"%d\">\n", nodes, 1);
    fprintf (xmf_file, "</Topology>\n");
    break;
  }

  fprintf (xmf_file, "<Geometry GeometryType=\"XYZ\">\n");
  fprintf (xmf_file, "<DataStructure Dimensions=\"%d 3\" NumberType=\"Float\" Presicion=\"8\" Format=\"HDF\">\n", nodes);
  fprintf (xmf_file, "%s:/%d/GEOM\n", h5file, output_frame);
  fprintf (xmf_file, "</DataStructure>\n");
  fprintf (xmf_file, "</Geometry>\n");

  switch (mode)
  {
  case OUT_MODE_MESH:

    if (ent & OUT_DISPL)
    {
      fprintf (xmf_file, "<Attribute Name=\"DISPL\" Center=\"Node\" AttributeType=\"Vector\">\n");
      fprintf (xmf_file, "<DataStructure Dimensions=\"%d 3\" NumberType=\"Float\" Presicion=\"8\" Format=\"HDF\">\n", nodes);
      fprintf (xmf_file, "%s:/%d/DISPL\n", h5file, output_frame);
      fprintf (xmf_file, "</DataStructure>\n");
      fprintf (xmf_file, "</Attribute>\n");
    }

    if (ent & OUT_LINVEL)
    {
      fprintf (xmf_file, "<Attribute Name=\"LINVEL\" Center=\"Node\" AttributeType=\"Vector\">\n");
      fprintf (xmf_file, "<DataStructure Dimensions=\"%d 3\" NumberType=\"Float\" Presicion=\"8\" Format=\"HDF\">\n", nodes);
      fprintf (xmf_file, "%s:/%d/LINVEL\n", h5file, output_frame);
      fprintf (xmf_file, "</DataStructure>\n");
      fprintf (xmf_file, "</Attribute>\n");
    }

    if (ent & OUT_NUMBER)
    {
      fprintf (xmf_file, "<Attribute Name=\"NUMBER\" Center=\"Cell\" AttributeType=\"Scalar\">\n");
      fprintf (xmf_file, "<DataStructure Dimensions=\"%d\" NumberType=\"Int\" Format=\"HDF\">\n", elements);
      fprintf (xmf_file, "%s:/%d/NUMBER\n", h5file, output_frame);
      fprintf (xmf_file, "</DataStructure>\n");
      fprintf (xmf_file, "</Attribute>\n");
    }

    if (ent & OUT_COLOR)
    {
      fprintf (xmf_file, "<Attribute Name=\"COLOR\" Center=\"Cell\" AttributeType=\"Scalar\">\n");
      fprintf (xmf_file, "<DataStructure Dimensions=\"%d\" NumberType=\"Int\" Format=\"HDF\">\n", elements);
      fprintf (xmf_file, "%s:/%d/COLOR\n", h5file, output_frame);
      fprintf (xmf_file, "</DataStructure>\n");
      fprintf (xmf_file, "</Attribute>\n");
    }

    if (ent & OUT_ANGVEL)
    {
      fprintf (xmf_file, "<Attribute Name=\"ANGVEL\" Center=\"Cell\" AttributeType=\"Vector\">\n");
      fprintf (xmf_file, "<DataStructure Dimensions=\"%d 3\" NumberType=\"Float\" Presicion=\"8\" Format=\"HDF\">\n", elements);
      fprintf (xmf_file, "%s:/%d/ANGVEL\n", h5file, output_frame);
      fprintf (xmf_file, "</DataStructure>\n");
      fprintf (xmf_file, "</Attribute>\n");
    }

    if (ent & OUT_FORCE)
    {
      fprintf (xmf_file, "<Attribute Name=\"FORCE\" Center=\"Cell\" AttributeType=\"Vector\">\n");
      fprintf (xmf_file, "<DataStructure Dimensions=\"%d 3\" NumberType=\"Float\" Presicion=\"8\" Format=\"HDF\">\n", elements);
      fprintf (xmf_file, "%s:/%d/FORCE\n", h5file, output_frame);
      fprintf (xmf_file, "</DataStructure>\n");
      fprintf (xmf_file, "</Attribute>\n");
    }

    if (ent & OUT_TORQUE)
    {
      fprintf (xmf_file, "<Attribute Name=\"TORQUE\" Center=\"Cell\" AttributeType=\"Vector\">\n");
      fprintf (xmf_file, "<DataStructure Dimensions=\"%d 3\" NumberType=\"Float\" Presicion=\"8\" Format=\"HDF\">\n", elements);
      fprintf (xmf_file, "%s:/%d/TORQUE\n", h5file, output_frame);
      fprintf (xmf_file, "</DataStructure>\n");
      fprintf (xmf_file, "</Attribute>\n");
    }
    break;

  case OUT_MODE_RB:

    if (ent & OUT_DISPL)
    {
      fprintf (xmf_file, "<Attribute Name=\"DISPL\" Center=\"Node\" AttributeType=\"Vector\">\n");
      fprintf (xmf_file, "<DataStructure Dimensions=\"%d 3\" NumberType=\"Float\" Presicion=\"8\" Format=\"HDF\">\n", nodes);
      fprintf (xmf_file, "%s:/%d/DISPL\n", h5file, output_frame);
      fprintf (xmf_file, "</DataStructure>\n");
      fprintf (xmf_file, "</Attribute>\n");
    }

    if (ent & OUT_LINVEL)
    {
      fprintf (xmf_file, "<Attribute Name=\"LINVEL\" Center=\"Node\" AttributeType=\"Vector\">\n");
      fprintf (xmf_file, "<DataStructure Dimensions=\"%d 3\" NumberType=\"Float\" Presicion=\"8\" Format=\"HDF\">\n", nodes);
      fprintf (xmf_file, "%s:/%d/LINVEL\n", h5file, output_frame);
      fprintf (xmf_file, "</DataStructure>\n");
      fprintf (xmf_file, "</Attribute>\n");
    }

    if (ent & OUT_NUMBER)
    {
      fprintf (xmf_file, "<Attribute Name=\"NUMBER\" Center=\"Node\" AttributeType=\"Scalar\">\n");
      fprintf (xmf_file, "<DataStructure Dimensions=\"%d\" NumberType=\"Int\" Format=\"HDF\">\n", nodes);
      fprintf (xmf_file, "%s:/%d/NUMBER\n", h5file, output_frame);
      fprintf (xmf_file, "</DataStructure>\n");
      fprintf (xmf_file, "</Attribute>\n");
    }

    if (ent & OUT_ANGVEL)
    {
      fprintf (xmf_file, "<Attribute Name=\"ANGVEL\" Center=\"Node\" AttributeType=\"Vector\">\n");
      fprintf (xmf_file, "<DataStructure Dimensions=\"%d 3\" NumberType=\"Float\" Presicion=\"8\" Format=\"HDF\">\n", nodes);
      fprintf (xmf_file, "%s:/%d/ANGVEL\n", h5file, output_frame);
      fprintf (xmf_file, "</DataStructure>\n");
      fprintf (xmf_file, "</Attribute>\n");
    }

    if (ent & OUT_FORCE)
    {
      fprintf (xmf_file, "<Attribute Name=\"FORCE\" Center=\"Node\" AttributeType=\"Vector\">\n");
      fprintf (xmf_file, "<DataStructure Dimensions=\"%d 3\" NumberType=\"Float\" Presicion=\"8\" Format=\"HDF\">\n", nodes);
      fprintf (xmf_file, "%s:/%d/FORCE\n", h5file, output_frame);
      fprintf (xmf_file, "</DataStructure>\n");
      fprintf (xmf_file, "</Attribute>\n");
    }

    if (ent & OUT_TORQUE)
    {
      fprintf (xmf_file, "<Attribute Name=\"TORQUE\" Center=\"Node\" AttributeType=\"Vector\">\n");
      fprintf (xmf_file, "<DataStructure Dimensions=\"%d 3\" NumberType=\"Float\" Presicion=\"8\" Format=\"HDF\">\n", nodes);
      fprintf (xmf_file, "%s:/%d/TORQUE\n", h5file, output_frame);
      fprintf (xmf_file, "</DataStructure>\n");
      fprintf (xmf_file, "</Attribute>\n");
    }

    if (ent & OUT_ORIENT)
    {
      fprintf (xmf_file, "<Attribute Name=\"ORIENT\" Center=\"Node\" AttributeType=\"Tensor\">\n");
      fprintf (xmf_file, "<DataStructure Dimensions=\"%d 9\" NumberType=\"Float\" Presicion=\"8\" Format=\"HDF\">\n", nodes);
      fprintf (xmf_file, "%s:/%d/ORIENT\n", h5file, output_frame);
      fprintf (xmf_file, "</DataStructure>\n");
      fprintf (xmf_file, "</Attribute>\n");
    }

    if (ent & OUT_ORIENT1)
    {
      fprintf (xmf_file, "<Attribute Name=\"ORIENT1\" Center=\"Node\" AttributeType=\"Vector\">\n");
      fprintf (xmf_file, "<DataStructure Dimensions=\"%d 3\" NumberType=\"Float\" Presicion=\"8\" Format=\"HDF\">\n", nodes);
      fprintf (xmf_file, "%s:/%d/ORIENT1\n", h5file, output_frame);
      fprintf (xmf_file, "</DataStructure>\n");
      fprintf (xmf_file, "</Attribute>\n");
    }

    if (ent & OUT_ORIENT2)
    {
      fprintf (xmf_file, "<Attribute Name=\"ORIENT2\" Center=\"Node\" AttributeType=\"Vector\">\n");
      fprintf (xmf_file, "<DataStructure Dimensions=\"%d 3\" NumberType=\"Float\" Presicion=\"8\" Format=\"HDF\">\n", nodes);
      fprintf (xmf_file, "%s:/%d/ORIENT2\n", h5file, output_frame);
      fprintf (xmf_file, "</DataStructure>\n");
      fprintf (xmf_file, "</Attribute>\n");
    }

    if (ent & OUT_ORIENT3)
    {
      fprintf (xmf_file, "<Attribute Name=\"ORIENT3\" Center=\"Node\" AttributeType=\"Vector\">\n");
      fprintf (xmf_file, "<DataStructure Dimensions=\"%d 3\" NumberType=\"Float\" Presicion=\"8\" Format=\"HDF\">\n", nodes);
      fprintf (xmf_file, "%s:/%d/ORIENT3\n", h5file, output_frame);
      fprintf (xmf_file, "</DataStructure>\n");
      fprintf (xmf_file, "</Attribute>\n");
    }
    break;
  case OUT_MODE_SD:

    if (ent & OUT_NUMBER)
    {
      fprintf (xmf_file, "<Attribute Name=\"NUMBER\" Center=\"Node\" AttributeType=\"Scalar\">\n");
      fprintf (xmf_file, "<DataStructure Dimensions=\"%d\" NumberType=\"Int\" Format=\"HDF\">\n", nodes);
      fprintf (xmf_file, "%s:/%d/NUMBER\n", h5file, output_frame);
      fprintf (xmf_file, "</DataStructure>\n");
      fprintf (xmf_file, "</Attribute>\n");
    }

    if (ent & OUT_DISPL)
    {
      fprintf (xmf_file, "<Attribute Name=\"DISPL\" Center=\"Node\" AttributeType=\"Scalar\">\n");
      fprintf (xmf_file, "<DataStructure Dimensions=\"%d\" NumberType=\"Float\" Presicion=\"8\" Format=\"HDF\">\n", nodes);
      fprintf (xmf_file, "%s:/%d/DISPL\n", h5file, output_frame);
      fprintf (xmf_file, "</DataStructure>\n");
      fprintf (xmf_file, "</Attribute>\n");
    }

    if (ent & OUT_LENGTH)
    {
      fprintf (xmf_file, "<Attribute Name=\"LENGTH\" Center=\"Node\" AttributeType=\"Scalar\">\n");
      fprintf (xmf_file, "<DataStructure Dimensions=\"%d\" NumberType=\"Float\" Presicion=\"8\" Format=\"HDF\">\n", nodes);
      fprintf (xmf_file, "%s:/%d/LENGTH\n", h5file, output_frame);
      fprintf (xmf_file, "</DataStructure>\n");
      fprintf (xmf_file, "</Attribute>\n");
    }

    if (ent & OUT_ORIENT)
    {
      fprintf (xmf_file, "<Attribute Name=\"ORIENT\" Center=\"Node\" AttributeType=\"Vector\">\n");
      fprintf (xmf_file, "<DataStructure Dimensions=\"%d 3\" NumberType=\"Float\" Presicion=\"8\" Format=\"HDF\">\n", nodes);
      fprintf (xmf_file, "%s:/%d/ORIENT\n", h5file, output_frame);
      fprintf (xmf_file, "</DataStructure>\n");
      fprintf (xmf_file, "</Attribute>\n");
    }

    if (ent & OUT_F)
    {
      fprintf (xmf_file, "<Attribute Name=\"F\" Center=\"Node\" AttributeType=\"Scalar\">\n");
      fprintf (xmf_file, "<DataStructure Dimensions=\"%d\" NumberType=\"Float\" Presicion=\"8\" Format=\"HDF\">\n", nodes);
      fprintf (xmf_file, "%s:/%d/F\n", h5file, output_frame);
      fprintf (xmf_file, "</DataStructure>\n");
      fprintf (xmf_file, "</Attribute>\n");
    }

    if (ent & OUT_SF)
    {
      fprintf (xmf_file, "<Attribute Name=\"SF\" Center=\"Node\" AttributeType=\"Scalar\">\n");
      fprintf (xmf_file, "<DataStructure Dimensions=\"%d\" NumberType=\"Float\" Presicion=\"8\" Format=\"HDF\">\n", nodes);
      fprintf (xmf_file, "%s:/%d/SF\n", h5file, output_frame);
      fprintf (xmf_file, "</DataStructure>\n");
      fprintf (xmf_file, "</Attribute>\n");
    }

    if (ent & OUT_SS)
    {
      fprintf (xmf_file, "<Attribute Name=\"SS\" Center=\"Node\" AttributeType=\"Scalar\">\n");
      fprintf (xmf_file, "<DataStructure Dimensions=\"%d\" NumberType=\"Float\" Presicion=\"8\" Format=\"HDF\">\n", nodes);
      fprintf (xmf_file, "%s:/%d/SS\n", h5file, output_frame);
      fprintf (xmf_file, "</DataStructure>\n");
      fprintf (xmf_file, "</Attribute>\n");
    }
    break;
  }

  fprintf (xmf_file, "</Grid>\n");
  fprintf (xmf_file, "</Grid>\n");
  fprintf (xmf_file, "</Domain>\n");
  fprintf (xmf_file, "</Xdmf>\n");

  fclose (xmf_file);
}

/* output XDMF files */
static void output_xdmf_files ()
{
  ostringstream h5_path, h5_text, xmf_path;
  FILE *xmf_file;
  hid_t h5_file;
  hid_t h5_step;

  if (trinum)
  {
    int i, j, num, ent, *set; /* number of and the set of triangles */

    ERRMEM (set = new int[trinum]);

    for (j = -1; j < outnum; j ++) /* for each output set */
    {
      num = 0;

      if (j < 0 && (outrest[1] & OUT_MODE_MESH)) /* output unselected triangles */
      {
	for (i = 0; i < trinum; i ++)
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

	ent = outrest[0];
      }
      else if (outmode[j] & OUT_MODE_MESH) /* output selected triangles */
      {
	num = find_triangle_set (&outpart[outidx[j]], &outpart[outidx[j+1]], set);

	ent = outent[j];
      }

      if (num)
      { 
	h5_path.str("");
	h5_path.clear();
	h5_path << output_path << j+1 << ".h5";

	if (curtime == 0.0)
	{
	  ASSERT ((h5_file = H5Fcreate(h5_path.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)) >= 0, "HDF5 file open error");
	}
	else
	{
	  ASSERT((h5_file = H5Fopen(h5_path.str().c_str(), H5F_ACC_RDWR, H5P_DEFAULT)) >= 0, "HDF5 file open error");
	}

	h5_text.str("");
	h5_text.clear();
	h5_text << output_frame;
	ASSERT ((h5_step = H5Gcreate (h5_file, h5_text.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) >= 0, "HDF5 file write error");
	double time = curtime;
        ASSERT (H5LTset_attribute_double (h5_step, ".", "TIME", &time, 1) >= 0, "HDF5 file write error"); 

	h5_triangle_dataset (num, set, ent, h5_step); /* append h5 file */

        xmf_path.str("");
	xmf_path.clear();
	xmf_path << output_path << j+1 << ".xmf"; /* append xmf file */

	int elements = num;
	int nodes = 3*num;
	int topo_size = 4*num;
	const char *label = "PARMEC triangles";
	string h5file = h5_path.str().substr(h5_path.str().find_last_of('/')+1);

        append_xmf_file (xmf_path.str().c_str(), OUT_MODE_MESH, elements, nodes, topo_size, label, h5file.c_str(), ent);

	H5Gclose (h5_step);
	H5Fclose (h5_file);
      }
    }

    delete [] set;
  }

  if (parnum)
  {
    int i, j, num, ent, *set, *pset; /* number of and the set of unselected particles */

    ERRMEM (set = new int[parnum]);

    for (j = -1; j < outnum; j ++)
    {
      num = 0;

      if (j < 0 && (outrest[1] & OUT_MODE_RB)) /* output unselected particles */
      {
	for (i = 0; i < parnum; i ++)
	{
	  if (flags[i] & OUTREST)
	  {
	    set[num ++] = i;
	  }
	}

	pset = set;
	ent = outrest[0];
      }
      else if (outmode[j] & OUT_MODE_RB) /* output selected particles */
      {
	num = outidx[j+1]-outidx[j];
	pset = &outpart[outidx[j]];
	ent = outent[j];
      }

      if (num)
      {
	h5_path.str("");
	h5_path.clear();
	h5_path << output_path << j+1 << "rb.h5";

	if (curtime == 0.0)
	{
	  ASSERT ((h5_file = H5Fcreate(h5_path.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)) >= 0, "HDF5 file open error");
	}
	else
	{
	  ASSERT((h5_file = H5Fopen(h5_path.str().c_str(), H5F_ACC_RDWR, H5P_DEFAULT)) >= 0, "HDF5 file open error");
	}

	h5_text.str("");
	h5_text.clear();
	h5_text << output_frame;
	ASSERT ((h5_step = H5Gcreate (h5_file, h5_text.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) >= 0, "HDF5 file write error");
	double time = curtime;
        ASSERT (H5LTset_attribute_double (h5_step, ".", "TIME", &time, 1) >= 0, "HDF5 file write error"); 
   
	h5_rb_dataset (num, pset, ent, h5_step); /* append h5 file */

        xmf_path.str("");
	xmf_path.clear();
	xmf_path << output_path << j+1 << "rb.xmf"; /* append xmf file */

	int elements = 0;
	int nodes = num;
	int topo_size = 0;
	const char *label = "PARMEC rigid bodies";
	string h5file = h5_path.str().substr(h5_path.str().find_last_of('/')+1);

        append_xmf_file (xmf_path.str().c_str(), OUT_MODE_RB, elements, nodes, topo_size, label, h5file.c_str(), ent);

	H5Gclose (h5_step);
	H5Fclose (h5_file);
      }
    }

    delete [] set;
  }

  /* TODO --> *cd.h5 amd *cd.xmf contact data output */

  if (sprnum)
  {
    int i, j, ent, num, *set; /* number of and the set of springs */
    MAP **map, *item; /* map of springs attached to particles, and iterator */
    MEM mem;

    MEM_Init (&mem, sizeof (MAP), 1024);
    ERRMEM (map = static_cast<MAP**>(MEM_CALLOC(parnum * sizeof(MAP*))));
    ERRMEM (set = new int[sprnum]);

    for (j = -1; j < outnum; j ++)
    {
      num = 0;

      if (j < 0 && (outrest[1] & OUT_MODE_SD)) /* output springs attached to unselected particles */
      {
	for (i = 0; i < sprnum; i ++)
	{
	  if (flags[sprpart[0][i]] && OUTREST || /* first or second particle is unselected */
	  (sprpart[1][i] >= 0 && (flags[sprpart[1][i]] & OUTREST)))
	  {
	    set[num ++] = i;
	  }

	  MAP_Insert (&mem, &map[sprpart[0][i]], (void*)(long)i, NULL, NULL); /* map springs to particles */
	  if (sprpart[1][i] >= 0) MAP_Insert (&mem, &map[sprpart[1][i]], (void*)(long)i, NULL, NULL);
	}

        ent = outrest[0];
      }
      else if (outmode[j] & OUT_MODE_SD) /* output springs attached to selected particles */
      {
	for (num = 0, i = outidx[j]; i < outidx[j+1]; i ++)
	{
	  for (item = MAP_First (map[outpart[i]]); item; item = MAP_Next(item))
	  {
	    set[num ++] = (int)(long)item->key;
	  }
	}

	ent = outent[j];
      }

      if (num) 
      {
	h5_path.str("");
	h5_path.clear();
	h5_path << output_path << j+1 << "sd.h5";

	if (curtime == 0.0)
	{
	  ASSERT ((h5_file = H5Fcreate(h5_path.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT)) >= 0, "HDF5 file open error");
	}
	else
	{
	  ASSERT((h5_file = H5Fopen(h5_path.str().c_str(), H5F_ACC_RDWR, H5P_DEFAULT)) >= 0, "HDF5 file open error");
	}

	h5_text.str("");
	h5_text.clear();
	h5_text << output_frame;
	ASSERT ((h5_step = H5Gcreate (h5_file, h5_text.str().c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT)) >= 0, "HDF5 file write error");
	double time = curtime;
        ASSERT (H5LTset_attribute_double (h5_step, ".", "TIME", &time, 1) >= 0, "HDF5 file write error"); 
 
	h5_spring_dataset (num, set, ent, h5_step); /* append h5 dataset */

        xmf_path.str("");
	xmf_path.clear();
	xmf_path << output_path << j+1 << "sd.xmf"; /* append xmf file */

	int elements = 0;
	int nodes = num;
	int topo_size = 0;
	const char *label = "PARMEC springs";
	string h5file = h5_path.str().substr(h5_path.str().find_last_of('/')+1);

        append_xmf_file (xmf_path.str().c_str(), OUT_MODE_SD, elements, nodes, topo_size, label, h5file.c_str(), ent);

	H5Gclose (h5_step);
	H5Fclose (h5_file);
      }
    }

    delete [] set;
    free (map);
    MEM_Release (&mem);
  }
};

/* output VTK files */
static void output_vtk_files ()
{
  ostringstream oss;
  ofstream out;
  int i, j;

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
	oss << output_path << j+1 << ".vtk." << output_frame;
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

	vtk_triangle_dataset (num, set, outrest[0], out);

	out.close();
      }
      else if (outmode[j] & OUT_MODE_MESH) /* output selected triangles */
      {
	oss.str("");
	oss.clear();
	oss << output_path << j+1 << ".vtk." << output_frame;
	out.open (oss.str().c_str());

	out << "# vtk DataFile Version 2.0\n";
	out << "PARMEC triangles output at time " << curtime << "\n";
	out << "ASCII\n";

	num = find_triangle_set (&outpart[outidx[j]], &outpart[outidx[j+1]], set);

	if (num) vtk_triangle_dataset (num, set, outent[j], out);

        out.close();
      }
    }

    delete [] set;
  }

  if (parnum)
  {
    int num, *set; /* number of and the set of unselected particles */

    ERRMEM (set = new int[parnum]);

    for (j = -1; j < outnum; j ++)
    {
      if (j < 0 && (outrest[1] & OUT_MODE_RB)) /* output unselected particles */
      {
	oss.str("");
	oss.clear();
	oss << output_path << j+1 << "rb.vtk." << output_frame;
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
	oss << output_path << j+1 << "rb.vtk." << output_frame;
	out.open (oss.str().c_str());

	out << "# vtk DataFile Version 2.0\n";
	out << "PARMEC rigid bodies output at time " << curtime << "\n";
	out << "ASCII\n";

	output_rb_dataset (outidx[j+1]-outidx[j], &outpart[outidx[j]], outent[j], out);

        out.close();
      }
    }

    delete [] set;
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
	oss << output_path << j+1 << "sd.vtk." << output_frame;
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
	oss << output_path << j+1 << "sd.vtk." << output_frame;
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

    delete [] set;
    free (map);
    MEM_Release (&mem);
  }
}

#if MED
/* output MED files; based on
 * https://gitlab.onelab.info/gmsh/gmsh/blob/master/Geo/GModelIO_MED.cpp */
static void output_med_files ()
{
  ostringstream oss;
  int i, j;

  if (trinum)
  {
    int num, *set; /* number of and the set of triangles */

    ERRMEM (set = new int[trinum]);

    for (j = -1; j < outnum; j ++)
    {
      if (j < 0 && (outrest[1] & OUT_MODE_MESH)) /* output unselected triangles */
      {
	if (med_fid < 0)
	{
	  oss.str("");
	  oss.clear();
	  oss << output_path << j+1 << ".med";
	  med_fid = MEDfileOpen((char*)oss.str().c_str(), MED_ACC_CREAT);
	  ASSERT (med_fid >= 0, "Unable to open file '%s'", oss.str().c_str());

	  oss.str("");
	  oss.clear();
	  oss << "PARMEC triangles output";
	  ASSERT(MEDfileCommentWr(med_fid, (char*)oss.str().c_str()) >= 0, "Unable to write MED descriptor");
	}

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

	med_triangle_dataset (num, set, outrest[0], med_fid);
      }
      else if (outmode[j] & OUT_MODE_MESH) /* output selected triangles */
      {
	if (med_fid < 0)
	{
	  oss.str("");
	  oss.clear();
	  oss << output_path << j+1 << ".med";
	  med_fid = MEDfileOpen((char*)oss.str().c_str(), MED_ACC_CREAT);
	  ASSERT (med_fid >= 0, "Unable to open file '%s'", oss.str().c_str());

	  oss.str("");
	  oss.clear();
	  oss << "PARMEC triangles output";
	  ASSERT(MEDfileCommentWr(med_fid, (char*)oss.str().c_str()) >= 0, "Unable to write MED descriptor");
	}

	num = find_triangle_set (&outpart[outidx[j]], &outpart[outidx[j+1]], set);

	if (num) med_triangle_dataset (num, set, outent[j], med_fid);
      }
    }

    delete [] set;
  }
}
#endif

/* output *?.dump files */
static void output_dump_files ()
{
  ostringstream oss;
  ofstream out;
  int i, j;

  if (ellnum)
  {
    int num, *set; /* number of and the set of unselected ellipsoids */

    ERRMEM (set = new int[ellnum]);

    for (j = -1; j < outnum; j ++)
    {
      if (j < 0 && (outrest[1] & OUT_MODE_SPH)) /* output unselected particles */
      {
	oss.str("");
	oss.clear();
	oss << output_path << j+1 << ".dump";
	if (output_frame) out.open (oss.str().c_str(), ios::app);
	else out.open (oss.str().c_str());

	for (num = i = 0; i < ellnum; i ++)
	{
	  if (flags[part[i]] & OUTREST)
	  {
	    set[num ++] = i;
	  }
	}

	output_dump_dataset (num, set, outrest[0], out);

        out.close();
      }
      else if (outmode[j] & OUT_MODE_SPH) /* output selected particles */
      {
	oss.str("");
	oss.clear();
	oss << output_path << j+1 << ".dump";
	if (output_frame) out.open (oss.str().c_str(), ios::app);
	else out.open (oss.str().c_str());

	num = find_ellipsoid_set (&outpart[outidx[j]], &outpart[outidx[j+1]], set);

	output_dump_dataset (num, set, outent[j], out);


        out.close();
      }
    }

    delete [] set;
  }
}

/* runtime file output */
void output_files ()
{
  if (outformat & OUT_FORMAT_DUMP) output_dump_files();

  if (outformat & OUT_FORMAT_VTK) output_vtk_files();

  if (outformat & OUT_FORMAT_XDMF) output_xdmf_files();

#if MED
  if (outformat & OUT_FORMAT_MED) output_med_files ();
#endif

  output_frame ++; 
}

/* runtime history output */
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
	  case HIS_LENGTH:
	  {
	    int l = sprmap[k];
	    REAL q[3] = {sprpnt[1][0][l]-sprpnt[0][0][l],
	                 sprpnt[1][1][l]-sprpnt[0][1][l],
	                 sprpnt[1][2][l]-sprpnt[0][2][l]};
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
	  case HIS_SS:
	  {
	    int l = sprmap[k];
	    value += (REAL)unspring[l];
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

/* read .h5 dataset */
static double* h5read (hid_t h5_step, const char *name, int *size)
{
  double *values = NULL;

  /* XXX
   * if (H5LTfind_dataset(h5_step, name)) is unreliable as
   * it passes for names being substrings of dataset names;
   * hence the below workaround */

  hsize_t cur, num;
  char buf[128];
  H5Gget_num_objs (h5_step, &num);
  for (cur = 0; cur < num; cur ++)
  {
    H5Gget_objname_by_idx (h5_step, cur, buf, 128);
    if (strcmp (buf, name) == 0) break;
  }

  if (cur < num)
  {
    int rank;
    H5LTget_dataset_ndims (h5_step, name, &rank);
    ASSERT (rank <= 2, "HDF5 file read error: rank > 2");
    hsize_t dims[2] = {0, 1};
    H5T_class_t class_id;
    size_t type_size;
    H5LTget_dataset_info (h5_step, name, dims, &class_id, &type_size);
    ASSERT (type_size == 8, "HDF5 file read error: expected double precision float");
    ERRMEM (values = (double*)malloc (type_size*dims[0]*dims[1]));
    ASSERT (H5LTread_dataset_double (h5_step, name, values) >= 0, "HDF5 file read error");
    if (size) *size = dims[0];
  }
  else if (size) *size = 0;

  return values;
}

/* output history from existing .h5 files */
void output_h5history ()
{
  MAP *h5map;
  int i;

  for (i = 0, h5map = NULL; i < hisnum; i ++)
  {
    if (h5file[i])
    {
      MAP *node = MAP_Find_Node (h5map, h5file[i], (MAP_Compare)strcmp);

      if (node) /* .h5 path already mapped */
      {
	MAP_Insert (NULL, (MAP**) &node->data, (void*) (long) i, NULL, NULL); /* append history index */
      }
      else  /* .h5 path not mapped */
      {
	MAP *h5set = NULL; /* empty index set */
	MAP_Insert (NULL, &h5set, (void*) (long) i, NULL, NULL); /* append history index */
	MAP_Insert (NULL, &h5map, h5file[i], h5set, (MAP_Compare)strcmp); /* map path to index set */
      }

      if (PyList_Size((PyObject*)history[i]) > 0)
      {
	Py_DECREF ((PyObject*)history[i]);
	history[i] = PyList_New(0); /* avoid appending to nonempty lists */
      }
    }
  }

  for (MAP *h5_path = MAP_First(h5map); h5_path; h5_path = MAP_Next (h5_path)) /* for each .h5 file */
  {
    MAP *ent2data = (MAP*)h5_path->data;
    hid_t h5_file, h5_step;
    hsize_t nsteps;
    double time;

    ASSERT((h5_file = H5Fopen((const char*)h5_path->key, H5F_ACC_RDONLY, H5P_DEFAULT)) >= 0, "HDF5 file open error");
    H5Gget_num_objs (h5_file, &nsteps);

    ASSERT ((h5_step = H5Gopen (h5_file, "/0", H5P_DEFAULT)) >= 0, "HDF5 file read error");
    int NUM0;
    double *GEOM0 = h5read (h5_step, "GEOM", &NUM0);
    H5Gclose (h5_step);

    for (int n = 0; n < nsteps; n ++)
    {
      char buf[1024];

      snprintf (buf, 1024, "/%d", n);
      ASSERT ((h5_step = H5Gopen (h5_file, buf, H5P_DEFAULT)) >= 0, "HDF5 file read error");

      int NUM = 0;
      double *ANGVEL = h5read (h5_step, "ANGVEL", &NUM);
      double *DISPL = h5read (h5_step, "DISPL", &NUM);
      double *FORCE = h5read (h5_step, "FORCE", &NUM);
      double *GEOM = h5read (h5_step, "GEOM", &NUM);
      double *LINVEL = h5read (h5_step, "LINVEL", &NUM);
      double *ORIENT = h5read (h5_step, "ORIENT", &NUM);
      double *TORQUE = h5read (h5_step, "TORQUE", &NUM);
      double *F = h5read (h5_step, "F", &NUM);
      double *LENGTH = h5read (h5_step, "LENGTH", &NUM);
      double *SF = h5read (h5_step, "SF", &NUM);
      double *SS = h5read (h5_step, "SS", &NUM);

      for (MAP *item = MAP_First (ent2data); item; item = MAP_Next (item))
      {
	i = (int) (long) item->key;
	
	if (hisent[i] == HIS_TIME)
	{
          ASSERT (H5LTget_attribute_double (h5_step, ".", "TIME", &time) >= 0, "HDF5 file read error");
          PyList_Append ((PyObject*)history[i], PyFloat_FromDouble(time));
	}
	else switch (hiskind[i]&(HIS_LIST|HIS_SPHERE|HIS_BOX))
	{
	case HIS_LIST:
	{
	  if (hiskind[i] & HIS_POINT) /* one particle point based */
	  {
	    int k = hislst[hisidx[i]];
	    REAL value;

	    switch(hisent[i])
	    {
	    case HIS_PX:
	    case HIS_PY:
	    case HIS_PZ:
	    case HIS_PL:
	    case HIS_DX:
	    case HIS_DY:
	    case HIS_DZ:
	    case HIS_DL:
	    {
	      ASSERT (GEOM0 && GEOM, "HDF5 file read error: GEOM dataset missing");
	      ASSERT (ORIENT, "HDF5 file read error: ORIENT dataset missing");
	      REAL x[3] = {GEOM[k*3], GEOM[k*3+1], GEOM[k*3+2]};
	      REAL X[3] = {GEOM0[k*3], GEOM0[k*3+1], GEOM0[k*3+2]};
	      REAL L[9] = {ORIENT[k*9], ORIENT[k*9+1], ORIENT[k*9+2],
			   ORIENT[k*9+3], ORIENT[k*9+4], ORIENT[k*9+5],
			   ORIENT[k*9+6], ORIENT[k*9+7], ORIENT[k*9+8]};
	      REAL P[3] = {source[0][i], source[1][i], source[2][i]};
	      REAL Q[3], p[3];

	      SUB (P, X, Q);
	      NVADDMUL (x, L, Q, p);

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
	      }
	    }
	    break;
	    case HIS_VX:
	    case HIS_VY:
	    case HIS_VZ:
	    case HIS_VL:
	    {
	      ASSERT (GEOM0 && GEOM, "HDF5 file read error: GEOM dataset missing");
	      ASSERT (ORIENT, "HDF5 file read error: ORIENT dataset missing");
	      ASSERT (LINVEL, "HDF5 file read error: LINVEL dataset missing");
	      ASSERT (ANGVEL, "HDF5 file read error: ANGVEL dataset missing");
	      REAL x[3] = {GEOM[k*3], GEOM[k*3+1], GEOM[k*3+2]};
	      REAL X[3] = {GEOM0[k*3], GEOM0[k*3+1], GEOM0[k*3+2]};
	      REAL v[3] = {LINVEL[k*3], LINVEL[k*3+1], LINVEL[k*3+2]};
	      REAL o[3] = {ANGVEL[k*3], ANGVEL[k*3+1], ANGVEL[k*3+2]};
	      REAL L[9] = {ORIENT[k*9], ORIENT[k*9+1], ORIENT[k*9+2],
			   ORIENT[k*9+3], ORIENT[k*9+4], ORIENT[k*9+5],
			   ORIENT[k*9+6], ORIENT[k*9+7], ORIENT[k*9+8]};
	      REAL P[3] = {source[0][i], source[1][i], source[2][i]};
	      REAL Q[3], p[3], a[3];

	      SUB (P, X, Q);
	      NVADDMUL (x, L, Q, p);
	      SUB (p, x, a);

	      switch (hisent[i])
	      {
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
	      }
	    }
	    break;
	    case HIS_OX:
	    case HIS_OY:
	    case HIS_OZ:
	    case HIS_OL:
	    {
	      ASSERT (ANGVEL, "HDF5 file read error: ANGVEL dataset missing");
	      REAL o[3] = {ANGVEL[k*3], ANGVEL[k*3+1], ANGVEL[k*3+2]};

	      switch (hisent[i])
	      {
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
	      }
	    }
	    break;
	    case HIS_FX:
	    case HIS_FY:
	    case HIS_FZ:
	    case HIS_FL:
	    {
	      ASSERT (FORCE, "HDF5 file read error: FORCE dataset missing");
	      REAL f[3] = {FORCE[k*3], FORCE[k*3+1], FORCE[k*3+2]};

	      switch (hisent[i])
	      {
	      case HIS_FX:
		value = f[0];
	      break;
	      case HIS_FY:
		value = f[1];
	      break;
	      case HIS_FZ:
		value = f[2];
	      break;
	      case HIS_FL:
	      {
		value = LEN(f);
	      }
	      break;
	      }
	    }
	    break;
	    case HIS_TX:
	    case HIS_TY:
	    case HIS_TZ:
	    case HIS_TL:
	    {
	      ASSERT (TORQUE, "HDF5 file read error: TORQUE dataset missing");
	      REAL t[3] = {TORQUE[k*3], TORQUE[k*3+1], TORQUE[k*3+2]};

	      switch (hisent[i])
	      {
	      case HIS_TX:
		value = t[0];
	      break;
	      case HIS_TY:
		value = t[1];
	      break;
	      case HIS_TZ:
		value = t[2];
	      break;
	      case HIS_TL:
	      {
		value = LEN(t);
	      }
	      break;
	      }
	    }
	    break;
	    }

	    PyList_Append ((PyObject*)history[i], PyFloat_FromDouble(value));
	  }
	  else /* particle or spring list based */
	  {
	    REAL value = 0.0;

	    for (int j = hisidx[i]; j < hisidx[i+1]; j ++)
	    {
	      int k = hislst[j];

	      switch (hisent[i])
	      {
	      case HIS_PX:
	        ASSERT (GEOM, "HDF5 file read error: GEOM dataset missing");
		value += GEOM[3*k];
	      break;
	      case HIS_PY:
	      {
	        ASSERT (GEOM, "HDF5 file read error: GEOM dataset missing");
		value += GEOM[3*k+1];
	      }
	      break;
	      case HIS_PZ:
	      {
	        ASSERT (GEOM, "HDF5 file read error: GEOM dataset missing");
		value += GEOM[3*k+2];
	      }
	      break;
	      case HIS_PL:
	      {
	        ASSERT (GEOM, "HDF5 file read error: GEOM dataset missing");
		REAL q[3] = {GEOM[3*k], GEOM[3*k+1], GEOM[3*k+2]};
		value += LEN(q);
	      }
	      break;
	      case HIS_DX:
	        ASSERT (GEOM0 && GEOM, "HDF5 file read error: GEOM dataset missing");
		value += GEOM[3*k] - GEOM0[3*k];
	      break;
	      case HIS_DY:
	        ASSERT (GEOM0 && GEOM, "HDF5 file read error: GEOM dataset missing");
		value += GEOM[3*k+1] - GEOM0[3*k+1];
	      break;
	      case HIS_DZ:
	        ASSERT (GEOM0 && GEOM, "HDF5 file read error: GEOM dataset missing");
		value += GEOM[3*k+2] - GEOM0[3*k+2];
	      break;
	      case HIS_DL:
	      {
	        ASSERT (GEOM0 && GEOM, "HDF5 file read error: GEOM dataset missing");
		REAL q[3] = {GEOM[3*k]-GEOM0[3*k], GEOM[3*k+1]-GEOM0[3*k+1], GEOM[3*k+2]-GEOM[3*k+2]};
		value += LEN(q);
	      }
	      break;
	      case HIS_VX:
	        ASSERT (LINVEL, "HDF5 file read error: LINVEL dataset missing");
		value += LINVEL[3*k];
	      break;
	      case HIS_VY:
	        ASSERT (LINVEL, "HDF5 file read error: LINVEL dataset missing");
		value += LINVEL[3*k+1];
	      break;
	      case HIS_VZ:
	        ASSERT (LINVEL, "HDF5 file read error: LINVEL dataset missing");
		value += LINVEL[3*k+2];
	      break;
	      case HIS_VL:
	      {
	        ASSERT (LINVEL, "HDF5 file read error: LINVEL dataset missing");
		REAL q[3] = {LINVEL[3*k], LINVEL[3*k+1], LINVEL[3*k+2]};
		value += LEN(q);
	      }
	      break;
	      case HIS_OX:
	        ASSERT (ANGVEL, "HDF5 file read error: ANGVEL dataset missing");
		value += ANGVEL[3*k];
	      break;
	      case HIS_OY:
	        ASSERT (ANGVEL, "HDF5 file read error: ANGVEL dataset missing");
		value += ANGVEL[3*k+1];
	      break;
	      case HIS_OZ:
	        ASSERT (ANGVEL, "HDF5 file read error: ANGVEL dataset missing");
		value += ANGVEL[3*k+2];
	      break;
	      case HIS_OL:
	      {
	        ASSERT (ANGVEL, "HDF5 file read error: ANGVEL dataset missing");
		REAL q[3] = {ANGVEL[3*k], ANGVEL[3*k+1], ANGVEL[3*k+2]};
		value += LEN(q);
	      }
	      break;
	      case HIS_FX:
	        ASSERT (FORCE, "HDF5 file read error: FORCE dataset missing");
		value += FORCE[3*k];
	      break;
	      case HIS_FY:
	        ASSERT (FORCE, "HDF5 file read error: FORCE dataset missing");
		value += FORCE[3*k+1];
	      break;
	      case HIS_FZ:
	        ASSERT (FORCE, "HDF5 file read error: FORCE dataset missing");
		value += FORCE[3*k+2];
	      break;
	      case HIS_FL:
	      {
	        ASSERT (FORCE, "HDF5 file read error: FORCE dataset missing");
		REAL q[3] = {FORCE[3*k], FORCE[3*k+1], FORCE[3*k+2]};
		value += LEN(q);
	      }
	      break;
	      case HIS_TX:
	        ASSERT (TORQUE, "HDF5 file read error: TORQUE dataset missing");
		value += TORQUE[3*k];
	      break;
	      case HIS_TY:
	        ASSERT (TORQUE, "HDF5 file read error: TORQUE dataset missing");
		value += TORQUE[3*k+1];
	      break;
	      case HIS_TZ:
	        ASSERT (TORQUE, "HDF5 file read error: TORQUE dataset missing");
		value += TORQUE[3*k+2];
	      break;
	      case HIS_TL:
	      {
	        ASSERT (TORQUE, "HDF5 file read error: TORQUE dataset missing");
		REAL q[3] = {TORQUE[3*k], TORQUE[3*k+1], TORQUE[3*k+2]};
		value += LEN(q);
	      }
	      break;
	      case HIS_LENGTH:
	      {
		int l = sprmap[k];
	        ASSERT (LENGTH, "HDF5 file read error: LENGTH dataset missing");
		value += LENGTH[l];
	      }
	      break;
	      case HIS_STROKE:
	      {
		int l = sprmap[k];
	        ASSERT (DISPL, "HDF5 file read error: DISPL dataset missing");
		value += DISPL[l];
	      }
	      break;
	      case HIS_STF:
	      {
		int l = sprmap[k];
	        ASSERT (F, "HDF5 file read error: F dataset missing");
		value += F[l];
	      }
	      break;
	      case HIS_SF:
	      {
		int l = sprmap[k];
	        ASSERT (SF, "HDF5 file read error: SF dataset missing");
		value += SF[l];
	      }
	      break;
	      case HIS_SS:
	      {
		int l = sprmap[k];
	        ASSERT (SS, "HDF5 file read error: SS dataset missing");
		value += SS[l];
	      }
	      break;
	      }
	    }

	    int j = hisidx[i+1]-hisidx[i];

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

      free (ANGVEL);
      free (DISPL);
      free (FORCE);
      free (GEOM);
      free (LINVEL);
      free (ORIENT);
      free (TORQUE);
      free (F);
      free (LENGTH);
      free (SF);
      free (SS);
      H5Gclose (h5_step);
    }

    free (GEOM0);
    H5Fclose (h5_file);
  }
}

/* close files and reset global output variables */
void output_reset ()
{
#if MED
  if (med_fid >= 0)
  {
    ASSERT (MEDfileClose (med_fid) >= 0, "Closing MED file has failed");
  }

  med_fid = (med_idt)-1;
#endif
}

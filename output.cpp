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

namespace parmec
{
int output_frame = 0; /* output files frame */
}

using namespace parmec;

/* output files */
void output_files ()
{
  using namespace std;
  ostringstream oss;
  ofstream out;
  int i;

  if (ellnum)
  {
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
  }

  if (trinum)
  {
    oss << outpath << output_frame << ".vtk";
    out.open (oss.str().c_str());

    out << "# vtk DataFile Version 2.0\n";
    out << "PARMEC triangles output\n";
    out << "ASCII\n";
    out << "DATASET UNSTRUCTURED_GRID\n";
    out << "POINTS " << 3*trinum << " float\n";
    for (i = 0; i < trinum; i ++)
    {
      out << tri [0][0][i] << " " << tri[0][1][i] << " " << tri[0][2][i] << " "
          << tri [1][0][i] << " " << tri[1][1][i] << " " << tri[1][2][i] << " "
          << tri [2][0][i] << " " << tri[2][1][i] << " " << tri[2][2][i] << "\n";
    }
    out << "CELLS " << trinum << " " << 4*trinum << "\n";
    for (i = 0; i < trinum; i ++)
    {
      out << 3 << " " << 3*i << " " << 3*i+1 << " " << 3*i+2 << "\n";
    }
    out << "CELL_TYPES " << trinum << "\n";
    for (i = 0; i < trinum; i ++)
    {
      out << 5 << "\n";
    }
    out << "CELL_DATA " << trinum << "\n";
    out << "SCALARS colors int\n";
    out << "LOOKUP_TABLE default\n";
    for (i = 0; i < trinum; i ++)
    {
      out << tricol[i] << "\n";
    }
 
    out.close();
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
	k = hispart[hisidx[i]];

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
      else /* particle list based */
      {
	REAL value = 0.0;

	for (j = hisidx[i]; j < hisidx[i+1]; j ++)
	{
	  k = hispart[j];

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

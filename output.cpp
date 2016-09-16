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
int output_frame = 0; /* output frame */
}

using namespace parmec;

/* output current state */
void output ()
{
  using namespace std;
  ostringstream oss;
  ofstream out;
  int i, j, k;

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

  for (i = 0; i < hisnum; i ++) /* append time histories */
  {
    switch (hiskind[i]&(HIS_LIST|HIS_SPHERE|HIS_BOX))
    {
    case HIS_LIST:
    {
      if (hiskind[i] & HIS_POINT) /* one particle point based */
      {
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
	  break;
	  case HIS_PY:
	  break;
	  case HIS_PZ:
	  break;
	  case HIS_PL:
	  break;
	  case HIS_TIME:
	  break;
	  /* TODO */
	  }
	}
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

  output_frame ++; 
}

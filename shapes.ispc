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

/* ellipsoids task */
task void ellipsoids_task (uniform int span, uniform int ellnum, uniform int part[],
  uniform REAL * uniform center[6], uniform REAL * uniform radii[3], uniform REAL * uniform orient[18],
  uniform REAL * uniform rotation[9], uniform REAL * uniform position[6])
{
  uniform int start = taskIndex*span;
  uniform int end = taskIndex == taskCount-1 ? ellnum : start+span;

  foreach (i = start ... end)
  {
    REAL L[9], X[3], x[3], C[3], c[3];

    int j = part[i];

    L[0] = rotation[0][j];
    L[1] = rotation[1][j];
    L[2] = rotation[2][j];
    L[3] = rotation[3][j];
    L[4] = rotation[4][j];
    L[5] = rotation[5][j];
    L[6] = rotation[6][j];
    L[7] = rotation[7][j];
    L[8] = rotation[8][j];

    X[0] = position[3][j];
    X[1] = position[4][j];
    X[2] = position[5][j];

    x[0] = position[0][j];
    x[1] = position[1][j];
    x[2] = position[2][j];

    C[0] = center[3][i];
    C[1] = center[4][i];
    C[2] = center[5][i];

    SCC (X, C);
    NVADDMUL (x, L, C, c);

    center[0][i] = c[0];
    center[1][i] = c[1];
    center[2][i] = c[2];

    if (radii[1][j] > 0.) /* ellipsoid */
    {
      REAL O[9], o[9];

      O[0] = orient[9][i];
      O[1] = orient[10][i];
      O[2] = orient[11][i];
      O[3] = orient[12][i];
      O[4] = orient[13][i];
      O[5] = orient[14][i];
      O[6] = orient[15][i];
      O[7] = orient[16][i];
      O[8] = orient[17][i];

      NNMUL (L, O, o);

      orient[0][i] = o[0];
      orient[1][i] = o[1];
      orient[2][i] = o[2];
      orient[3][i] = o[3];
      orient[4][i] = o[4];
      orient[5][i] = o[5];
      orient[6][i] = o[6];
      orient[7][i] = o[7];
      orient[8][i] = o[8];
    }
  }
}

/* nodes task */
task void nodes_task (uniform int span, uniform int nodnum, uniform REAL * uniform nodes [6], uniform int nodpart[],
  uniform int * uniform analytical, uniform REAL * uniform rotation[9], uniform REAL * uniform position[6])
{
  uniform int start = taskIndex*span;
  uniform int end = taskIndex == taskCount-1 ? nodnum : start+span;

  if (analytical)
  {
    foreach (i = start ... end)
    {
      REAL L[9], X[3], x[3], C[3], c[3];

      cif (analytical[i] == 0) /* non-analytical update */
      {
	int j = nodpart[i];

	L[0] = rotation[0][j];
	L[1] = rotation[1][j];
	L[2] = rotation[2][j];
	L[3] = rotation[3][j];
	L[4] = rotation[4][j];
	L[5] = rotation[5][j];
	L[6] = rotation[6][j];
	L[7] = rotation[7][j];
	L[8] = rotation[8][j];

	X[0] = position[3][j];
	X[1] = position[4][j];
	X[2] = position[5][j];

	x[0] = position[0][j];
	x[1] = position[1][j];
	x[2] = position[2][j];

	C[0] = nodes[3][i];
	C[1] = nodes[4][i];
	C[2] = nodes[5][i];

	SCC (X, C);
	NVADDMUL (x, L, C, c);

	nodes[0][i] = c[0];
	nodes[1][i] = c[1];
	nodes[2][i] = c[2];
      }
    }
  }
  else /* full update */
  {
    foreach (i = start ... end)
    {
      REAL L[9], X[3], x[3], C[3], c[3];

      int j = nodpart[i];

      L[0] = rotation[0][j];
      L[1] = rotation[1][j];
      L[2] = rotation[2][j];
      L[3] = rotation[3][j];
      L[4] = rotation[4][j];
      L[5] = rotation[5][j];
      L[6] = rotation[6][j];
      L[7] = rotation[7][j];
      L[8] = rotation[8][j];

      X[0] = position[3][j];
      X[1] = position[4][j];
      X[2] = position[5][j];

      x[0] = position[0][j];
      x[1] = position[1][j];
      x[2] = position[2][j];

      C[0] = nodes[3][i];
      C[1] = nodes[4][i];
      C[2] = nodes[5][i];

      SCC (X, C);
      NVADDMUL (x, L, C, c);

      nodes[0][i] = c[0];
      nodes[1][i] = c[1];
      nodes[2][i] = c[2];
    }
  }
}

/* faces task */
task void faces_task (uniform int span, uniform int facnum, uniform int * uniform facnod[3],
  uniform int factri[], uniform REAL * uniform nodes [6], uniform REAL * uniform tri[3][3])
{
  uniform int start = taskIndex*span;
  uniform int end = taskIndex == taskCount-1 ? facnum : start+span;

  foreach (i = start ... end)
  {
    int j = factri[i];

    int k[3] = {facnod[0][i], facnod[1][i], facnod[2][i]};

    tri[0][0][j] = nodes[0][k[0]];
    tri[0][1][j] = nodes[1][k[0]];
    tri[0][2][j] = nodes[2][k[0]];
    tri[1][0][j] = nodes[0][k[1]];
    tri[1][1][j] = nodes[1][k[1]];
    tri[1][2][j] = nodes[2][k[1]];
    tri[2][0][j] = nodes[0][k[2]];
    tri[2][1][j] = nodes[1][k[2]];
    tri[2][2][j] = nodes[2][k[2]];
  }
}

/* update shapes */
export void shapes (uniform int threads, uniform int ellnum, uniform int part[],
  uniform REAL * uniform center[6], uniform REAL * uniform radii[3], uniform REAL * uniform orient[18],
  uniform int nodnum, uniform REAL * uniform nodes[6], uniform int nodpart[], uniform int * uniform analytical,
  uniform int facnum, uniform int * uniform facnod[3], uniform int factri[], uniform REAL * uniform tri[3][3],
  uniform REAL * uniform rotation[9], uniform REAL * uniform position[6])
{
  launch [threads] ellipsoids_task (ellnum/threads, ellnum, part, center, radii, orient, rotation, position);

  launch [threads] nodes_task (nodnum/threads, nodnum, nodes, nodpart, analytical, rotation, position);
  sync;

  launch [threads] faces_task (facnum/threads, facnum, facnod, factri, nodes, tri);
}

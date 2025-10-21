/*
   The MIT License (MIT)

   Copyright (c) 2016 Tomasz Koziara

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

#include "mem.h"

#ifndef __mesh__
#define __mesh__

typedef struct element ELEMENT;
typedef struct face FACE;
typedef struct mesh_data MESH_DATA;

/* triangular or quadrilateral face */
struct face
{
  REAL normal [3]; /* current normal */

  int type, /* 3, 4 => triangle, quadrilateral */
      nodes [4], /* node numbers */
      index, /* index within the element */
      color; /* surface color */

  ELEMENT *ele;

  FACE *next, *n; /* element, mesh list */
};

/* finite element */
struct element
{
  short type, /* 4, 5, 6, 8 => tetrahedron, pyramid, wedge, hexahedron */
        neighs, /* number of neighbours */
        flag;  /* auxiliary flag used internally */

  int nodes [8], /* node numbers */
      material; /* material number */

  ELEMENT *prev,
          *next,
          *adj [6]; /* neighbouring elements */

  FACE *faces; /* corresponding surface faces */
};

/* general mesh */
struct mesh_data
{
  MEM facmem,
      elemem,
      mapmem;

  REAL (*nodes) [3];

  ELEMENT *surfeles,
          *bulkeles;

  FACE *faces;

  int  surfeles_count,
       bulkeles_count,
       nodes_count;
};

/* create mesh from a vector of nodes, element list in format =>
 * {nuber of nodes, node0, node1, ..., material}, {REPEAT}, ..., 0 (end of list); and surface colors in format =>
 * global surface, {number of nodes, node0, node1, ..., color}, {REPEAT}, ..., 0 (end of list); */
MESH_DATA* MESH_Create (REAL (*nodes) [3], int *elements, int *surfaces);

/* calculate mass characteristics: scalar mass, mass center, inertia tensor */
void MESH_Char (MESH_DATA *msh, REAL *mass, REAL *center, REAL *inertia);

/* free mesh memory */
void MESH_Destroy (MESH_DATA *msh);

#endif

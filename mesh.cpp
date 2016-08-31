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

#include <limits.h>
#include <float.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include "macros.h"
#include "parmec.h"
#include "mem.h"
#include "map.h"
#include "mesh.h"
#include "simplex.h"

/* used in some pools */
#define MEMCHUNK 128

static int tet [][4] =  /* 1-based indexing as the element lists start with the number of nodes */
  {{1,3,2,0},
   {1,2,4,0},
   {2,3,4,0},
   {3,1,4,0}}, pyr [][4] = 
  {{1,4,3,2},
   {1,2,5,0},
   {2,3,5,0},
   {3,4,5,0},
   {4,1,5,0}}, wed [][4] =
  {{1,3,2,0},
   {4,5,6,0},
   {1,2,5,4},
   {2,3,6,5},
   {3,1,4,6}}, hex [][4] =
  {{1,4,3,2},
   {1,2,6,5},
   {2,3,7,6},
   {3,4,8,7},
   {1,5,8,4},
   {5,6,7,8}};

/* maximal number of neighbours */
inline static int neighs (int type)
{
  switch (type)
  {
    case 4: return 4;
    case 5: return 5;
    case 6: return 5;
    case 8: return 6;
  }
  return 0;
}

inline static void swap (int *a, int *b)
{
  int c = *a; *a = *b; *b = c;
}

/* quick sort of ints */
static void sort (int* begin, int* end)
{
  int *lower = begin, *upper = end,
    bound = *(begin+(end-begin)/2);
  
  while (lower <= upper)
  {
    while (*lower < bound) lower++;
    while (bound < *upper) upper--;

    if (lower < upper) swap (lower++, upper--);
    else lower++;
  }

  if (begin < upper) sort (begin, upper);
  if (upper < end) sort (upper+1, end);
}

static int lexcmp (int *a, int *b, int m)
{
  int n;
  for (n = 0; n < m; n ++)
  {
    if (a [n] < b [n]) return -1;
    else if (a [n] > b [n]) return 1;
  }
  return 0;
}

/* comparison used in face mapping */
static int face_compare (void *pone, void *ptwo)
{
  FACE *one = static_cast<FACE*>(pone), *two = static_cast<FACE*>(ptwo);
  if (one->type < two->type) return -1;
  else if (one->type > two->type) return 1;
  else return lexcmp (one->nodes, two->nodes, one->type);
}

static ELEMENT* create_element (MEM *elemem, int *element)
{
  ELEMENT *ele;
  int n;

  ele = static_cast<ELEMENT*>(MEM_Alloc (elemem));
  ele->type = element [0];
  for (n = 1; n <= element [0]; n ++)
    ele->nodes [n-1] = element [n];
  ele->material = element [n]; /* material number */

  return ele;
}

static void setup_face (ELEMENT *ele, int n, FACE *fac, int dosort)
{
  switch (ele->type)
  {
    case 4:
      fac->type = 3;
      fac->nodes [0] = ele->nodes [tet [n][0]-1]; /* shift due to the 1-based indexing */
      fac->nodes [1] = ele->nodes [tet [n][1]-1];
      fac->nodes [2] = ele->nodes [tet [n][2]-1];
      if (dosort) sort (fac->nodes, fac->nodes+2);
    break;
    case 5:
    if (n == 0)
    { fac->type = 4;
      fac->nodes [0] = ele->nodes [pyr [n][0]-1];
      fac->nodes [1] = ele->nodes [pyr [n][1]-1];
      fac->nodes [2] = ele->nodes [pyr [n][2]-1];
      fac->nodes [3] = ele->nodes [pyr [n][3]-1];
      if (dosort) sort (fac->nodes, fac->nodes+3); }
    else
    { fac->type = 3;
      fac->nodes [0] = ele->nodes [pyr [n][0]-1];
      fac->nodes [1] = ele->nodes [pyr [n][1]-1];
      fac->nodes [2] = ele->nodes [pyr [n][2]-1];
      if (dosort) sort (fac->nodes, fac->nodes+2); }
    break;
    case 6:
    if (n < 2)
    { fac->type = 3;
      fac->nodes [0] = ele->nodes [wed [n][0]-1];
      fac->nodes [1] = ele->nodes [wed [n][1]-1];
      fac->nodes [2] = ele->nodes [wed [n][2]-1];
      if (dosort) sort (fac->nodes, fac->nodes+2); }
    else
    { fac->type = 4;
      fac->nodes [0] = ele->nodes [wed [n][0]-1];
      fac->nodes [1] = ele->nodes [wed [n][1]-1];
      fac->nodes [2] = ele->nodes [wed [n][2]-1];
      fac->nodes [3] = ele->nodes [wed [n][3]-1];
      if (dosort) sort (fac->nodes, fac->nodes+3); }
    break;
    case 8:
      fac->type = 4;
      fac->nodes [0] = ele->nodes [hex [n][0]-1];
      fac->nodes [1] = ele->nodes [hex [n][1]-1];
      fac->nodes [2] = ele->nodes [hex [n][2]-1];
      fac->nodes [3] = ele->nodes [hex [n][3]-1];
      if (dosort) sort (fac->nodes, fac->nodes+3);
    break;
  }
}

static void setup_normal (REAL (*nodes) [3], FACE *fac)
{
  int n0, n1, n2;
  REAL *normal;

  n0 = fac->nodes [0];
  n1 = fac->nodes [1];
  n2 = fac->nodes [2];
  normal = fac->normal;
  NORMAL (nodes [n0], nodes [n1], nodes [n2], normal);
  NORMALIZE (normal);
}

/* get new element, old face list and create the element faces => return the new face list */
static FACE* create_faces (MEM *facmem, MEM *mapmem, MAP **faces, ELEMENT *ele, FACE *list)
{
  FACE *fac, tmp;
  int n, m;

  m = neighs (ele->type); 
      
  for (n = 0; n < m; n ++)
  {
    /* set up temporary face for to
     * be used as the map search key */
    setup_face (ele, n, &tmp, 1); /* nodes sorted for map key comparisons */
    fac = static_cast<FACE*>(MAP_Find (*faces, &tmp, face_compare));  /* is it there ? */

    if (fac) /* was mapped already */
    {
      /* set up element adjacency */
      ele->adj [ele->neighs] = fac->ele;
      fac->ele->adj [fac->ele->neighs] = ele;
      fac->ele->neighs ++;
      ele->neighs ++;

      fac->ele = NULL; /* mark as the inner face (***) */
    }
    else
    {
      fac = static_cast<FACE*>(MEM_Alloc (facmem));
      fac->ele = ele;
      setup_face (ele, n, fac, 1);
      fac->index = n; /* local index */
      MAP_Insert (mapmem, faces, fac, /* map by the type/nodes key */
	fac, face_compare);
      fac->next = list;
      list = fac;
    }
  }

  return list;
}

static int maximal (int *element)
{
  int n, ret = 0;
  for (n = 1; n <= element [0]; n++)
    if (element [n] > ret) ret = element [n];
  return ret;
}

static int minimal (int *element)
{
  int n, ret = INT_MAX;
  for (n = 1; n <= element [0]; n++)
    if (element [n] < ret) ret = element [n];
  return ret;
}

static void element_char_add (MESH_DATA *msh, ELEMENT *ele, REAL *me, REAL *sx, REAL *sy, REAL *sz, REAL *euler)
{
  REAL rho = parmec::mparam[DENSITY][ele->material];
  REAL zero [3] = {0, 0, 0}, J, *a, *b, *c;
  int (*ver) [4], nv[8], i, j;

  switch (ele->type)
  {
  case 4:
    ver = tet;
    nv[0] = nv[1] = nv[2] = nv[3] = 4;
  break;
  case 5:
    ver = pyr;
    nv[0] = 4; nv[1] = nv[2] = nv[3] = nv[4] = 3;
  break;
  case 6:
    ver = wed;
    nv[0] = nv[1] = 3; nv[2] = nv[3] = nv[4] = 4;
  break;
  case 8:
    ver = hex;
    nv[0] = nv[1] = nv[2] = nv[3] = 
    nv[4] = nv[5] = nv[6] = nv[7] = 4;
  break;
  }

  for (i = 0; i < ele->type; i++)
  {
    a = msh->nodes[ele->nodes[ver[i][0]]];

    for (j = 1; j < nv[i]-1; j ++)
    {
      b = msh->nodes[ele->nodes[ver[i][j]]];
      c = msh->nodes[ele->nodes[ver[i][j+1]]];

      J = rho * simplex_J (zero, a, b, c);
      *me += simplex_1 (J, zero, a, b, c);
      *sx += simplex_x (J, zero, a, b, c);
      *sy += simplex_y (J, zero, a, b, c);
      *sz += simplex_z (J, zero, a, b, c);
      euler [0] += simplex_xx (J, zero, a, b, c);
      euler [3] += simplex_xy (J, zero, a, b, c);
      euler [4] += simplex_yy (J, zero, a, b, c);
      euler [6] += simplex_xz (J, zero, a, b, c);
      euler [7] += simplex_yz (J, zero, a, b, c);
      euler [8] += simplex_zz (J, zero, a, b, c);
    }
  }
}

/* create mesh from a vector of nodes, element list in format =>
 * {nuber of nodes, node0, node1, ..., material}, {REPEAT}, ..., 0 (end of list); and surface colors in format =>
 * global surface, {number of nodes, node0, node1, ..., surface}, {REPEAT}, ..., 0 (end of list); */
MESH_DATA* MESH_Create (REAL (*nodes) [3], int *elements, int *surfaces)
{
  int maximal_node,
      minimal_node,
      elements_count,
      faces_count,
      temp, *eleptr, n;
  REAL (*node) [3];
  MEM *elemem,
      facmem,
      mapmem;
  ELEMENT *ele, *enx, *elist;
  FACE *fac, *cac, *gac, *flist;
  MAP *faces, *smap;
  MESH_DATA *msh;
  
  maximal_node = 0;
  minimal_node = INT_MAX;
  elements_count = 0;
  faces_count = 0;

  /* create mesh storage */
  ERRMEM (msh = static_cast<MESH_DATA*>(MEM_CALLOC (sizeof (MESH_DATA))));
  elemem = &msh->elemem;
 
  /* calculate elements */ 
  for (eleptr = elements; eleptr [0]; eleptr += (eleptr [0]+2)) elements_count ++;

  MEM_Init (elemem, sizeof (ELEMENT), elements_count);
  MEM_Init (&facmem, sizeof (FACE), MEMCHUNK);
  MEM_Init (&mapmem, sizeof (MAP), MEMCHUNK);
  MEM_Init (&msh->mapmem, sizeof (MAP), MIN (elements_count, MEMCHUNK));

  elist = NULL;
  flist = NULL;
  faces = NULL;

  /* create elements list & face adjacency map */
  for (eleptr = elements; eleptr [0]; eleptr += (eleptr [0]+2))
  {
    ASSERT (
      eleptr [0] == 4 || /* tetrahedron */
      eleptr [0] == 5 || /* pyramid */
      eleptr [0] == 6 || /* wedge */
      eleptr [0] == 8,   /* hexahedron */
      "ERROR: unsupported element type");

    ele = create_element (elemem, eleptr);
    flist = create_faces (&facmem, &mapmem, &faces, ele, flist);
    ele->next = elist;
    elist = ele;

    /* node number extrema */
    temp = maximal (eleptr);
    if (temp > maximal_node)
      maximal_node = temp;
    temp = minimal (eleptr);
    if (temp < minimal_node)
      minimal_node = temp;
  }

  /* calculate faces */
  for (fac = flist; fac; fac = fac->next)
    if (fac->ele) faces_count ++;

  /* alocate additional storage */
  MEM_Init (&msh->facmem, sizeof (FACE), faces_count);
  msh->nodes_count = (maximal_node - minimal_node + 1);
  ERRMEM (msh->nodes = static_cast<REAL(*)[3]>(malloc (sizeof (REAL [3]) * (msh->nodes_count))));
  msh->surfeles_count = msh->bulkeles_count = 0;
  msh->surfeles = msh->bulkeles = NULL;

  /* set up elements */
  for (ele = elist; ele; ele = enx)
  {
    enx = ele->next;

    if (minimal_node > 0) /* impose 0-based indexing */
    {
      for (temp = 0; temp < ele->type; temp ++)
	ele->nodes [temp] -= minimal_node;
    }

    ele->prev = NULL;
   
    if (ele->neighs < neighs (ele->type)) /* surface element */
    {
      msh->surfeles_count ++;
      ele->next = msh->surfeles;
      if (msh->surfeles) msh->surfeles->prev = ele;
      msh->surfeles = ele;
    }
    else /* bulk element */
    {
      msh->bulkeles_count ++;
      ele->next = msh->bulkeles;
      if (msh->bulkeles) msh->bulkeles->prev = ele;
      msh->bulkeles = ele;
    }
  }

  /* create surfaces map => skip first element of 'surfaces' == the global surface kind */
  for (eleptr = (surfaces + 1), smap = NULL, temp = 0;
    eleptr [0]; eleptr += (eleptr [0]+2), temp ++)
  {
    fac = static_cast<FACE*>(MEM_Alloc (&facmem));
    
    ASSERT (
      eleptr [0] == 3 || /* triangle */
      eleptr [0] == 4,   /* quad */
      "ERROR: unsupported face type");

    fac->type = eleptr [0];
    for (n = 0; n < eleptr [0]; n ++)
      fac->nodes [n] = eleptr [n+1];
    sort (fac->nodes, fac->nodes+fac->type-1);

    fac->color = eleptr [eleptr [0] + 1];
    MAP_Insert (&mapmem, &smap, fac, /* map by the type/nodes key */
      fac, face_compare);
  }

  /* set up nodes */
  for (temp = minimal_node,
       node = msh->nodes;
       temp <= maximal_node;
       temp ++, node ++)
  {
    COPY (nodes [temp], *node);
  }

  /* set up faces */
  for (fac = flist; fac; fac = fac->next)
  {
    if (fac->ele) /* see (***) */
    {
      ele = fac->ele;

      cac = static_cast<FACE*>(MEM_Alloc (&msh->facmem));
      setup_face (ele, fac->index, cac, 0); /* setup face nodes without sorting them */
      cac->index = fac->index;
      cac->ele = fac->ele;
      setup_normal (msh->nodes, cac); /* calculate outer spatial normal */
      cac->next = ele->faces; /* append element face list */
      ele->faces = cac;

      /* set the mapped surface kind if possible => otherwise the global one */
      gac = static_cast<FACE*>(MAP_Find (smap, fac, face_compare));
      cac->color = (gac ? gac->color : surfaces [0]);
    }
  }

  /* create mesh face list */
  for (ele = msh->surfeles; ele; ele = ele->next)
  {
    for (fac = ele->faces; fac; fac = fac->next)
    {
      fac->n = msh->faces;
      msh->faces = fac;
    }
  }

  /* clean up */
  MEM_Release (&facmem);
  MEM_Release (&mapmem);

  return msh;
}

/* calculate mass characteristics: scalar mass, mass center, inertia tensor */
void MESH_Char (MESH_DATA *msh, REAL *mass, REAL *center, REAL *inertia)
{
  REAL me, sx, sy, sz, euler [9];
  ELEMENT *ele;

  me = sx = sy = sz = 0.0;
  SET9 (euler, 0.0);

  for (ele = msh->bulkeles; ele; ele = ele->next)
    element_char_add (msh, ele, &me, &sx, &sy, &sz, euler);

  for (ele = msh->surfeles; ele; ele = ele->next)
    element_char_add (msh, ele, &me, &sx, &sy, &sz, euler);

  center [0] = sx / me;
  center [1] = sy / me;
  center [2] = sz / me;

  euler [0] -= (2*sx - center[0]*me)*center[0];
  euler [4] -= (2*sy - center[1]*me)*center[1];
  euler [8] -= (2*sz - center[2]*me)*center[2];
  euler [3] -= center[0]*sy + center[1]*sx - center[0]*center[1]*me;
  euler [6] -= center[0]*sz + center[2]*sx - center[0]*center[2]*me;
  euler [7] -= center[1]*sz + center[2]*sy - center[1]*center[2]*me;
  euler [1] = euler[3];
  euler [2] = euler[6];
  euler [5] = euler[7];

  /* convert Euler tensor to the inertia tensor */
  REAL trace = TRACE (euler);
  IDENTITY (inertia);
  SCALE9 (inertia, trace);
  NNSUB (inertia, euler, inertia); /* inertia = tr(euler)*one - euler */
}

/* free mesh memory */
void MESH_Destroy (MESH_DATA *msh)
{
  MEM_Release (&msh->facmem);
  MEM_Release (&msh->elemem);
  MEM_Release (&msh->mapmem);
  free (msh->nodes);
  free (msh);
}

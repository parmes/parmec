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
#include <float.h>
#include <algorithm>
#include <vector>
#include <set>
#include <map>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <omp.h>
#include "macros.h"
#include "parmec.h"
#include "timer.h"
#include "mesh.h"
#include "parmec_ispc.h"
#include "partition_ispc.h"
#include "forces_ispc.h"
#include "dynamics_ispc.h"
#include "shapes_ispc.h"
#include "obstacles_ispc.h"

using namespace parmec;

#ifndef Py_RETURN_FALSE
#define Py_RETURN_FALSE return Py_INCREF(Py_False), Py_False
#endif

#ifndef Py_RETURN_TRUE
#define Py_RETURN_TRUE return Py_INCREF(Py_True), Py_True
#endif

#ifndef Py_RETURN_NONE
#define Py_RETURN_NONE return Py_INCREF(Py_None), Py_None
#endif

/* string buffer length */
#define BUFLEN 512

/* minimal type initialization */
#define TYPEINIT(typedesc, type, name, flags, dealloc, new, methods, members, getset)\
memset (&(typedesc), 0, sizeof (PyTypeObject));\
(typedesc).tp_basicsize = sizeof (type);\
(typedesc).tp_name = name;\
(typedesc).tp_flags = flags;\
(typedesc).tp_dealloc = (destructor)dealloc;\
(typedesc).tp_new = new;\
(typedesc).tp_methods = methods;\
(typedesc).tp_members = members;\
(typedesc).tp_getset = getset

/* string test */
static int is_string (PyObject *obj, const char *var)
{
  if (obj)
  {
    if (!PyString_Check (obj))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be a string", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }
  }

  return 1;
}

/* positive test */
static int is_positive (double num, const char *var)
{
  if (num <= 0)
  {
    char buf [BUFLEN];
    sprintf (buf, "'%s' must be positive", var);
    PyErr_SetString (PyExc_ValueError, buf);
    return 0;
  }

  return 1;
}

/* non negative test */
static int is_non_negative (double num, const char *var)
{
  if (num < 0)
  {
    char buf [BUFLEN];
    sprintf (buf, "'%s' must be non negative", var);
    PyErr_SetString (PyExc_ValueError, buf);
    return 0;
  }

  return 1;
}

/* in [lo, hi) range test */
static int is_in_range (double num, double lo, double hi, const char *var)
{
  if (num < lo || num >= hi)
  {
    char buf [BUFLEN];
    sprintf (buf, "'%s' must belong to [%g, %g)", var, lo, hi);
    PyErr_SetString (PyExc_ValueError, buf);
    return 0;
  }

  return 1;
}

/* tuple test */
static int is_tuple (PyObject *obj, const char *var, int len)
{
  if (obj)
  {
    if (!PyTuple_Check (obj))
    {
      char buf [BUFLEN];
      snprintf (buf, BUFLEN, "'%s' must be a tuple", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }

    if (len > 0 && PyTuple_Size (obj) != len)
    {
      char buf [BUFLEN];
      snprintf (buf, BUFLEN, "'%s' must have %d elements", var, len);
      PyErr_SetString (PyExc_ValueError, buf);
      return 0;
    }
  }

  return 1;
}

/* list test */
static int is_list (PyObject *obj, const char *var, int len)
{
  if (obj)
  {
    if (!PyList_Check (obj))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be a list object", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }

    if (len > 0 && PyList_Size (obj) != len)
    {
      char buf [BUFLEN];
      snprintf (buf, BUFLEN, "'%s' must have %d items", var, len);
      PyErr_SetString (PyExc_ValueError, buf);
      return 0;
    }
  }

  return 1;
}

/* list of tuples test => returns length of the list or zero */
static int is_list_of_tuples (PyObject *obj, const char *var, int min_length, int tuple_length)
{
  if (obj)
  {
    if (!PyList_Check (obj))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be a list object", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }

    int i, j, n = PyList_Size (obj);

    if (n < min_length)
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must have at least %d items", var, min_length);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }

    for (i = 0; i < n; i ++)
    {
      PyObject *item = PyList_GetItem (obj, i);

      if (!PyTuple_Check (item))
      {
	char buf [BUFLEN];
	sprintf (buf, "'%s' must be a list of tuples: item %d is not a tuple", var, i);
	PyErr_SetString (PyExc_ValueError, buf);
	return 0;
      }

      j = PyTuple_Size (item);

      if (j != tuple_length)
      {
	char buf [BUFLEN];
	sprintf (buf, "'%s' list items must be tuples of length %d: item %d has length %d", var, tuple_length, i, j);
	PyErr_SetString (PyExc_ValueError, buf);
	return 0;
      }
    }

    return n;
  }

  return 1;
}

/* list made of one integer folowed by tuples => returns length of the list or zero */
static int is_integer_and_tuples (PyObject *obj, const char *var, int tuple_length)
{
  if (obj)
  {
    if (!PyList_Check (obj))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be a list object", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }

    int i, j, n = PyList_Size (obj);

    if (n < 1)
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must have at least one item", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }

    if (PyInt_AsLong (PyList_GetItem (obj, 0)) <= 0)
    {
	char buf [BUFLEN];
	sprintf (buf, "item 0 must be a positive integer");
	PyErr_SetString (PyExc_ValueError, buf);
	return 0;
    }

    for (i = 1; i < n; i ++)
    {
      PyObject *item = PyList_GetItem (obj, i);

      if (!PyTuple_Check (item))
      {
	char buf [BUFLEN];
	sprintf (buf, "'%s' must be a list of tuples: item %d is not a tuple", var, i);
	PyErr_SetString (PyExc_ValueError, buf);
	return 0;
      }

      j = PyTuple_Size (item);

      if (j != tuple_length)
      {
	char buf [BUFLEN];
	sprintf (buf, "'%s' list items must be tuples of length %d: item %d has length %d", var, tuple_length, i, j);
	PyErr_SetString (PyExc_ValueError, buf);
	return 0;
      }
    }

    return n;
  }

  return 1;
}

/* positive or list of given length test */
static int is_positive_or_list (PyObject *obj, const char *var, int len)
{
  if (PyList_Check (obj))
  {
    if (len > 0 && PyList_Size (obj) != len)
    {
      char buf [BUFLEN];
      snprintf (buf, BUFLEN, "'%s' must have %d items", var, len);
      PyErr_SetString (PyExc_ValueError, buf);
      return 0;
    }
  }
  else
  {
    double num = PyFloat_AsDouble (obj);

    if (num <= 0)
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be positive", var);
      PyErr_SetString (PyExc_ValueError, buf);
      return 0;
    }
  }

  return 1;
}

/* test whether obj is a Python function pointer */
static int is_callable (PyObject *obj, const char *var)
{
  if (obj)
  {
    if (!PyCallable_Check (obj))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be callable", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }
  }

  return 1;
}

/* test whether an object is a list (details as above) or a number */
static int is_list_or_number (PyObject *obj, const char *var, int len)
{
  if (obj)
  {
    if (!(PyList_Check (obj) || PyNumber_Check (obj)))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be a list or a number object", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }

    if (PyList_Check (obj))
    {
      if (!(len > 0 && PyList_Size (obj) != len))
      {
	char buf [BUFLEN];
	sprintf (buf, "'%s' must have %d items", var, len);
	PyErr_SetString (PyExc_ValueError, buf);
	return 0;
      }
    }
  }

  return 1;
}

/* define keywords */
#define KEYWORDS(...) const char *kwl [] = {__VA_ARGS__, NULL}

/* parse arguments with keywords */
#define PARSEKEYS(fmt, ...) if (!PyArg_ParseTupleAndKeywords (args, kwds, fmt, (char**)kwl, __VA_ARGS__)) return NULL

/* parse arguments without keywords */
#define PARSE(fmt, ...) if (!PyArg_ParseTuple (args, fmt, __VA_ARGS__)) return NULL

/* object types assertion */
#define TYPETEST(test) if(!(test)) return NULL

/* string argument if block comparison */
#define IFIS(obj, val) if (strcmp (PyString_AsString (obj), val) == 0)
#define ELIF(obj, val) else if (strcmp (PyString_AsString (obj), val) == 0)
#define ELSE else

/* temporary surface pairing */
struct pair
{
  int color1;
  int color2;
  REAL iparam[NIPARAM];
};

/* material comparison by color pair */
struct cmp
{
  bool operator() (const pair& a, const pair& b)
  {
    if (a.color1 == b.color1) return a.color2 < b.color2;
    else return a.color1 < b.color1;
  }
};

/* sort materials according to surface color pairs */
static void sort_materials ()
{
  std::vector<pair> v;

  v.reserve (pairnum);

  for (int i = 0; i < pairnum; i++)
  {
    pair x;

    x.color1 = pairs[2*i];
    x.color2 = pairs[2*i+1];
    for (int j = 0; j < NIPARAM; j ++) x.iparam[j] = iparam[j][i];

    v.push_back (x);
  }

  std::sort (v.begin(), v.end(), cmp());

  int i = 0;

  for (std::vector<pair>::const_iterator p = v.begin(); p != v.end(); ++p, ++i)
  {
    pairs[2*i] = p->color1;
    pairs[2*i+1] = p->color2;
    for (int j = 0; j < NIPARAM; j ++) iparam[j][i] = p->iparam[j];
  }
}

/* reset simulation */
static PyObject* RESET (PyObject *self, PyObject *args, PyObject *kwds)
{
  master_free (master, parnum);

  slave_free (slave, parnum);

  reset_all_data ();

  Py_RETURN_NONE;
}

/* create material */
static PyObject* MATERIAL (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("density", "young", "poisson");
  double density, young, poisson;

  PARSEKEYS ("ddd", &density, &young, &poisson);

  TYPETEST (is_positive (density, kwl[0]) && is_positive (young, kwl[1]) && is_positive (poisson, kwl[2]));

  if (matnum >= material_buffer_size) material_buffer_grow ();

  int i = matnum ++;

  mparam[DENSITY][i] = density;
  mparam[YOUNG][i] = young;
  mparam[POISSON][i] = poisson;

  return PyInt_FromLong (i);
}

/* create spherical particle */
static PyObject* SPHERE (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("center", "radius", "material", "color");
  int material, color;
  PyObject *cen;
  double rad;

  PARSEKEYS ("Odii", &cen, &rad, &material, &color);

  TYPETEST (is_tuple (cen, kwl[0], 3) && is_positive (rad, kwl[1]) &&
            is_in_range (material, 0, matnum, kwl[2]) && is_positive (color, kwl[3]));

  if (ellnum >= ellipsoid_buffer_size) ellipsoid_buffer_grow ();

  int j = ellnum ++;

  if (parnum >= particle_buffer_size) particle_buffer_grow ();

  int i = parnum ++;

  parmat[i] = material;

  part[j] = i;

  angular[0][i] = 0.0;
  angular[1][i] = 0.0;
  angular[2][i] = 0.0;

  linear[0][i] = 0.0;
  linear[1][i] = 0.0;
  linear[2][i] = 0.0;

  center[0][j] = position[0][i] = PyFloat_AsDouble (PyTuple_GetItem (cen, 0));
  center[1][j] = position[1][i] = PyFloat_AsDouble (PyTuple_GetItem (cen, 1));
  center[2][j] = position[2][i] = PyFloat_AsDouble (PyTuple_GetItem (cen, 2));

  center[3][j] = center[0][j];
  center[4][j] = center[1][j];
  center[5][j] = center[2][j];

  position[3][i] = position[0][i];
  position[4][i] = position[1][i];
  position[5][i] = position[2][i];

  radii[0][j] = rad;
  radii[1][j] = -1.0;
  radii[2][j] = -1.0;

  orient[0][j] = orient[4][j] = orient[8][j] =
  orient[9][j] = orient[13][j] = orient[17][j] = 1.0;
  orient[1][j] = orient[2][j] = orient[3][j] =
  orient[5][j] = orient[6][j] = orient[7][j] =
  orient[10][j] = orient[11][j] = orient[12][j] =
  orient[14][j] = orient[15][j] = orient[16][j] = 0.0;

  ellcol[j] = color;

  double volume = (4./3.)*M_PI*rad*rad*rad;

  mass[i] = volume*mparam[DENSITY][material];

  rotation[0][i] = rotation[4][i] = rotation[8][i] = 1.0;
  rotation[1][i] = rotation[2][i] = rotation[3][i] =
  rotation[5][i] = rotation[6][i] = rotation[7][i] = 0.0;

  inertia[0][i] = inertia[4][i] = inertia[8][i] = 0.4*mass[i]*radii[0][j]*radii[0][j];
  inertia[1][i] = inertia[2][i] = inertia[3][i] =
  inertia[5][i] = inertia[6][i] = inertia[7][i] = 0.0;

  /* return regular particle */
  analytical[i] = 0;

  return PyInt_FromLong (i);
}

/* create meshed particle */
static PyObject* MESH (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("nodes", "elements", "material", "colors");
  int material, *lele, *lsur, i, j, k, l, n, m, nn, o;
  PyObject *nodes, *elements, *colors;
  REAL (*lnod) [3], mi, ci[3], ii[9];
  MESH_DATA *msh;

  PARSEKEYS ("OOdO", &nodes, &elements, &material, &colors);

  TYPETEST (is_list (nodes, kwl[0], 0) && is_list (elements, kwl[1], 0) && is_list_or_number (colors, kwl[3], 0));

  /* test element definitions */
  l = PyList_Size (elements);
  for (i = 0; i < l; i ++)
  {
    k = PyInt_AsLong (PyList_GetItem (elements, i));

    if (!(k == 4 || k == 5 || k == 6 || k == 8))
    {
      PyErr_SetString (PyExc_ValueError, "An element must have 4, 5, 6, or 8 nodes");
      return NULL;
    }

    /* add one more for the
     * material definition */
    i += (k + 1);

    if (i >= l) /* incomplete */
    {
      PyErr_SetString (PyExc_ValueError, "The last element definition is incomplete");
      return NULL;
    }
  }

  /* test nodes list */
  if (PyList_Size(nodes) % 3)
  {
    PyErr_SetString (PyExc_ValueError, "Nodes list length must be a multiple of 3");
    return NULL;
  }

  /* read elements */
  ERRMEM (lele = (int*)malloc ((l + 1) * sizeof (int)));
  nn = PyList_Size (nodes) / 3; /* nodes count */

  for (m = n = i = 0; i < l; m ++)
  {
    lele [n ++] = k = PyInt_AsLong (PyList_GetItem (elements, i ++));

    for (j = 0; j < k; j ++)
    {
      lele [n ++] = PyInt_AsLong (PyList_GetItem (elements, i ++));

      if (lele  [n-1] < 0 || lele [n-1] >= nn) /* must be within the right range */
      {
	char buf [BUFLEN];
	sprintf (buf, "Node %d in element %d is outside of range [0, %d]",j , m, nn-1);
	PyErr_SetString (PyExc_ValueError, buf);
	return NULL;
      }
    }

    /* test for repeated nodes in element definition */
    for (j = 1; j <= k; j ++)
    {
      for (o = j + 1; o <= k; o ++)
      {
	if (lele [n-j] == lele [n-o])
	{
	  char buf [BUFLEN];
	  sprintf (buf, "Nodes %d and %d in element %d are the same", k-j, k-o, m);
	  PyErr_SetString (PyExc_ValueError, buf);
	  return NULL;
	}
      }
    }

    lele [n ++] = PyInt_AsLong (PyList_GetItem (elements, i ++)); /* volid */
  }
  lele [n] = 0; /* end of list */

 if (PyList_Check (colors))
 {
    /* test color definitions */
    l = PyList_Size (colors);
    for (i = 1; i < l; i ++)
    {
      k = PyInt_AsLong (PyList_GetItem (colors, i));

      if (!(k == 3 || k == 4))
      {
	PyErr_SetString (PyExc_ValueError, "A face must have 3 or 4 nodes");
	return NULL;
      }

      /* add one more for the
       * color definition */
      i += (k + 1);

      if (i >= l) /* incomplete */
      {
	PyErr_SetString (PyExc_ValueError, "The last face definition is incomplete");
	return NULL;
      }
    }

    /* read colors */
    ERRMEM (lsur = (int*)malloc ((l + 1) * sizeof (int)));
    lsur [0] = PyInt_AsLong (PyList_GetItem (colors, 0)); /* gcolor */

    for (m = 0, n = i = 1; i < l; m ++)
    {
      lsur [n ++] = k = PyInt_AsLong (PyList_GetItem (colors, i ++));

      for (j = 0; j < k; j ++)
      {
	lsur [n ++] = PyInt_AsLong (PyList_GetItem (colors, i ++));

	if (lsur [n-1] < 0 || lsur [n-1] >= nn) /* must be within the right range */
	{
	  char buf [BUFLEN];
	  sprintf (buf, "Node %d in face %d is outside of range [0, %d]", j, m, nn-1);
	  PyErr_SetString (PyExc_ValueError, buf);
	  return NULL;
	}
      }

      /* test for repeated nodes in face definition */
      for (j = 1; j <= k; j ++)
      {
	for (o = j + 1; o <= k; o ++)
	{
	  if (lsur [n-j] == lsur [n-o])
	  {
	    char buf [BUFLEN];
	    sprintf (buf, "Nodes %d and %d in face %d are the same", k-j, k-o, m);
	    PyErr_SetString (PyExc_ValueError, buf);
	    return NULL;
	  }
	}
      }

      lsur [n ++] = PyInt_AsLong (PyList_GetItem (colors, i ++)); /* surfid */
    }
    lsur [n] = 0; /* end of list */
 }
 else
 {
    ERRMEM (lsur = (int*)malloc (2 * sizeof (int)));
    lsur [0] = PyInt_AsLong (colors);
    lsur [1] = 0; /* end of list */
 }

  /* nodes */
  ERRMEM (lnod = (REAL(*)[3])malloc (nn * sizeof (REAL [3])));
  for (i = 0; i < nn; i ++) 
    for (j = 0; j < 3; j ++)
      lnod [i][j] = PyFloat_AsDouble (PyList_GetItem (nodes, 3*i + j));

  /* create temporary mesh */
  msh = MESH_Create (lnod, lele, lsur);

  /* insert mesh data */

  int node_count = msh->nodes_count;
  int element_node_count = 0;
  ELEMENT *ele;
  for (ele = msh->surfeles; ele; ele = ele->next) element_node_count += ele->type;
  int element_count = msh->surfeles_count;
  int triangle_count = 0;
  FACE *fac;
  for (fac = msh->faces; fac; fac = fac->n) triangle_count += fac->type == 3 ? 1 : 2;
  element_buffer_grow (node_count, element_node_count, element_count, triangle_count);

  for (ele = msh->surfeles, k = elenum, j = eleidx[k]; ele; ele = ele->next, k ++)
  {
    eletype[k] = ele->type;
    for (i = 0; i < ele->type; i ++)
    {
      elenod[j++] = nodnum + ele->nodes[i];
    }
    eleidx[k+1] = eleidx[k]+ele->type;
    elepart[k] = parnum;
    elemat[k] = ele->material;
  }
  elenum += element_count;

  for (k = facnum, fac = msh->faces; fac; fac = fac->n)
  {
    /* insert triangle */
    if (trinum >= triangle_buffer_size) triangle_buffer_grow ();
    l = trinum ++;
    tri[0][0][l] = msh->nodes[fac->nodes[0]][0];
    tri[0][1][l] = msh->nodes[fac->nodes[0]][1];
    tri[0][2][l] = msh->nodes[fac->nodes[0]][2];
    tri[1][0][l] = msh->nodes[fac->nodes[1]][0];
    tri[1][1][l] = msh->nodes[fac->nodes[1]][1];
    tri[1][2][l] = msh->nodes[fac->nodes[1]][2];
    tri[2][0][l] = msh->nodes[fac->nodes[2]][0];
    tri[2][1][l] = msh->nodes[fac->nodes[2]][1];
    tri[2][2][l] = msh->nodes[fac->nodes[2]][2];
    tricol[l] = fac->color;
    triobs[l] = parnum;

    /* insert face */
    facnod[0][k] = nodnum + fac->nodes[0];
    facnod[1][k] = nodnum + fac->nodes[1];
    facnod[2][k] = nodnum + fac->nodes[2];
    facpart[k] = parnum;
    factri[k] = l;
    k ++;

    if (fac->type == 4) /* insert second triangle and face */
    {
      if (trinum >= triangle_buffer_size) triangle_buffer_grow ();
      l = trinum ++;
      tri[0][0][l] = msh->nodes[fac->nodes[0]][0];
      tri[0][1][l] = msh->nodes[fac->nodes[0]][1];
      tri[0][2][l] = msh->nodes[fac->nodes[0]][2];
      tri[1][0][l] = msh->nodes[fac->nodes[2]][0];
      tri[1][1][l] = msh->nodes[fac->nodes[2]][1];
      tri[1][2][l] = msh->nodes[fac->nodes[2]][2];
      tri[2][0][l] = msh->nodes[fac->nodes[3]][0];
      tri[2][1][l] = msh->nodes[fac->nodes[3]][1];
      tri[2][2][l] = msh->nodes[fac->nodes[3]][2];
      tricol[l] = fac->color;
      triobs[l] = parnum;

      facnod[0][k] = nodnum + fac->nodes[0];
      facnod[1][k] = nodnum + fac->nodes[2];
      facnod[2][k] = nodnum + fac->nodes[3];
      facpart[k] = parnum;
      factri[k] = l;
      k ++;
    }
  }
  facnum += triangle_count;

  for (j = 0; j < node_count; j ++)
  {
    k = nodnum+j;
    parmec::nodes[0][k] = parmec::nodes[3][k] = msh->nodes[j][0];
    parmec::nodes[1][k] = parmec::nodes[4][k] = msh->nodes[j][1];
    parmec::nodes[2][k] = parmec::nodes[5][k] = msh->nodes[j][2];
  }
  nodnum += node_count;

  /* insert particle data */

  i = parnum ++;

  parmat[i] = material;

  angular[0][i] = 0.0;
  angular[1][i] = 0.0;
  angular[2][i] = 0.0;

  linear[0][i] = 0.0;
  linear[1][i] = 0.0;
  linear[2][i] = 0.0;

  MESH_Char (msh, &mi, ci, ii);

  position[0][i] = ci[0];
  position[1][i] = ci[1];
  position[2][i] = ci[2];

  position[3][i] = position[0][i];
  position[4][i] = position[1][i];
  position[5][i] = position[2][i];

  mass[i] = mi;

  rotation[0][i] = rotation[4][i] = rotation[8][i] = 1.0;
  rotation[1][i] = rotation[2][i] = rotation[3][i] =
  rotation[5][i] = rotation[6][i] = rotation[7][i] = 0.0;

  inertia[0][i] = ii[0];
  inertia[1][i] = ii[1];
  inertia[2][i] = ii[2];
  inertia[3][i] = ii[3];
  inertia[4][i] = ii[4];
  inertia[5][i] = ii[5];
  inertia[6][i] = ii[6];
  inertia[7][i] = ii[7];
  inertia[8][i] = ii[8];

  /* clean up */
  MESH_Destroy (msh);
  free (lnod);
  free (lele);
  free (lsur);

  /* return regular particle */
  analytical[i] = 0;

  return PyInt_FromLong (i);
}

/* create analytical particle */
static PyObject* ANALYTICAL (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("inertia", "mass", "rotation", "position", "material", "particle");
  PyObject *inertia, *rotation, *position;
  int material, particle, i;
  REAL vinertia[6], vrotation[9], vposition[3];
  double mass;

  inertia = NULL;
  mass = 0.0;
  rotation = NULL;
  position = NULL;
  material = 0;
  particle = -1;

  PARSEKEYS ("|OdOOi", &inertia, &mass, &rotation, &position, &particle);

  TYPETEST (is_list (inertia, kwl[0], 6) && is_positive (mass, kwl[1]) &&
            is_list (rotation, kwl[2], 9) && is_tuple (position, kwl[3], 3) &&
	    is_non_negative (material, kwl[4]) && is_non_negative (particle, kwl[4]));

  if (inertia)
  {
    vinertia[0] = PyFloat_AsDouble (PyList_GetItem (inertia, 0));
    vinertia[1] = PyFloat_AsDouble (PyList_GetItem (inertia, 1));
    vinertia[2] = PyFloat_AsDouble (PyList_GetItem (inertia, 2));
    vinertia[3] = PyFloat_AsDouble (PyList_GetItem (inertia, 3));
    vinertia[4] = PyFloat_AsDouble (PyList_GetItem (inertia, 4));
    vinertia[5] = PyFloat_AsDouble (PyList_GetItem (inertia, 5));
  }
  else
  {
    vinertia[0] = 1.0;
    vinertia[1] = 1.0;
    vinertia[2] = 1.0;
    vinertia[3] = 0.0;
    vinertia[4] = 0.0;
    vinertia[5] = 0.0;
  }

  if (rotation)
  {
    vrotation[0] = PyFloat_AsDouble (PyList_GetItem (rotation, 0));
    vrotation[1] = PyFloat_AsDouble (PyList_GetItem (rotation, 1));
    vrotation[2] = PyFloat_AsDouble (PyList_GetItem (rotation, 2));
    vrotation[3] = PyFloat_AsDouble (PyList_GetItem (rotation, 3));
    vrotation[4] = PyFloat_AsDouble (PyList_GetItem (rotation, 4));
    vrotation[5] = PyFloat_AsDouble (PyList_GetItem (rotation, 5));
    vrotation[6] = PyFloat_AsDouble (PyList_GetItem (rotation, 6));
    vrotation[7] = PyFloat_AsDouble (PyList_GetItem (rotation, 7));
    vrotation[8] = PyFloat_AsDouble (PyList_GetItem (rotation, 8));
  }
  else
  {
    vrotation[0] = 1.0;
    vrotation[1] = 0.0;
    vrotation[2] = 0.0;
    vrotation[3] = 0.0;
    vrotation[4] = 1.0;
    vrotation[5] = 0.0;
    vrotation[6] = 0.0;
    vrotation[7] = 0.0;
    vrotation[8] = 1.0;
  }

  if (position)
  {
    vposition[0] = PyFloat_AsDouble (PyTuple_GetItem (position, 0));
    vposition[1] = PyFloat_AsDouble (PyTuple_GetItem (position, 1));
    vposition[2] = PyFloat_AsDouble (PyTuple_GetItem (position, 2));
  }
  else
  {
    vposition[0] = 0.0;
    vposition[1] = 0.0;
    vposition[2] = 0.0;
  }

  if (particle < 0) /* create new analytical particle */
  {
    i = particle = parnum ++;

    parmat[i] = material;

    angular[0][i] = 0.0;
    angular[1][i] = 0.0;
    angular[2][i] = 0.0;

    linear[0][i] = 0.0;
    linear[1][i] = 0.0;
    linear[2][i] = 0.0;

    parmec::position[0][i] = vposition[0];
    parmec::position[1][i] = vposition[1];
    parmec::position[2][i] = vposition[2];
    parmec::position[3][i] = vposition[0];
    parmec::position[4][i] = vposition[1];
    parmec::position[5][i] = vposition[2];

    if (mass > 0.0) parmec::mass[i] = mass;

    parmec::rotation[0][i] = vrotation[0];
    parmec::rotation[1][i] = vrotation[1];
    parmec::rotation[2][i] = vrotation[2]; 
    parmec::rotation[3][i] = vrotation[3];
    parmec::rotation[4][i] = vrotation[4];
    parmec::rotation[5][i] = vrotation[5];
    parmec::rotation[6][i] = vrotation[6];
    parmec::rotation[7][i] = vrotation[7];
    parmec::rotation[8][i] = vrotation[8];

    parmec::inertia[0][i] = vinertia[0];
    parmec::inertia[4][i] = vinertia[1];
    parmec::inertia[8][i] = vinertia[2];
    parmec::inertia[1][i] = parmec::inertia[3][i] = vinertia[3];
    parmec::inertia[2][i] = parmec::inertia[6][i] = vinertia[4];
    parmec::inertia[5][i] = parmec::inertia[7][i] = vinertia[5];
  }
  else /* redefine an existing particle as analytical */
  {
    if (particle >= parnum)
    {
      PyErr_SetString (PyExc_ValueError, "Particle index is out of range");
      return NULL;
    }

    declare_analytical (particle);

    if (position)
    {
      parmec::position[0][i] = vposition[0];
      parmec::position[1][i] = vposition[1];
      parmec::position[2][i] = vposition[2];
      parmec::position[3][i] = vposition[0];
      parmec::position[4][i] = vposition[1];
      parmec::position[5][i] = vposition[2];
    }

    if (mass > 0.0) parmec::mass[i] = mass;

    if (rotation)
    {
      parmec::rotation[0][i] = vrotation[0];
      parmec::rotation[1][i] = vrotation[1];
      parmec::rotation[2][i] = vrotation[2]; 
      parmec::rotation[3][i] = vrotation[3];
      parmec::rotation[4][i] = vrotation[4];
      parmec::rotation[5][i] = vrotation[5];
      parmec::rotation[6][i] = vrotation[6];
      parmec::rotation[7][i] = vrotation[7];
      parmec::rotation[8][i] = vrotation[8];
    }

    if (inertia)
    {
      parmec::inertia[0][i] = vinertia[0];
      parmec::inertia[4][i] = vinertia[1];
      parmec::inertia[8][i] = vinertia[2];
      parmec::inertia[1][i] = parmec::inertia[3][i] = vinertia[3];
      parmec::inertia[2][i] = parmec::inertia[6][i] = vinertia[4];
      parmec::inertia[5][i] = parmec::inertia[7][i] = vinertia[5];
    }
  }

  /* return analytical particle */
  analytical[particle] = 1;

  return PyInt_FromLong (particle);
}

/* create obstacle */
static PyObject* OBSTACLE (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("triangles", "color", "point", "linear", "angular");
  PyObject *triangles, *color, *point, *lin, *ang;

  point = NULL;
  lin = NULL;
  ang = NULL;

  PARSEKEYS ("OO|OOO", &triangles, &color, &point, &lin, &ang);

  TYPETEST (is_list_of_tuples (triangles, kwl[0], 1, 9) &&
            is_positive_or_list (color, kwl[1], PyList_Size(triangles)) &&
            is_tuple (point, kwl[2], 3) && is_callable (lin, kwl[3]) &&
	    is_callable (ang, kwl[4]));

  int m = PyList_Size (triangles);

  int haveobs = 0;

  if (point && (lin || ang))
  {
    if (obsnum >= obstacle_buffer_size) obstacle_buffer_grow ();

    int i = obsnum ++;

    trirng[2*i] = trinum;
    trirng[2*i+1] = trinum+m;

    obspnt[3*i] = PyFloat_AsDouble (PyTuple_GetItem (point, 0)); 
    obspnt[3*i+1] = PyFloat_AsDouble (PyTuple_GetItem (point, 1)); 
    obspnt[3*i+2] = PyFloat_AsDouble (PyTuple_GetItem (point, 2)); 

    linhis[i] = lin;
    anghis[i] = ang;

    haveobs = 1;
  }

  for (int n = 0; n < m; n ++)
  {
    PyObject *t = PyList_GetItem (triangles, n);

    if (trinum >= triangle_buffer_size) triangle_buffer_grow ();

    int i = trinum ++;

    tri[0][0][i] = PyFloat_AsDouble (PyTuple_GetItem (t, 0)); 
    tri[0][1][i] = PyFloat_AsDouble (PyTuple_GetItem (t, 1)); 
    tri[0][2][i] = PyFloat_AsDouble (PyTuple_GetItem (t, 2)); 
    tri[1][0][i] = PyFloat_AsDouble (PyTuple_GetItem (t, 3)); 
    tri[1][1][i] = PyFloat_AsDouble (PyTuple_GetItem (t, 4)); 
    tri[1][2][i] = PyFloat_AsDouble (PyTuple_GetItem (t, 5)); 
    tri[2][0][i] = PyFloat_AsDouble (PyTuple_GetItem (t, 6)); 
    tri[2][1][i] = PyFloat_AsDouble (PyTuple_GetItem (t, 7)); 
    tri[2][2][i] = PyFloat_AsDouble (PyTuple_GetItem (t, 8)); 

    if (PyList_Check (color))
    {
      tricol[i] = PyInt_AsLong (PyList_GetItem (color, n));
    }
    else
    {
      tricol[i] = PyInt_AsLong (color);
    }

    triobs[i] = haveobs ? -obsnum-2 : -1; /* <0 - moving obstalce (-index-2), -1 - static obstacle, >= 0 - triangulated particle */
  }

  Py_RETURN_NONE;
}

/* create translational spring constraint */
static PyObject* SPRING (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("part1", "point1", "part2", "point2", "spring", "dashpot", "direction");
  PyObject *point1, *point2, *spring, *dashpot, *direction;
  int part1, part2;

  direction = NULL;
  dashpot = NULL;

  PARSEKEYS ("dOdOO|OO", &part1, &point1, &part2, &point2, &spring, &dashpot, &direction);

  TYPETEST (is_non_negative (part1, kwl[0]) && is_tuple (point1, kwl[1], 3) &&
            is_tuple (point2, kwl[3], 3) && is_list (spring, kwl[4], 0) &&
	    is_list (dashpot, kwl[4], 0) && is_tuple (direction, kwl[5], 3));

  if (part2 < -1)
  {
    PyErr_SetString (PyExc_ValueError, "Particle two index is out of range [-1, 0, ...]");
    return NULL;
  }

  if (PyList_Size (spring) < 4 || PyList_Size (spring) % 2)
  {
    PyErr_SetString (PyExc_ValueError, "Invalid spring lookup table list length");
    return NULL;
  }

  if (PyList_Size (dashpot) < 4 || PyList_Size (dashpot) % 2)
  {
    PyErr_SetString (PyExc_ValueError, "Invalid dashpot lookup table list length");
    return NULL;
  }

  int spring_lookup = PyList_Size (spring);

  int dashpot_lookup = dashpot ? PyList_Size (dashpot) : 4;

  spring_buffer_grow (spring_lookup, dashpot_lookup);

  int i = sprnum ++;

  sprpart[0][i] = part1;
  sprpart[1][i] = part2;

  sprpnt[0][0][i] = PyFloat_AsDouble (PyTuple_GetItem (point1,0));
  sprpnt[0][1][i] = PyFloat_AsDouble (PyTuple_GetItem (point1,1));
  sprpnt[0][2][i] = PyFloat_AsDouble (PyTuple_GetItem (point1,2));
  sprpnt[0][3][i] = sprpnt[0][0][i];
  sprpnt[0][4][i] = sprpnt[0][1][i];
  sprpnt[0][5][i] = sprpnt[0][2][i];
  sprpnt[1][0][i] = PyFloat_AsDouble (PyTuple_GetItem (point2,0));
  sprpnt[1][1][i] = PyFloat_AsDouble (PyTuple_GetItem (point2,1));
  sprpnt[1][2][i] = PyFloat_AsDouble (PyTuple_GetItem (point2,2));
  sprpnt[1][3][i] = sprpnt[1][0][i];
  sprpnt[1][4][i] = sprpnt[1][1][i];
  sprpnt[1][5][i] = sprpnt[1][2][i];

  int j = 0, k = spridx[i];

  for (; j < spring_lookup/2; j ++, k ++)
  {
    REAL stroke = PyFloat_AsDouble(PyList_GetItem(spring,2*j));
    REAL force = PyFloat_AsDouble(PyList_GetItem(spring,2*j+1));
    parmec::spring[0][k] = stroke;
    parmec::spring[1][k] = force;
  }
  spridx[sprnum] = k;

  if (dashpot)
  {
    for (j = 0, k = dashidx[i]; j < dashpot_lookup/2; j ++, k ++)
    {
      REAL velocity = PyFloat_AsDouble(PyList_GetItem(dashpot,2*j));
      REAL force = PyFloat_AsDouble(PyList_GetItem(dashpot,2*j+1));
      parmec::dashpot[0][k] = velocity;
      parmec::dashpot[1][k] = force;
    }
    dashidx[sprnum] = k;
  }
  else /* default zero force */
  {
    k = dashidx[i];
    parmec::dashpot[0][k] = -REAL_MAX;
    parmec::dashpot[1][k] = 0.0;
    parmec::dashpot[0][k+1] = +REAL_MAX;
    parmec::dashpot[1][k+1] = 0.0;
    dashidx[sprnum] = k+2;
  }

  if (direction)
  {
    REAL dir[3];

    dir[0] = PyFloat_AsDouble(PyTuple_GetItem(direction,0));
    dir[1] = PyFloat_AsDouble(PyTuple_GetItem(direction,1));
    dir[2] = PyFloat_AsDouble(PyTuple_GetItem(direction,2));

    REAL len = LEN(dir);

    if (len == 0.0)
    {
      PyErr_SetString (PyExc_ValueError, "Invalid zero direction");
      return NULL;
    }

    REAL inv = 1.0/len;

    dir[0] *= inv;
    dir[1] *= inv;
    dir[2] *= inv;

    sprdir[0][i] = dir[0];
    sprdir[1][i] = dir[1];
    sprdir[2][i] = dir[2];

    sprdirup[i] = 0; /* constant direction */

    REAL dif[3] = {sprpnt[1][0][i] - sprpnt[0][0][i],
                   sprpnt[1][1][i] - sprpnt[0][1][i],
                   sprpnt[1][2][i] - sprpnt[0][2][i]};

    stroke0[i] = DOT (dif, dir);
  }
  else /* direction = (p2 - p1)/|p2 - p1| */
  {
    sprdirup[i] = 1;

    REAL dif[3] = {sprpnt[1][0][i] - sprpnt[0][0][i],
                   sprpnt[1][1][i] - sprpnt[0][1][i],
                   sprpnt[1][2][i] - sprpnt[0][2][i]};

    stroke0[i] = LEN (dif);
  }

  Py_RETURN_NONE;
}

/* define surface pairing for the granular interaction model */
static PyObject* GRANULAR (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("color1", "color2", "spring", "damper", "friction", "rolling", "drilling", "kskn");
  int color1, color2;
  double spring, damper, rolling, drilling, kskn;
  PyObject *friction;
  double fri[2];

  damper = 1.0;
  rolling = 0.0;
  drilling = 0.0;
  kskn = 0.5;
  friction = NULL;
  fri[0] = 0.0;
  fri[1] = 0.0;

  PARSEKEYS ("iid|dOddd", &color1, &color2, &spring, &damper, &friction, &rolling, &drilling, &kskn);

  TYPETEST (is_positive (spring, kwl[2]) && is_non_negative (damper, kwl[3]) &&
            is_non_negative (rolling, kwl[5]) && is_non_negative (drilling, kwl[6]) &&
	    is_positive (kskn, kwl[7]));

  if (PyTuple_Check (friction))
  {
    TYPETEST (is_tuple (friction, kwl[5], 2));

    fri[0] = PyFloat_AsDouble (PyTuple_GetItem (friction,0));
    fri[1] = PyFloat_AsDouble (PyTuple_GetItem (friction,1));

    TYPETEST (is_non_negative (fri[0], kwl[4]) && is_non_negative (fri[1], kwl[4]));
  }
  else
  {
    fri[0] = PyFloat_AsDouble (friction);

    TYPETEST (is_non_negative (fri[0], kwl[4]));

    fri[1] = fri[0];
  }

  int i = -1;

  if (color1 == 0 && color2 == 0) /* default */
  {
    i = 0;
  }
  else if (color1 > 0 && color2 > 0)
  {
    if (pairnum >= pair_buffer_size) pair_buffer_grow ();

    int i = pairnum ++;

    pairs[2*i] = std::min(color1, color2);
    pairs[2*i+1] = std::max(color1, color2);
  }
  else
  {
    PyErr_SetString (PyExc_ValueError, "invalid color values");
    return NULL;
  }

  ikind[i] = GRANULAR_FORCE;
  iparam[parmec::SPRING][i] = spring;
  iparam[DAMPER][i] = damper;
  iparam[FRISTAT][i] = fri[0];
  iparam[FRIDYN][i] = fri[1];
  iparam[FRIROL][i] = rolling;
  iparam[FRIDRIL][i] = drilling;
  iparam[KSKN][i] = kskn;

  Py_RETURN_NONE;
}

/* set particle velocity */
static PyObject* VELOCITY (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("parnum", "linear", "angular");
  PyObject *lin, *ang;
  int i;

  ang = NULL;

  PARSEKEYS ("iO|O", &i, &lin, &ang);

  TYPETEST (is_in_range (i, 0, parnum, kwl[0]) && is_tuple (lin, kwl[1], 3) && is_tuple (ang, kwl[2], 3));

  linear[0][i] = PyFloat_AsDouble (PyTuple_GetItem (lin, 0)); 
  linear[1][i] = PyFloat_AsDouble (PyTuple_GetItem (lin, 1)); 
  linear[2][i] = PyFloat_AsDouble (PyTuple_GetItem (lin, 2)); 

  if (ang)
  {
    angular[0][i] = PyFloat_AsDouble (PyTuple_GetItem (ang, 0)); 
    angular[1][i] = PyFloat_AsDouble (PyTuple_GetItem (ang, 1)); 
    angular[2][i] = PyFloat_AsDouble (PyTuple_GetItem (ang, 2)); 
  }

  Py_RETURN_NONE;
}

/* set gravity */
static PyObject* GRAVITY (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("gx", "gy", "gz");
  double gx, gy, gz;

  PARSEKEYS ("ddd", &gx, &gy, &gz);

  gravity[0] = gx;
  gravity[1] = gy;
  gravity[2] = gz;

  Py_RETURN_NONE;
}

/* progress bar by Ross Hemsley;
 * http://www.rosshemsley.co.uk/2011/02/creating-a-progress-bar-in-c-or-any-other-console-app/ */
static void progressbar (unsigned int x, unsigned int n, unsigned int w = 50)
{
  if (n < 100)
  {
    x *= 100/n;
    n = 100;
  }

  if ((x != n) && (x % (n/100) != 0)) return;

  using namespace std;
       
  float ratio  =  x/(float)n;
  int c =  ratio * w;
	        
  cout << setw(3) << (int)(ratio*100) << "% [";
  for (int x=0; x<c; x++) cout << "=";
  for (int x=c; x<w; x++) cout << " ";
  cout << "]\r" << flush;
}

/* ISPC calls are used below */
using namespace ispc;

/* estimate critical time step */
static PyObject* CRITICAL (PyObject *self, PyObject *args, PyObject *kwds)
{
  REAL h = critical (parnum, mass, pairnum, iparam);

  return PyFloat_FromDouble (h);
}

/* update obstacles time histories from callbacks */
static void obstaclev (int obsnum, REAL *obsang, REAL *obslin, callback_t anghis[], callback_t linhis[], REAL time)
{
  PyObject *result, *args;
  int i;

  args = Py_BuildValue ("(d)", time);

  for (i = 0; i < obsnum; i ++, obsang += 3, obslin += 3)
  {
    if (anghis[i])
    {
      result = PyObject_CallObject ((PyObject*)anghis[i], args);

      ASSERT (is_tuple (result, "Returned value", 3), "Obstacle angular velocity callback did not return a (ox, oy, oz) tuple");

      obsang[0] = PyFloat_AsDouble(PyTuple_GetItem (result, 0));
      obsang[1] = PyFloat_AsDouble(PyTuple_GetItem (result, 1));
      obsang[2] = PyFloat_AsDouble(PyTuple_GetItem (result, 2));

      Py_DECREF (result);
    }
    else
    {
      SET (obsang, 0.0);
    }

    if (linhis[i])
    {
      result = PyObject_CallObject ((PyObject*)linhis[i], args);

      ASSERT (is_tuple (result, "Returned value", 3), "Obstacle linear velocity callback did not return a (vx, vy, vz) tuple");

      obslin[0] = PyFloat_AsDouble(PyTuple_GetItem (result, 0));
      obslin[1] = PyFloat_AsDouble(PyTuple_GetItem (result, 1));
      obslin[2] = PyFloat_AsDouble(PyTuple_GetItem (result, 2));

      Py_DECREF (result);
    }
    else
    {
      SET (obslin, 0.0);
    }
  }

  Py_DECREF (args);
}

/* run DEM simulation */
static PyObject* DEM (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("duration", "step", "interval", "prefix");
  double duration, step, interval, time, t0;
  PyObject *prefix;
  timing tt;

  prefix = NULL;

  PARSEKEYS ("dd|dO", &duration, &step, &interval, &prefix);

  TYPETEST (is_positive (duration, kwl[0]) && is_positive (step, kwl[1]) &&
            is_positive (interval, kwl[2]) && is_string (prefix, kwl[3]));

  timerstart (&tt);

  if (prefix)
  {
    char *pre, *out;
    int len, i;
    len = strlen (outpath);
    for (i = len-1; i >= 0; i --)
    {
      if (outpath[i] == '/' || /* unix */
	  outpath[i] == '\\') break; /* windows */
    }
    pre = PyString_AsString (prefix);
    len = (i+1) + strlen (pre) + 1;
    ERRMEM (out = new char [len]);
    outpath[i+1] = '\0';
    strcpy (out, outpath);
    strcpy (out+i+1, pre);
    delete outpath;
    outpath = out;
  }

  sort_materials ();

  if (curtime == 0.0)
  {
    output ();

    euler (threads, parnum, angular, linear, rotation, position, 0.5*step);

    shapes (threads, ellnum, part, center, radii, orient, nodnum, nodes,
            nodpart, NULL, facnum, facnod, factri, tri, rotation, position);

    obstaclev (obsnum, obsang, obslin, anghis, linhis, 0.5*step);

    obstacles (obsnum, trirng, obspnt, obsang, obslin, tri, step);
  }

  invert_inertia (threads, parnum, inertia, inverse, mass, invm);

  zero_force_and_torque (parnum, force, torque);

  partitioning *tree = partitioning_create (threads, ellnum, center);

  /* time stepping */
  for (t0 = time = 0.0; time < duration; time += step)
  {
    REAL *icenter[6] = {center[0]+ellcon, center[1]+ellcon, center[2]+ellcon, center[3]+ellcon, center[4]+ellcon, center[5]+ellcon};
    REAL *iradii[3] = {radii[0]+ellcon, radii[1]+ellcon, radii[2]+ellcon};
    REAL *iorient[18] = {orient[0]+ellcon, orient[1]+ellcon, orient[2]+ellcon, orient[3]+ellcon, orient[4]+ellcon, orient[5]+ellcon,
                         orient[6]+ellcon, orient[7]+ellcon, orient[8]+ellcon, orient[9]+ellcon, orient[10]+ellcon, orient[11]+ellcon,
			 orient[12]+ellcon, orient[13]+ellcon, orient[14]+ellcon, orient[15]+ellcon, orient[16]+ellcon, orient[17]+ellcon};
    REAL *itri[3][3] = {{tri[0][0]+tricon, tri[0][1]+tricon, tri[0][2]+tricon},
                        {tri[1][0]+tricon, tri[1][1]+tricon, tri[1][2]+tricon},
			{tri[2][0]+tricon, tri[2][1]+tricon, tri[2][2]+tricon}};
    int *ifacnod[3] = {facnod[0]+faccon, facnod[1]+faccon, facnod[2]+faccon};

    if (partitioning_store (threads, tree, ellnum-ellcon, ellcol+ellcon, part+ellcon, icenter, iradii, iorient) > 0)
    {
      partitioning_destroy (tree);

      tree = partitioning_create (threads, ellnum-ellcon, icenter);

      ASSERT (partitioning_store (threads, tree, ellnum-ellcon, ellcol+ellcon, part+ellcon, icenter, iradii, iorient) == 0, "Repartitioning failed");
    }

    condet (threads, tree, master, parnum, ellnum-ellcon, ellcol+ellcon, part+ellcon,
            icenter, iradii, iorient, trinum-tricon, tricol+tricon, triobs+tricon, itri);

    contacts (threads, master, slave, parnum, angular, linear, rotation, position, inertia, inverse,
              mass, invm, obspnt, obslin, obsang, parmat, mparam, pairnum, pairs, ikind, iparam, step);

    springs (threads, sprnum, sprpart, sprpnt, spring, spridx, dashpot, dashidx, sprdir,
             sprdirup, stroke0, angular, linear, rotation, position, force, torque);

    dynamics (threads, master, slave, parnum, angular, linear, rotation,
      position, inertia, inverse, mass, invm, force, torque, gravity, step);

    if (time >= t0 + interval) /* full update, due to output */
    {
      shapes (threads, ellnum, part, center, radii, orient,
	      nodnum, nodes, nodpart, NULL, facnum, facnod,
	      factri, tri, rotation, position);
    }
    else /* partial update, skipping analytical particles */
    {
      shapes (threads, ellnum-ellcon, part+ellcon, icenter, iradii, iorient,
	      nodnum, nodes, nodpart, analytical, facnum-faccon, ifacnod,
	      factri+faccon, tri, rotation, position);
    }

    obstaclev (obsnum, obsang, obslin, anghis, linhis, time+step);

    obstacles (obsnum, trirng, obspnt, obsang, obslin, tri, step);

    if (time >= t0 + interval)
    {
      output ();

      t0 += interval;
    }

    progressbar (time/step, duration/step);

    curtime += step;
  }

  partitioning_destroy (tree);

  double dt = timerend (&tt);

  printf("[ ===             %10.3f sec                    === ]\n", dt);

  return Py_BuildValue ("d", dt); /* PyFloat_FromDouble (dt) */
}

static PyMethodDef methods [] =
{
  {"RESET", (PyCFunction)RESET, METH_NOARGS, "Reset simulation"},
  {"MATERIAL", (PyCFunction)MATERIAL, METH_VARARGS|METH_KEYWORDS, "Create material"},
  {"SPHERE", (PyCFunction)SPHERE, METH_VARARGS|METH_KEYWORDS, "Create spherical particle"},
  {"MESH", (PyCFunction)MESH, METH_VARARGS|METH_KEYWORDS, "Create meshed particle"},
  {"ANALYTICAL", (PyCFunction)ANALYTICAL, METH_VARARGS|METH_KEYWORDS, "Create analytical particle"},
  {"OBSTACLE", (PyCFunction)OBSTACLE, METH_VARARGS|METH_KEYWORDS, "Create obstacle"},
  {"SPRING", (PyCFunction)::SPRING, METH_VARARGS|METH_KEYWORDS, "Create translational spring"},
  {"GRANULAR", (PyCFunction)GRANULAR, METH_VARARGS|METH_KEYWORDS, "Define surface pairing for the granular interaction model"},
  {"VELOCITY", (PyCFunction)VELOCITY, METH_VARARGS|METH_KEYWORDS, "Set particle velocity"},
  {"GRAVITY", (PyCFunction)GRAVITY, METH_VARARGS|METH_KEYWORDS, "Set gravity"},
  {"CRITICAL", (PyCFunction)CRITICAL, METH_NOARGS, "Estimate critical time step"},
  {"DEM", (PyCFunction)DEM, METH_VARARGS|METH_KEYWORDS, "Run DEM simulation"},
  {NULL, 0, 0, NULL}
};

/* interpret an input file (return 0 on success) */
int input (const char *path)
{
  char *line;
  int error;

  Py_Initialize();

  if (!Py_InitModule3 ("parmec", methods, "parmec module")) return -1;

  PyRun_SimpleString ("from parmec import RESET\n"
                      "from parmec import MATERIAL\n"
                      "from parmec import SPHERE\n"
                      "from parmec import OBSTACLE\n"
                      "from parmec import SPRING\n"
                      "from parmec import GRANULAR\n"
                      "from parmec import VELOCITY\n"
                      "from parmec import GRAVITY\n"
                      "from parmec import CRITICAL\n"
                      "from parmec import DEM\n");

  ERRMEM (line = new char [128 + strlen (path)]);
  sprintf (line, "execfile ('%s')", path);

  error = PyRun_SimpleString (line); /* we do not run a file directly because FILE destriptors differe
					between WIN32 and UNIX while Python is often provided in binary form */
  delete line;

  return error;
}

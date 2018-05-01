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
#include "timeseries.h"
#include "macros.h"
#include "parmec.h"
#include "output.h"
#include "timer.h"
#include "mesh.h"
#include "constants.h"
#include "parmec_ispc.h"
#include "partition_ispc.h"
#include "forces_ispc.h"
#include "dynamics_ispc.h"
#include "shapes_ispc.h"
#include "obstacles_ispc.h"
#include "constrain_ispc.h"

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

/* bool test */
static int is_bool (PyObject *obj, const char *var)
{
  if (obj)
  {
    if (!PyBool_Check (obj))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be Boolean (True/False)", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }
  }

  return 1;
}

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

/* file path test */
static int is_file (PyObject *obj, const char *var)
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

    FILE *f = fopen (PyString_AsString (obj), "r");
    if (!f)
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s = %s' must be an existing file path", var, PyString_AsString(obj));
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }
    fclose (f);
  }

  return 1;
}

/* string or list test */
static int is_string_or_list (PyObject *obj, const char *var)
{
  if (obj)
  {
    if (!(PyString_Check (obj) || PyList_Check(obj)))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be a string or a list", var);
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
static int is_ge_lt (double num, double lo, double hi, const char *var)
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

/* in [lo, hi] range test */
static int is_ge_le (double num, double lo, double hi, const char *var)
{
  if (num < lo || num > hi)
  {
    char buf [BUFLEN];
    sprintf (buf, "'%s' must belong to [%g, %g]", var, lo, hi);
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

/* tuple or list of tuples check */
static int is_tuple_or_list_of_tuples (PyObject *obj, const char *var, int tuple_length, int min_length)
{
  if (obj)
  {
    if (PyTuple_Check (obj))
    {
      if (tuple_length > 0 && PyTuple_Size (obj) != tuple_length)
      {
	char buf [BUFLEN];
	snprintf (buf, BUFLEN, "tuple '%s' must have %d elements", var, tuple_length);
	PyErr_SetString (PyExc_ValueError, buf);
	return 0;
      }
    }
    else if (PyList_Check (obj))
    {
      int i, j, n = PyList_Size (obj);

      if (n < min_length)
      {
	char buf [BUFLEN];
	sprintf (buf, "list '%s' must have at least %d items", var, min_length);
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
    else
    {
      char buf [BUFLEN];
      snprintf (buf, BUFLEN, "'%s' must be a tuple or a list of tuples", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }
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
      if (len > 0 && PyList_Size (obj) != len)
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

/* test whether an object is a number or a list (details as above) or a string */
static int is_number_or_list_or_string (PyObject *obj, const char *var, int div, int len)
{
  if (obj)
  {
    if (!(PyNumber_Check(obj) || PyString_Check (obj) || PyList_Check (obj)))
    {
      char buf [BUFLEN];
      sprintf (buf, "'%s' must be a number or a list or a string object", var);
      PyErr_SetString (PyExc_TypeError, buf);
      return 0;
    }

    if (PyList_Check (obj))
    {
      if (!(PyList_Size (obj) % div == 0 && PyList_Size (obj) >= len))
      {
	char buf [BUFLEN];
	sprintf (buf, "'%s' must have N * %d elements, where N >= %d", var, div, len / div);
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

/* test if a string has a specific ending */
static int endswith (const char *string, const char *ending)
{
  if (strlen (string) >= strlen (ending) &&
      strcmp (string+strlen(string)-strlen(ending), ending) == 0) return 1;
  else return 0;
}

/* command line arguments */
static PyObject* ARGV (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("nonparmec");
  PyObject *nonparmec, *list;

  nonparmec = Py_True;

  PARSEKEYS ("|O", &nonparmec);

  TYPETEST (is_bool (nonparmec, kwl [0]));

  if (!(list = PyList_New (0))) return NULL;

  if (nonparmec == Py_True)
  {
    for (int i = 0; i < argc; i ++)
    {
      if (endswith (argv[i], "parmec4")) continue;
      else if (endswith (argv[i], "parmec8")) continue;
      else if (strcmp (argv[i], "-ntasks") == 0)
      {
	i ++;
	continue;
      }
      else if (strlen (argv[i]) > 3 && strcmp(argv[i]+strlen(argv[i])-3, ".py") == 0) continue;
      else PyList_Append (list, PyString_FromString (argv[i]));
    }
  }
  else
  {
    for (int i = 0; i < argc; i ++)
    {
      if (endswith (argv[i], "parmec4")) continue;
      else if (endswith (argv[i], "parmec8")) continue;
      else PyList_Append (list, PyString_FromString (argv[i]));
    }
  }

  return list;
}

/* reset simulation */
static PyObject* RESET (PyObject *self, PyObject *args, PyObject *kwds)
{
  ispc::master_free (master, parnum);

  ispc::slave_free (slave, parnum);

  master = ispc::master_alloc (NULL, 0, particle_buffer_size);

  slave = ispc::slave_alloc (NULL, 0, particle_buffer_size);

  reset ();

  Py_RETURN_NONE;
}

/* create time series */
static PyObject* TSERIES (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("points");
  PyObject *points;
  TMS *ts;

  PARSEKEYS ("O", &points);

  TYPETEST (is_number_or_list_or_string (points, kwl [0], 1, 2));

  if (PyNumber_Check (points))
  {
    ts = TMS_Constant (PyFloat_AsDouble (points));
  }
  else if (PyString_Check (points))
  {
    if (!(ts = TMS_File (PyString_AsString (points))))
    {
      PyErr_SetString (PyExc_ValueError, "Could not open file");
      return NULL;
    }
  }
  else if (PyList_Check (points))
  {
    REAL *times, *values;
    int i, n;

    if (PyList_Check (PyList_GetItem (points, 0)))
    {
      n = PyList_Size (points);

      ERRMEM (times = (REAL*)malloc (sizeof (REAL [n])));
      ERRMEM (values = (REAL*)malloc (sizeof (REAL [n])));

      for (i = 0; i < n; i ++)
      {
	PyObject *pv = PyList_GetItem (points, i);

	TYPETEST (is_list (pv, "[t,v]", 2));

	times [i] = PyFloat_AsDouble (PyList_GetItem (pv, 0));
	values [i] = PyFloat_AsDouble (PyList_GetItem (pv, 1));
      }
    }
    else if (PyTuple_Check (PyList_GetItem (points, 0)))
    {
      n = PyList_Size (points);

      ERRMEM (times = (REAL*)malloc (sizeof (REAL [n])));
      ERRMEM (values = (REAL*)malloc (sizeof (REAL [n])));

      for (i = 0; i < n; i ++)
      {
	PyObject *pv = PyList_GetItem (points, i);

	TYPETEST (is_tuple (pv, "(t,v)", 2));

	times [i] = PyFloat_AsDouble (PyTuple_GetItem (pv, 0));
	values [i] = PyFloat_AsDouble (PyTuple_GetItem (pv, 1));
      }
    }
    else
    {
      if (PyList_Size(points) < 4)
      {
	PyErr_SetString (PyExc_ValueError, "Time series must have at least two points");
	return NULL;
      }

      n = PyList_Size (points) / 2;

      ERRMEM (times = (REAL*)malloc (sizeof (REAL [n])));
      ERRMEM (values = (REAL*)malloc (sizeof (REAL [n])));

      for (i = 0; i < n; i ++)
      {
	times [i] = PyFloat_AsDouble (PyList_GetItem (points, 2*i));
	values [i] = PyFloat_AsDouble (PyList_GetItem (points, 2*i + 1));
      }
    }

    ts = TMS_Create (n, times, values);

    free (times);
    free (values);
  }
  else
  {
    PyErr_SetString (PyExc_TypeError, "Invalid points format");
    return NULL;
  }

  if (tmsnum >= time_series_buffer_size) time_series_buffer_grow ();

  int i = tmsnum ++;

  parmec::tms[i] = ts;

  return PyInt_FromLong (i);
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
            is_ge_lt (material, 0, matnum, kwl[2]) && is_positive (color, kwl[3]));

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
  invm[i] = 1.0/mass[i];

  rotation[0][i] = rotation[4][i] = rotation[8][i] = 1.0;
  rotation[1][i] = rotation[2][i] = rotation[3][i] =
  rotation[5][i] = rotation[6][i] = rotation[7][i] = 0.0;

  inertia[0][i] = inertia[4][i] = inertia[8][i] = 0.4*mass[i]*radii[0][j]*radii[0][j];
  inertia[1][i] = inertia[2][i] = inertia[3][i] =
  inertia[5][i] = inertia[6][i] = inertia[7][i] = 0.0;

  REAL J[9] = {inertia[0][i], inertia[1][i], inertia[2][i],
               inertia[3][i], inertia[4][i], inertia[5][i],
	       inertia[6][i], inertia[7][i], inertia[8][i]}, Jiv[9], det;
  INVERT (J, Jiv, det);
  inverse[0][i] = Jiv[0];
  inverse[1][i] = Jiv[1];
  inverse[2][i] = Jiv[2];
  inverse[3][i] = Jiv[3];
  inverse[4][i] = Jiv[4];
  inverse[5][i] = Jiv[5];
  inverse[6][i] = Jiv[6];
  inverse[7][i] = Jiv[7];
  inverse[8][i] = Jiv[8];

  /* return regular particle */
  flags[i] = OUTREST;

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

  PARSEKEYS ("OOiO", &nodes, &elements, &material, &colors);

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
    nodpart[k] = parnum;
  }
  nodnum += node_count;

  /* insert particle data */

  if (parnum >= particle_buffer_size) particle_buffer_grow ();

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
  invm[i] = 1.0/mi;

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

  REAL J[9] = {inertia[0][i], inertia[1][i], inertia[2][i],
               inertia[3][i], inertia[4][i], inertia[5][i],
	       inertia[6][i], inertia[7][i], inertia[8][i]}, Jiv[9], det;
  INVERT (J, Jiv, det);
  inverse[0][i] = Jiv[0];
  inverse[1][i] = Jiv[1];
  inverse[2][i] = Jiv[2];
  inverse[3][i] = Jiv[3];
  inverse[4][i] = Jiv[4];
  inverse[5][i] = Jiv[5];
  inverse[6][i] = Jiv[6];
  inverse[7][i] = Jiv[7];
  inverse[8][i] = Jiv[8];

  /* clean up */
  MESH_Destroy (msh);
  free (lnod);
  free (lele);
  free (lsur);

  /* return regular particle */
  flags[i] = OUTREST;

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

  PARSEKEYS ("|OdOOii", &inertia, &mass, &rotation, &position, &material, &particle);

  TYPETEST (is_list (inertia, kwl[0], 6) && is_non_negative (mass, kwl[1]) &&
            is_list (rotation, kwl[2], 9) && is_tuple (position, kwl[3], 3) &&
	    is_non_negative (material, kwl[4]));

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
    if (parnum >= particle_buffer_size) particle_buffer_grow ();

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

    parmec::rotation[0][i] = vrotation[0];
    parmec::rotation[1][i] = vrotation[1];
    parmec::rotation[2][i] = vrotation[2]; 
    parmec::rotation[3][i] = vrotation[3];
    parmec::rotation[4][i] = vrotation[4];
    parmec::rotation[5][i] = vrotation[5];
    parmec::rotation[6][i] = vrotation[6];
    parmec::rotation[7][i] = vrotation[7];
    parmec::rotation[8][i] = vrotation[8];

    parmec::mass[i] = (mass == 0.0 ? 1.0 : mass);
    parmec::invm[i] = 1.0/parmec::mass[i];

    parmec::inertia[0][i] = vinertia[0];
    parmec::inertia[4][i] = vinertia[1];
    parmec::inertia[8][i] = vinertia[2];
    parmec::inertia[1][i] = parmec::inertia[3][i] = vinertia[3];
    parmec::inertia[2][i] = parmec::inertia[6][i] = vinertia[4];
    parmec::inertia[5][i] = parmec::inertia[7][i] = vinertia[5];

    REAL J[9] = {parmec::inertia[0][i], parmec::inertia[1][i], parmec::inertia[2][i],
		 parmec::inertia[3][i], parmec::inertia[4][i], parmec::inertia[5][i],
		 parmec::inertia[6][i], parmec::inertia[7][i], parmec::inertia[8][i]}, Jiv[9], det;
    INVERT (J, Jiv, det);
    parmec::inverse[0][i] = Jiv[0];
    parmec::inverse[1][i] = Jiv[1];
    parmec::inverse[2][i] = Jiv[2];
    parmec::inverse[3][i] = Jiv[3];
    parmec::inverse[4][i] = Jiv[4];
    parmec::inverse[5][i] = Jiv[5];
    parmec::inverse[6][i] = Jiv[6];
    parmec::inverse[7][i] = Jiv[7];
    parmec::inverse[8][i] = Jiv[8];
  }
  else /* redefine an existing particle as analytical */
  {
    if (particle >= parnum)
    {
      PyErr_SetString (PyExc_ValueError, "Particle index is out of range");
      return NULL;
    }

    declare_analytical (particle);

    i = particle;

    if (position)
    {
      parmec::position[0][i] = vposition[0];
      parmec::position[1][i] = vposition[1];
      parmec::position[2][i] = vposition[2];
      parmec::position[3][i] = vposition[0];
      parmec::position[4][i] = vposition[1];
      parmec::position[5][i] = vposition[2];
    }

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

    if (mass > 0.0)
    {
      parmec::mass[i] = mass;
      parmec::invm[i] = 1.0/parmec::mass[i];
    }

    if (inertia)
    {
      parmec::inertia[0][i] = vinertia[0];
      parmec::inertia[4][i] = vinertia[1];
      parmec::inertia[8][i] = vinertia[2];
      parmec::inertia[1][i] = parmec::inertia[3][i] = vinertia[3];
      parmec::inertia[2][i] = parmec::inertia[6][i] = vinertia[4];
      parmec::inertia[5][i] = parmec::inertia[7][i] = vinertia[5];

      REAL J[9] = {parmec::inertia[0][i], parmec::inertia[1][i], parmec::inertia[2][i],
		   parmec::inertia[3][i], parmec::inertia[4][i], parmec::inertia[5][i],
		   parmec::inertia[6][i], parmec::inertia[7][i], parmec::inertia[8][i]}, Jiv[9], det;
      INVERT (J, Jiv, det);
      parmec::inverse[0][i] = Jiv[0];
      parmec::inverse[1][i] = Jiv[1];
      parmec::inverse[2][i] = Jiv[2];
      parmec::inverse[3][i] = Jiv[3];
      parmec::inverse[4][i] = Jiv[4];
      parmec::inverse[5][i] = Jiv[5];
      parmec::inverse[6][i] = Jiv[6];
      parmec::inverse[7][i] = Jiv[7];
      parmec::inverse[8][i] = Jiv[8];
    }
  }

  /* return analytical particle */
  flags[particle] = parmec::ANALYTICAL|OUTREST;

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
  KEYWORDS ("part1", "point1", "part2", "geom2", "spring", "dashpot", "direction", "planar", "unload", "ylim", "inactive", "offset");
  PyObject *point1, *geom2, *spring, *dashpot, *direction, *planar, *unload, *ylim, *inactive;
  int part1, part2, offset;

  direction = NULL;
  dashpot = NULL;
  planar = NULL;
  unload = NULL;
  ylim = NULL;
  inactive = Py_False;
  offset = -1; /* unspecified spring stroke offest */

  PARSEKEYS ("iOiOO|OOOOOOi", &part1, &point1, &part2, &geom2, &spring, &dashpot, &direction, &planar, &unload, &ylim, &inactive, &offset);

  TYPETEST (is_non_negative (part1, kwl[0]) && is_tuple (point1, kwl[1], 3) &&
            is_tuple_or_list_of_tuples (geom2, kwl[3], 3, 2) &&
	    is_list (spring, kwl[4], 0) && is_list_or_number (dashpot, kwl[5], 0) &&
	    is_tuple (direction, kwl[6], 3) && is_string (planar, kwl[7]) &&
	    is_list (unload, kwl[8], 0) && is_tuple (ylim, kwl[9], 2) && is_bool (inactive, kwl[10]));

  if (offset < -1 || offset >= tmsnum)
  {
    PyErr_SetString (PyExc_ValueError, "offset TSERIES index is out of range");
    return NULL;
  }

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

  if (dashpot && PyList_Check (dashpot) && (PyList_Size (dashpot) < 4 || PyList_Size (dashpot) % 2))
  {
    PyErr_SetString (PyExc_ValueError, "Invalid dashpot lookup table list length");
    return NULL;
  }

  if (unload && (PyList_Size (unload) < 4 || PyList_Size (unload) % 2))
  {
    PyErr_SetString (PyExc_ValueError, "Invalid unload lookup table list length");
    return NULL;
  }

  int spring_lookup = PyList_Size (spring);

  int dashpot_lookup = dashpot &&  PyList_Check(dashpot) ? PyList_Size (dashpot) : 4;

  int unload_lookup = unload ? PyList_Size (unload) : 0;

  spring_buffer_grow (spring_lookup, dashpot_lookup, unload_lookup);

  int i = sprnum ++;

  sprid[i] = sprmap[i] = i;

  springs_changed = 1;

  sprpart[0][i] = part1;
  sprpart[1][i] = part2;

  sprpnt[0][0][i] = PyFloat_AsDouble (PyTuple_GetItem (point1,0));
  sprpnt[0][1][i] = PyFloat_AsDouble (PyTuple_GetItem (point1,1));
  sprpnt[0][2][i] = PyFloat_AsDouble (PyTuple_GetItem (point1,2));
  sprpnt[0][3][i] = sprpnt[0][0][i];
  sprpnt[0][4][i] = sprpnt[0][1][i];
  sprpnt[0][5][i] = sprpnt[0][2][i];

  if (PyTuple_Check(geom2))
  {
    sprpnt[1][0][i] = PyFloat_AsDouble (PyTuple_GetItem (geom2,0));
    sprpnt[1][1][i] = PyFloat_AsDouble (PyTuple_GetItem (geom2,1));
    sprpnt[1][2][i] = PyFloat_AsDouble (PyTuple_GetItem (geom2,2));
    sprpnt[1][3][i] = sprpnt[1][0][i];
    sprpnt[1][4][i] = sprpnt[1][1][i];
    sprpnt[1][5][i] = sprpnt[1][2][i];
  }
  else
  {
    PyObject *point = PyList_GetItem (geom2, 0);
    PyObject *normal = PyList_GetItem (geom2, 1);

    sprpnt[1][0][i] = PyFloat_AsDouble (PyTuple_GetItem (point,0));
    sprpnt[1][1][i] = PyFloat_AsDouble (PyTuple_GetItem (point,1));
    sprpnt[1][2][i] = PyFloat_AsDouble (PyTuple_GetItem (point,2));
    sprpnt[1][3][i] = sprpnt[1][0][i];
    sprpnt[1][4][i] = sprpnt[1][1][i];
    sprpnt[1][5][i] = sprpnt[1][2][i];

    REAL dir[3] = {PyFloat_AsDouble (PyTuple_GetItem (normal,0)),
                   PyFloat_AsDouble (PyTuple_GetItem (normal,1)),
                   PyFloat_AsDouble (PyTuple_GetItem (normal,2))};

    NORMALIZE (dir);

    sprdir[0][i] = dir[0];
    sprdir[1][i] = dir[1];
    sprdir[2][i] = dir[2];
    sprdir[3][i] = dir[0];
    sprdir[4][i] = dir[1];
    sprdir[5][i] = dir[2];

    REAL p[3] = {sprpnt[0][0][i], sprpnt[0][1][i], sprpnt[0][2][i]};
    REAL q[3] = {sprpnt[1][0][i], sprpnt[1][1][i], sprpnt[1][2][i]};
    REAL w[3], len;
    SUB (p, q, w);
    len = DOT (dir, w);
    SUBMUL (p, len, dir, q); /* q = projection of p onto the spring plane */

    sprpnt[1][0][i] = q[0];
    sprpnt[1][1][i] = q[1];
    sprpnt[1][2][i] = q[2];
    sprpnt[1][3][i] = q[0];
    sprpnt[1][4][i] = q[1];
    sprpnt[1][5][i] = q[2];
  }

  sprflg[i] = 0;
  sproffset[i] = offset < 0 ? offset : lcurve_from_time_series (offset);

  int j = 0, k = spridx[i];

  for (; j < spring_lookup/2; j ++, k ++)
  {
    REAL stroke = PyFloat_AsDouble(PyList_GetItem(spring,2*j));
    REAL force = PyFloat_AsDouble(PyList_GetItem(spring,2*j+1));
    parmec::spring[0][k] = stroke;
    parmec::spring[1][k] = force;
  }
  spridx[sprnum] = k;

  if (dashpot && PyList_Check(dashpot))
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
  else if (PyNumber_Check(dashpot)) /* critical damping ratio */
  {
    REAL ratio = PyFloat_AsDouble (dashpot);
    if (ratio < 0.0)
    {
      PyErr_SetString (PyExc_ValueError, "Critical damping ratio not in [0.0, +Inf) interval");
      return NULL;
    }
    k = dashidx[i];
    parmec::dashpot[0][k] = ratio;
    parmec::dashpot[1][k] = ratio;
    dashidx[sprnum] = k+1;
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

  if (unload)
  {
    for (j = 0, k = unidx[i]; j < unload_lookup/2; j ++, k ++)
    {
      REAL stroke = PyFloat_AsDouble(PyList_GetItem(unload,2*j));
      REAL force = PyFloat_AsDouble(PyList_GetItem(unload,2*j+1));
      parmec::unload[0][k] = stroke;
      parmec::unload[1][k] = force;

      if (j)
      {
	REAL slope = (force-parmec::unload[1][k-1])/(stroke-parmec::unload[0][k-1]);
	if (slope <= 0.0) /* unloading curve must be monotonic --> (x,y) slope search */
	{
	  PyErr_SetString (PyExc_ValueError, "Unloading curve must be monotonically increasing");
	  return NULL;
	}
      }
    }
    unidx[sprnum] = k;
    parmec::sprtype[i] = SPRING_GENERAL_NONLINEAR;
  }
  else
  {
    parmec::unidx[sprnum] = parmec::unidx[sprnum-1];
    parmec::sprtype[i] = SPRING_NONLINEAR_ELASTIC;
  }

  if (inactive == Py_True)
  {
    parmec::unspring[i] = -1; /* zero force */
  }
  else
  {
    parmec::unspring[i] = -3; /* unsued by UNSPRING */
  }

  if (PyList_Check (geom2)) /* project point1 onto plane */
  {
    sprflg[i] |= SPRDIR_PROJECT;

    REAL dif[3] = {sprpnt[0][0][i] - sprpnt[1][0][i],
		   sprpnt[0][1][i] - sprpnt[1][1][i],
		   sprpnt[0][2][i] - sprpnt[1][2][i]};

    REAL nor[3] = {sprdir[0][i], sprdir[1][i], sprdir[2][i]};

    stroke0[i] = 0.0; /* point-to-plane contact type spring - no initial stroke in this case */
  }
  else
  {
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
      sprdir[3][i] = dir[0];
      sprdir[4][i] = dir[1];
      sprdir[5][i] = dir[2];

      if (planar)
      {
	IFIS (planar, "ON") /* direction in plane */
	{
	  sprflg[i] |= SPRDIR_PLANAR;
	}
	ELIF (planar, "OFF") /* constant direction */
	{
	  sprflg[i] |= SPRDIR_CONSTANT;
	}
	ELSE
	{
	  PyErr_SetString (PyExc_ValueError, "Invalid planar switch");
	  return NULL;
	}
      }
      else /* constant direction */
      {
	sprflg[i] |= SPRDIR_CONSTANT;
      }

      REAL dif[3] = {sprpnt[1][0][i] - sprpnt[0][0][i],
		     sprpnt[1][1][i] - sprpnt[0][1][i],
		     sprpnt[1][2][i] - sprpnt[0][2][i]};

      if (sprflg[i] & SPRDIR_CONSTANT)
      {
	stroke0[i] = DOT (dif, dir); /* stroke(0) = projection along direction */
      }
      else /* stroke(0) = length in orthogonal plane */
      {
	REAL dot = DOT (dif, dir);

	dif[0] -= dot*dir[0];
	dif[1] -= dot*dir[1];
	dif[2] -= dot*dir[2];

	stroke0[i] = LEN (dif);
      }
    }
    else /* direction = (p2 - p1)/|p2 - p1| */
    {
      sprflg[i] |= SPRDIR_FOLLOWER;

      REAL dif[3] = {sprpnt[1][0][i] - sprpnt[0][0][i],
		     sprpnt[1][1][i] - sprpnt[0][1][i],
		     sprpnt[1][2][i] - sprpnt[0][2][i]};

      stroke0[i] = LEN (dif);
    }
  }

  if (ylim)
  {
    yield[0][i] = PyFloat_AsDouble (PyTuple_GetItem (ylim,0));
    yield[1][i] = PyFloat_AsDouble (PyTuple_GetItem (ylim,1));

    if (yield[0][i] > 0.0)
    {
      PyErr_SetString (PyExc_ValueError, "Compressive yield limit must be <= 0.0");
      return NULL;
    }

    if (yield[1][i] < 0.0)
    {
      PyErr_SetString (PyExc_ValueError, "Tensile yield limit must be >= 0.0");
      return NULL;
    }
  }
  else
  {
    parmec::yield[0][i] = 0.0;
    parmec::yield[1][i] = 0.0;
  }

  return PyInt_FromLong (i);
}

/* add collective spring undoing condition */
static PyObject* UNSPRING (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("tsprings", "msprings", "limits", "entity", "operator", "abs", "nsteps", "nfreq", "unload", "activate");
  PyObject *tsprings, *msprings, *limits, *entity, *operat, *abs, *activate;
  int nsteps, nfreq, unload;

  entity = NULL;
  operat = NULL;
  abs = Py_False;
  activate = NULL;
  nsteps = 1;
  nfreq = 1;
  unload = -1;

  PARSEKEYS ("OOO|OOOiiOO", &tsprings, &msprings, &limits, &entity, &operat, &abs, &nsteps, &nfreq, &unload, &activate);

  TYPETEST (is_list (tsprings, kwl[0], 0) && is_list (msprings, kwl[1], 0) && is_tuple (limits, kwl[2], 2) &&
            is_string (entity, kwl[3]) && is_string (operat, kwl[4]) && is_bool (abs, kwl[5]) &&
	    is_positive (nsteps, kwl[6]) && is_positive (nfreq, kwl[7]) && is_list (activate, kwl[9], 0));

  if (unload < -1 || unload >= tmsnum)
  {
    PyErr_SetString (PyExc_ValueError, "unload TSERIES index is out of range");
    return NULL;
  }
  else if (unload >= 0)
  {
    TMS *ts = (TMS*)tms[unload];

    for (int j = 0; j < ts->size; j ++)
    {
      REAL stroke = ts->points[0][j];
      REAL force = ts->points[1][j];

      if (j)
      {
	REAL slope = (force - ts->points[1][j-1])/(stroke - ts->points[0][j-1]);
	if (slope <= 0.0) /* unloading curve must be monotonic --> (x,y) slope search */
	{
	  PyErr_SetString (PyExc_ValueError, "Unloading curve must be monotonically increasing");
	  return NULL;
	}
      }
    }
  }

  for (int i = 0; i < PyList_Size (tsprings); i ++)
  {
    int j = PyInt_AsLong(PyList_GetItem(tsprings, i));

    if (j < 0 || j >= sprnum)
    {
      PyErr_SetString (PyExc_ValueError, "tsprings SPRING index is out of range");
      return NULL;
    }
  }

  for (int i = 0; i < PyList_Size (msprings); i ++)
  {
    int j = PyInt_AsLong(PyList_GetItem(msprings, i));

    if (j < 0 || j >= sprnum)
    {
      PyErr_SetString (PyExc_ValueError, "msprings SPRING index is out of range");
      return NULL;
    }
    else if (parmec::unspring[j] == -2) /* spring index already used in msprings by an earlier UNSPRING call */
    {
      char message [1024];
      snprintf (message, 1024, "in msprings, spring index %d already in use by an earlier UNSPRING call", j);
      PyErr_SetString (PyExc_ValueError, message);
      return NULL;
    }
    else if (parmec::unspring[j] == -1) /* inactctive SPRING */
    {
      char message [1024];
      snprintf (message, 1024, "in msprings, trying to use an inactive SPRING with index %d", j);
      PyErr_SetString (PyExc_ValueError, message);
      return NULL;
    }
  }

  int actsize = activate ? PyList_Size (activate) : 0;

  unspring_buffer_grow (PyList_Size (tsprings), PyList_Size (msprings), actsize);

  int i = unsprnum ++;

  tspridx[i+1] = tspridx[i];
  for (int k = 0; k < PyList_Size (tsprings); k ++)
  {
    int j = PyInt_AsLong(PyList_GetItem(tsprings, k));
    parmec::tsprings[tspridx[i+1]] = j;
    tspridx[i+1] ++;
  }

  mspridx[i+1] = mspridx[i];
  for (int k = 0; k < PyList_Size (msprings); k ++)
  {
    int j = PyInt_AsLong(PyList_GetItem(msprings, k));
    parmec::msprings[mspridx[i+1]] = j;
    parmec::unspring[j] = -2; /* mark as used by UNSPRING */
    mspridx[i+1] ++;
  }

  PyObject *lim0 = PyTuple_GetItem (limits, 0);

  if (lim0 == Py_None) unlim[0][i] = -REAL_MAX;
  else if (PyNumber_Check (lim0))
  {
    unlim[0][i] = PyFloat_AsDouble (lim0);
  }
  else 
  {
    PyErr_SetString (PyExc_ValueError, "Invalid lower limit");
    return NULL;
  }

  PyObject *lim1 = PyTuple_GetItem (limits, 1);

  if (lim1 == Py_None) unlim[1][i] = REAL_MAX;
  else if (PyNumber_Check (lim1))
  {
    unlim[1][i] = PyFloat_AsDouble (lim1);
  }
  else 
  {
    PyErr_SetString (PyExc_ValueError, "Invalid upper limit");
    return NULL;
  }

  if (unlim[0][i] >= unlim[1][i])
  {
    PyErr_SetString (PyExc_ValueError, "Lower limit >= upper limit");
    return NULL;
  }

  if (entity)
  {
    IFIS (entity, "STROKE")
    {
      unent[i] = HIS_STROKE;
    }
    ELIF (entity, "STF") /* total force */
    {
      unent[i] = HIS_STF;
    }
    ELIF (entity, "SF") /* elastic force */
    {
      unent[i] = HIS_SF;
    }
    ELSE
    {
      PyErr_SetString (PyExc_ValueError, "Invalid entity");
      return NULL;
    }
  }
  else /* default */
  {
    unent[i] = HIS_SF;
  }

  if (operat)
  {
    IFIS (operat, "SUM")
    {
      unop[i] = OP_SUM;
    }
    ELIF (operat, "MIN")
    {
      unop[i] = OP_MIN;
    }
    ELIF (operat, "MAX")
    {
      unop[i] = OP_MAX;
    }
    ELSE
    {
      PyErr_SetString (PyExc_ValueError, "Invalid operator");
      return NULL;
    }
  }
  else /* default */
  {
    unent[i] = HIS_SF;
  }

  unabs[i] = abs == Py_True ? 1 : 0;
  parmec::nsteps[i] = nsteps;
  parmec::nfreq[i] = nfreq;
  unaction[i] = unload < 0 ? unload : lcurve_from_time_series (unload);

  actidx[i+1] = actidx[i];
  for (int k = 0; activate && k < PyList_Size (activate); k ++)
  {
    int j = PyInt_AsLong(PyList_GetItem(activate, k));
    if (parmec::unspring[j] != -1) /* test if inactive */
    {
      char buff[1024];
      snprintf (buff, 1024, "in activate list: SPRING %d is already active", j);
      PyErr_SetString (PyExc_ValueError, buff);
      return NULL;
    }
    parmec::activate[actidx[i+1]] = j;
    actidx[i+1] ++;
  }

  return PyInt_FromLong (i);
}

/* Calculate equivalent point mass from particle inertia and mass properties */
static PyObject* EQM (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("part1", "point1", "part2", "point2", "direction");
  PyObject *point1, *point2, *direction;
  int part1, part2;

  part2 = -1;
  point2 = NULL;
  direction = NULL;

  PARSEKEYS ("iO|iOO", &part1, &point1, &part2, &point2, &direction);

  TYPETEST (is_ge_lt (part1, 0, parnum, kwl[0]) && is_tuple (point1, kwl[1], 3) &&
            is_tuple (point2, kwl[3], 3) && is_tuple (direction, kwl[4], 3));

  REAL p1[3], p2[3], d[3], *pd = NULL;

  p1[0] = PyFloat_AsDouble(PyTuple_GetItem (point1, 0));
  p1[1] = PyFloat_AsDouble(PyTuple_GetItem (point1, 1));
  p1[2] = PyFloat_AsDouble(PyTuple_GetItem (point1, 2));

  if (part2 >= 0)
  {
    if (part2 >= parnum)
    {
      PyErr_SetString (PyExc_ValueError, "Particle number 'part2' is out of range");
      return NULL;
    }

    p2[0] = PyFloat_AsDouble(PyTuple_GetItem (point2, 0));
    p2[1] = PyFloat_AsDouble(PyTuple_GetItem (point2, 1));
    p2[2] = PyFloat_AsDouble(PyTuple_GetItem (point2, 2));
  }

  if (direction)
  {
    d[0] = PyFloat_AsDouble(PyTuple_GetItem (direction, 0));
    d[1] = PyFloat_AsDouble(PyTuple_GetItem (direction, 1));
    d[2] = PyFloat_AsDouble(PyTuple_GetItem (direction, 2));
    pd = d;
  }

  REAL mass = ispc::eqm (part1, p1, part2, p2, pd, parmec::inverse, parmec::invm, parmec::position);

  return PyFloat_FromDouble (mass);
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
    PyErr_SetString (PyExc_ValueError, "Invalid color values");
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

/* restrain particle  motion */
static PyObject* RESTRAIN (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("parnum", "linear", "angular");
  PyObject *lin, *ang;
  REAL dir[3][3], dot;
  int i, j;

  lin = NULL;
  ang = NULL;

  PARSEKEYS ("i|OO", &j, &lin, &ang);

  TYPETEST (is_ge_lt (j, 0, parnum, kwl[0]) && is_list (lin, kwl[1], 0) && is_list (ang, kwl[2], 0));

  if (lin || ang)
  {
    if (cnsnum >= constrain_buffer_size) constrain_buffer_grow ();

    i = cnsnum ++;

    cnspart[i] = j;
  }

  if (lin)
  {
    if (PyList_Size (lin) == 3)
    {
      dir[0][0] = PyFloat_AsDouble (PyList_GetItem (lin, 0)); 
      dir[0][1] = PyFloat_AsDouble (PyList_GetItem (lin, 1)); 
      dir[0][2] = PyFloat_AsDouble (PyList_GetItem (lin, 2)); 
      dir[1][0] = 0.0;
      dir[1][1] = 0.0;
      dir[1][2] = 0.0;
      dir[2][0] = 0.0;
      dir[2][1] = 0.0;
      dir[2][2] = 0.0;
      NORMALIZE (dir[0]);
    }
    else if (PyList_Size (lin) == 6)
    {
      dir[0][0] = PyFloat_AsDouble (PyList_GetItem (lin, 0)); 
      dir[0][1] = PyFloat_AsDouble (PyList_GetItem (lin, 1)); 
      dir[0][2] = PyFloat_AsDouble (PyList_GetItem (lin, 2)); 
      dir[1][0] = PyFloat_AsDouble (PyList_GetItem (lin, 3));  
      dir[1][1] = PyFloat_AsDouble (PyList_GetItem (lin, 4)); 
      dir[1][2] = PyFloat_AsDouble (PyList_GetItem (lin, 5)); 
      dir[2][0] = 0.0;
      dir[2][1] = 0.0;
      dir[2][2] = 0.0;
      NORMALIZE (dir[0]);
      dot = DOT (dir[0], dir[1]);
      SUBMUL (dir[1], dot, dir[0], dir[1]);
      NORMALIZE (dir[1]);
    }
    else if (PyList_Size (lin) == 9)
    {
      dir[0][0] = PyFloat_AsDouble (PyList_GetItem (lin, 0)); 
      dir[0][1] = PyFloat_AsDouble (PyList_GetItem (lin, 1)); 
      dir[0][2] = PyFloat_AsDouble (PyList_GetItem (lin, 2)); 
      dir[1][0] = PyFloat_AsDouble (PyList_GetItem (lin, 3));  
      dir[1][1] = PyFloat_AsDouble (PyList_GetItem (lin, 4)); 
      dir[1][2] = PyFloat_AsDouble (PyList_GetItem (lin, 5)); 
      dir[2][0] = PyFloat_AsDouble (PyList_GetItem (lin, 6));  
      dir[2][1] = PyFloat_AsDouble (PyList_GetItem (lin, 7));  
      dir[2][2] = PyFloat_AsDouble (PyList_GetItem (lin, 8));  
      NORMALIZE (dir[0]);
      dot = DOT (dir[0], dir[1]);
      SUBMUL (dir[1], dot, dir[0], dir[1]);
      NORMALIZE (dir[1]);
      dot = DOT (dir[0], dir[2]);
      SUBMUL (dir[2], dot, dir[0], dir[2]);
      dot = DOT (dir[1], dir[2]);
      SUBMUL (dir[2], dot, dir[1], dir[2]);
      NORMALIZE (dir[2]);
    }
    else
    {
      PyErr_SetString (PyExc_ValueError, "Invalid linear constrained directions list");
      return NULL;
    }

    cnslin[0][i] = dir[0][0];
    cnslin[1][i] = dir[0][1];
    cnslin[2][i] = dir[0][2];
    cnslin[3][i] = dir[1][0];
    cnslin[4][i] = dir[1][1];
    cnslin[5][i] = dir[1][2];
    cnslin[6][i] = dir[2][0];
    cnslin[7][i] = dir[2][1];
    cnslin[8][i] = dir[2][2];
  }
  else if (ang)
  {
    cnslin[0][i] = 0;
    cnslin[1][i] = 0;
    cnslin[2][i] = 0;
    cnslin[3][i] = 0;
    cnslin[4][i] = 0;
    cnslin[5][i] = 0;
    cnslin[6][i] = 0;
    cnslin[7][i] = 0;
    cnslin[8][i] = 0;
  }

  if (ang)
  {
    if (PyList_Size (ang) == 3)
    {
      dir[0][0] = PyFloat_AsDouble (PyList_GetItem (ang, 0)); 
      dir[0][1] = PyFloat_AsDouble (PyList_GetItem (ang, 1)); 
      dir[0][2] = PyFloat_AsDouble (PyList_GetItem (ang, 2)); 
      dir[1][0] = 0.0;
      dir[1][1] = 0.0;
      dir[1][2] = 0.0;
      dir[2][0] = 0.0;
      dir[2][1] = 0.0;
      dir[2][2] = 0.0;
      NORMALIZE (dir[0]);
    }
    else if (PyList_Size (ang) == 6)
    {
      dir[0][0] = PyFloat_AsDouble (PyList_GetItem (ang, 0)); 
      dir[0][1] = PyFloat_AsDouble (PyList_GetItem (ang, 1)); 
      dir[0][2] = PyFloat_AsDouble (PyList_GetItem (ang, 2)); 
      dir[1][0] = PyFloat_AsDouble (PyList_GetItem (ang, 3));  
      dir[1][1] = PyFloat_AsDouble (PyList_GetItem (ang, 4)); 
      dir[1][2] = PyFloat_AsDouble (PyList_GetItem (ang, 5)); 
      dir[2][0] = 0.0;
      dir[2][1] = 0.0;
      dir[2][2] = 0.0;
      NORMALIZE (dir[0]);
      dot = DOT (dir[0], dir[1]);
      SUBMUL (dir[1], dot, dir[0], dir[1]);
      NORMALIZE (dir[1]);
    }
    else if (PyList_Size (ang) == 9)
    {
      dir[0][0] = PyFloat_AsDouble (PyList_GetItem (ang, 0)); 
      dir[0][1] = PyFloat_AsDouble (PyList_GetItem (ang, 1)); 
      dir[0][2] = PyFloat_AsDouble (PyList_GetItem (ang, 2)); 
      dir[1][0] = PyFloat_AsDouble (PyList_GetItem (ang, 3));  
      dir[1][1] = PyFloat_AsDouble (PyList_GetItem (ang, 4)); 
      dir[1][2] = PyFloat_AsDouble (PyList_GetItem (ang, 5)); 
      dir[2][0] = PyFloat_AsDouble (PyList_GetItem (ang, 6));  
      dir[2][1] = PyFloat_AsDouble (PyList_GetItem (ang, 7));  
      dir[2][2] = PyFloat_AsDouble (PyList_GetItem (ang, 8));  
      NORMALIZE (dir[0]);
      dot = DOT (dir[0], dir[1]);
      SUBMUL (dir[1], dot, dir[0], dir[1]);
      NORMALIZE (dir[1]);
      dot = DOT (dir[0], dir[2]);
      SUBMUL (dir[2], dot, dir[0], dir[2]);
      dot = DOT (dir[1], dir[2]);
      SUBMUL (dir[2], dot, dir[1], dir[2]);
      NORMALIZE (dir[2]);
    }
    else
    {
      PyErr_SetString (PyExc_ValueError, "Invalid angular constrained directions list");
      return NULL;
    }

    cnsang[0][i] = dir[0][0];
    cnsang[1][i] = dir[0][1];
    cnsang[2][i] = dir[0][2];
    cnsang[3][i] = dir[1][0];
    cnsang[4][i] = dir[1][1];
    cnsang[5][i] = dir[1][2];
    cnsang[6][i] = dir[2][0];
    cnsang[7][i] = dir[2][1];
    cnsang[8][i] = dir[2][2];
  }
  else if (lin)
  {
    cnsang[0][i] = 0;
    cnsang[1][i] = 0;
    cnsang[2][i] = 0;
    cnsang[3][i] = 0;
    cnsang[4][i] = 0;
    cnsang[5][i] = 0;
    cnsang[6][i] = 0;
    cnsang[7][i] = 0;
    cnsang[8][i] = 0;
  }

  Py_RETURN_NONE;
}

/* prescribe particle moion */
static PyObject* PRESCRIBE (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("parnum", "linear", "angular", "kind");
  PyObject *lin, *ang, *kind;
  char *kindvalue;
  int i, j;

  lin = NULL;
  ang = NULL;
  kind = NULL;

  PARSEKEYS ("i|OOO", &j, &lin, &ang, &kind);

  TYPETEST (is_ge_lt (j, 0, parnum, kwl[0]) && is_string (kind, kwl[3]));

  if (lin || ang)
  {
    if (prsnum >= prescribe_buffer_size) prescribe_buffer_grow ();

    i = prsnum ++;

    prspart[i] = j;
  }

  if (lin)
  {
    if (PyCallable_Check(lin))
    {
      prslin[i] = lin;
      tmslin[0][i] = -1;
      tmslin[1][i] = -1;
      tmslin[2][i] = -1;
    }
    else if (PyTuple_Check(lin))
    {
      if (PyTuple_Size(lin) != 3) 
      {
	PyErr_SetString (PyExc_ValueError, "'linear' tuple size != 3");
	return NULL;
      }

      prslin[i] = NULL;
      tmslin[0][i] = PyInt_AsLong(PyTuple_GetItem(lin, 0));
      tmslin[1][i] = PyInt_AsLong(PyTuple_GetItem(lin, 1));
      tmslin[2][i] = PyInt_AsLong(PyTuple_GetItem(lin, 2));

      if (tmslin[0][i] < 0 || tmslin[0][i] >= tmsnum ||
          tmslin[1][i] < 0 || tmslin[1][i] >= tmsnum ||
          tmslin[2][i] < 0 || tmslin[2][i] >= tmsnum)
      {
	PyErr_SetString (PyExc_ValueError, "'linear' time series number out of range");
	return NULL;
      }
    }
    else
    {
      PyErr_SetString (PyExc_ValueError, "'linear' is neither a tuple nor a callback");
      return NULL;
    }
  }
  else
  {
    prslin[i] = NULL;
    tmslin[0][i] = -1;
    tmslin[1][i] = -1;
    tmslin[2][i] = -1;
  }

  if (ang)
  {
    if (PyCallable_Check(ang))
    {
      prsang[i] = ang;
      tmsang[0][i] = -1;
      tmsang[1][i] = -1;
      tmsang[2][i] = -1;
    }
    else if (PyTuple_Check(ang))
    {
      if (PyTuple_Size(ang) != 3) 
      {
	PyErr_SetString (PyExc_ValueError, "'angular' tuple size != 3");
	return NULL;
      }

      prsang[i] = NULL;
      tmsang[0][i] = PyInt_AsLong(PyTuple_GetItem(ang, 0));
      tmsang[1][i] = PyInt_AsLong(PyTuple_GetItem(ang, 1));
      tmsang[2][i] = PyInt_AsLong(PyTuple_GetItem(ang, 2));

      if (tmsang[0][i] < 0 || tmsang[0][i] >= tmsnum ||
          tmsang[1][i] < 0 || tmsang[1][i] >= tmsnum ||
          tmsang[2][i] < 0 || tmsang[2][i] >= tmsnum)
      {
	PyErr_SetString (PyExc_ValueError, "'angular' time series number out of range");
	return NULL;
      }
    }
    else
    {
      PyErr_SetString (PyExc_ValueError, "'angular' is neither a time series number nor a callback");
      return NULL;
    }
  }
  else
  {
    prsang[i] = NULL;
    tmsang[0][i] = -1;
    tmsang[1][i] = -1;
    tmsang[2][i] = -1;
  }

  if (kind)
  {
    IFIS (kind, "vv")
    {
      kindvalue = (char*)"vv";
    }
    ELIF (kind, "va")
    {
      kindvalue = (char*)"va";
    }
    ELIF (kind, "av")
    {
      kindvalue = (char*)"av";
    }
    ELIF (kind, "aa")
    {
      kindvalue = (char*)"aa";
    }
    ELSE
    {
      PyErr_SetString (PyExc_ValueError, "Invalid time history kind");
      return NULL;
    }
  }
  else
  {
    kindvalue = (char*)"vv";
  }

  if (lin || ang)
  {
    linkind[i] = kindvalue[0] == 'v' ? 0 : 1;
    angkind[i] = kindvalue[1] == 'v' ? 0 : 1;
  }

  Py_RETURN_NONE;
}

/* set particle velocity */
static PyObject* VELOCITY (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("parnum", "linear", "angular");
  PyObject *lin, *ang;
  int i;

  lin = NULL;
  ang = NULL;

  PARSEKEYS ("i|OO", &i, &lin, &ang);

  TYPETEST (is_ge_lt (i, 0, parnum, kwl[0]) && is_tuple (lin, kwl[1], 3) && is_tuple (ang, kwl[2], 3));

  if (lin)
  {
    linear[0][i] = PyFloat_AsDouble (PyTuple_GetItem (lin, 0)); 
    linear[1][i] = PyFloat_AsDouble (PyTuple_GetItem (lin, 1)); 
    linear[2][i] = PyFloat_AsDouble (PyTuple_GetItem (lin, 2)); 
  }

  if (ang)
  {
    angular[3][i] = PyFloat_AsDouble (PyTuple_GetItem (ang, 0)); 
    angular[4][i] = PyFloat_AsDouble (PyTuple_GetItem (ang, 1)); 
    angular[5][i] = PyFloat_AsDouble (PyTuple_GetItem (ang, 2)); 

    REAL o[3] = {angular[3][i], angular[4][i], angular[5][i]}, O[3];
    REAL L[9] = {rotation[0][i], rotation[1][i], rotation[2][i],
                 rotation[3][i], rotation[4][i], rotation[5][i],
		 rotation[6][i], rotation[7][i], rotation[8][i]};

    TVMUL (L,o,O);

    angular[0][i] = O[0];
    angular[1][i] = O[1];
    angular[2][i] = O[2];
  }

  Py_RETURN_NONE;
}

/* set gravity */
static PyObject* GRAVITY (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("gx", "gy", "gz");
  PyObject *gx, *gy, *gz;

  PARSEKEYS ("OOO", &gx, &gy, &gz);

  if (PyCallable_Check (gx))
  {
    gravfunc[0] = gx;

    gravtms[0] = -1;
  }
  else if (PyFloat_Check(gx))
  {
    gravity[0] = PyFloat_AsDouble(gx);
  }
  else if (PyInt_Check(gx))
  {
    gravtms[0] = PyInt_AsLong(gx);

    gravfunc[0] = NULL;

    if (gravtms[0] < 0 || gravtms[0] >= tmsnum)
    {
      PyErr_SetString (PyExc_ValueError, "Invalid TSERIES number for gx");
      return NULL;
    }
  }
  else
  {
    PyErr_SetString (PyExc_ValueError, "Invalid gx component");
    return NULL;
  }

  if (PyCallable_Check (gy))
  {
    gravfunc[1] = gy;

    gravtms[1] = -1;
  }
  else if (PyFloat_Check(gy))
  {
    gravity[1] = PyFloat_AsDouble(gy);
  }
  else if (PyInt_Check(gy))
  {
    gravtms[1] = PyInt_AsLong(gy);

    gravfunc[1] = NULL;

    if (gravtms[1] < 0 || gravtms[1] >= tmsnum)
    {
      PyErr_SetString (PyExc_ValueError, "Invalid TSERIES number for gy");
      return NULL;
    }
  }
  else
  {
    PyErr_SetString (PyExc_ValueError, "Invalid gy component");
    return NULL;
  }

  if (PyCallable_Check (gz))
  {
    gravfunc[2] = gz;

    gravtms[2] = -1;
  }
  else if (PyFloat_Check(gz))
  {
    gravity[2] = PyFloat_AsDouble(gz);
  }
  else if (PyInt_Check(gz))
  {
    gravtms[2] = PyInt_AsLong(gz);

    gravfunc[2] = NULL;

    if (gravtms[2] < 0 || gravtms[2] >= tmsnum)
    {
      PyErr_SetString (PyExc_ValueError, "Invalid TSERIES number for gz");
      return NULL;
    }
  }

  else
  {
    PyErr_SetString (PyExc_ValueError, "Invalid gz component");
    return NULL;
  }

  Py_RETURN_NONE;
}

/* set global damping */
static PyObject* DAMPING (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("linear", "angular");
  PyObject *linear, *angular;

  PARSEKEYS ("OO", &linear, &angular);

  if (PyCallable_Check (linear))
  {
    lindamp = linear;

    lindamptms[0] = lindamptms[1] = lindamptms[2] = -1;
  }
  else if (PyTuple_Check(linear))
  {
    if (PyTuple_Size(linear) != 3)
    {
      PyErr_SetString (PyExc_ValueError, "Invalid linear damping tuple size");
      return NULL;
    }

    lindamptms[0] = PyInt_AsLong(PyTuple_GetItem(linear, 0));
    lindamptms[1] = PyInt_AsLong(PyTuple_GetItem(linear, 1));
    lindamptms[2] = PyInt_AsLong(PyTuple_GetItem(linear, 2));

    lindamp = NULL;

    if (lindamptms[0] < 0 || lindamptms[0] >= tmsnum ||
        lindamptms[1] < 0 || lindamptms[1] >= tmsnum ||
        lindamptms[2] < 0 || lindamptms[2] >= tmsnum)
    {
      PyErr_SetString (PyExc_ValueError, "Invalid TSERIES number in linear damping tuple");
      return NULL;
    }
  }
  else
  {
    PyErr_SetString (PyExc_ValueError, "Invalid linear damping argument type");
    return NULL;
  }

  if (PyCallable_Check (angular))
  {
    angdamp = angular;

    angdamptms[0] = angdamptms[1] = angdamptms[2] = -1;
  }
  else if (PyTuple_Check(angular))
  {
    if (PyTuple_Size(angular) != 3)
    {
      PyErr_SetString (PyExc_ValueError, "Invalid angular damping tuple size");
      return NULL;
    }

    angdamptms[0] = PyInt_AsLong(PyTuple_GetItem(angular, 0));
    angdamptms[1] = PyInt_AsLong(PyTuple_GetItem(angular, 1));
    angdamptms[2] = PyInt_AsLong(PyTuple_GetItem(angular, 2));

    angdamp = NULL;

    if (angdamptms[0] < 0 || angdamptms[0] >= tmsnum ||
        angdamptms[1] < 0 || angdamptms[1] >= tmsnum ||
        angdamptms[2] < 0 || angdamptms[2] >= tmsnum)
    {
      PyErr_SetString (PyExc_ValueError, "Invalid TSERIES number in angular damping tuple");
      return NULL;
    }
  }
  else
  {
    PyErr_SetString (PyExc_ValueError, "Invalid angular damping argument type");
    return NULL;
  }

  Py_RETURN_NONE;
}

/* temporary critical step */
struct cristep
{
  REAL step;
  int index;
  REAL omega;
  REAL ratio;
};

/* spring step comparison by step size */
struct cmp_cristep
{
  bool operator() (const cristep& a, const cristep& b)
  {
    if (a.step == b.step) return a.index < b.index;
    else return a.step < b.step;
  }
};

/* estimate critical time step */
static PyObject* CRITICAL (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("perspring", "perparticle");
  int perspring, perparticle;

  perspring = -1;
  perparticle = -1;

  PARSEKEYS ("|ii", &perspring, &perparticle);

  if (perspring < -1 || perspring > sprnum)
  {
    PyErr_SetString (PyExc_ValueError, "Parameter 'perspring' is out of range");
    return NULL;
  }

  if (perparticle < -1 || perparticle > parnum)
  {
    PyErr_SetString (PyExc_ValueError, "Parameter 'perparticle' is out of range");
    return NULL;
  }

  PyObject *perspringout = NULL;

  if (perspring >= 0)
  {
    /* estimate per-spring steps */
    REAL *hcri = new REAL[sprnum];
    REAL *ocri = new REAL[sprnum];
    REAL *rcri = new REAL[sprnum];
    ispc::critical_perspring (sprnum, sprpart, sprpnt, sprdir, spring, spridx,
                 dashpot, dashidx, inverse, invm, position, hcri, ocri, rcri);

    /* sort per-spring steps */
    std::vector<cristep> v;
    v.reserve (sprnum);
    for (int i = 0; i < sprnum; i ++)
    {
      struct cristep x = {hcri[i], i, ocri[i], rcri[i]};
      v.push_back(x);
    }
    std::sort (v.begin(), v.end(), cmp_cristep());

    /* create output list */
    int i = 0;
    perspringout = PyList_New(perspring);
    for (std::vector<cristep>::const_iterator p = v.begin(); p != v.end() && i < perspring; ++p, ++i)
    {
      PyList_SetItem (perspringout, i, Py_BuildValue ("(d, i, d, d)", p->step, p->index, p->omega, p->ratio));
    }
    delete [] hcri;
    delete [] ocri;
    delete [] rcri;
  }

  PyObject *perparticleout = NULL;

  if (perparticle >= 0)
  {
    /* estimate per-particle steps */
    REAL *hcri = new REAL[parnum];
    REAL *ocri = new REAL[parnum];
    REAL *rcri = new REAL[parnum];
    ispc::critical_perparticle (ntasks, master, slave, parnum, mass, inertia, invm, inverse,
                                rotation, position, sprnum, sprflg, sprpart, sprpnt, sprdir,
				spring, spridx, dashpot, dashidx, kact, kmax, emax, krot, hcri,
				ocri, rcri);

    /* sort per-particle steps */
    std::vector<cristep> v;
    v.reserve (parnum);
    for (int i = 0; i < parnum; i ++)
    {
      struct cristep x = {hcri[i], i, ocri[i], rcri[i]};
      v.push_back(x);
    }
    std::sort (v.begin(), v.end(), cmp_cristep());

    /* create output list */
    int i = 0;
    perparticleout = PyList_New(perparticle);
    for (std::vector<cristep>::const_iterator p = v.begin(); p != v.end() && i < perparticle; ++p, ++i)
    {
      PyList_SetItem (perparticleout, i, Py_BuildValue ("(d, i, d, d)", p->step, p->index, p->omega, p->ratio));
    }
    delete [] hcri;
    delete [] ocri;
    delete [] rcri;
  }

  if (perspringout && perparticleout) return Py_BuildValue ("(O, O)", perspringout, perparticleout);
  else if (perspringout) return perspringout;
  else if (perparticleout) return perparticleout;
  else
  {
    REAL h = ispc::critical (parnum, mass, pairnum, iparam, sprnum, sprpart, spring, spridx, dashpot, dashidx);

    return PyFloat_FromDouble (h);
  }
}

/* time history output */
static PyObject* HISTORY (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("entity", "source", "point", "h5file", "h5last");
  PyObject *entity, *source, *point, *h5file, *h5last;
  REAL s[6] = {0., 0., 0., 0., 0., 0.};
  int list_size, *list, kind, hisent;
  int srckind = 0; /* 0: partilce, 1: spring */

  point = NULL;
  source = NULL;
  h5file = NULL;
  h5last = Py_False;

  PARSEKEYS ("O|OOOO", &entity, &source, &point, &h5file, &h5last);

  TYPETEST (is_string (entity, kwl[0]) && is_tuple (point, kwl[2], 3) &&
            is_file (h5file, kwl[3]) && is_bool (h5last, kwl[4]));

  IFIS (entity, "TIME")
  {
    hisent = HIS_TIME;
  }
  ELIF (entity, "PX")
  {
    hisent = HIS_PX;
  }
  ELIF (entity, "PY")
  {
    hisent = HIS_PY;
  }
  ELIF (entity, "PZ")
  {
    hisent = HIS_PZ;
  }
  ELIF (entity, "|P|")
  {
    hisent = HIS_PL;
  }
  ELIF (entity, "DX")
  {
    hisent = HIS_DX;
  }
  ELIF (entity, "DY")
  {
    hisent = HIS_DY;
  }
  ELIF (entity, "DZ")
  {
    hisent = HIS_DZ;
  }
  ELIF (entity, "|D|")
  {
    hisent = HIS_DL;
  }
  ELIF (entity, "VX")
  {
    hisent = HIS_VX;
  }
  ELIF (entity, "VY")
  {
    hisent = HIS_VY;
  }
  ELIF (entity, "VZ")
  {
    hisent = HIS_VZ;
  }
  ELIF (entity, "|V|")
  {
    hisent = HIS_VL;
  }
  ELIF (entity, "OX")
  {
    hisent = HIS_OX;
  }
  ELIF (entity, "OY")
  {
    hisent = HIS_OY;
  }
  ELIF (entity, "OZ")
  {
    hisent = HIS_OZ;
  }
  ELIF (entity, "|O|")
  {
    hisent = HIS_OL;
  }
  ELIF (entity, "FX")
  {
    hisent = HIS_FX;
  }
  ELIF (entity, "FY")
  {
    hisent = HIS_FY;
  }
  ELIF (entity, "FZ")
  {
    hisent = HIS_FZ;
  }
  ELIF (entity, "|F|")
  {
    hisent = HIS_FL;
  }
  ELIF (entity, "TX")
  {
    hisent = HIS_TX;
  }
  ELIF (entity, "TY")
  {
    hisent = HIS_TY;
  }
  ELIF (entity, "TZ")
  {
    hisent = HIS_TZ;
  }
  ELIF (entity, "|T|")
  {
    hisent = HIS_TL;
  }
  ELIF (entity, "LENGTH")
  {
    hisent = HIS_LENGTH;
    srckind = 1;
  }
  ELIF (entity, "STROKE")
  {
    hisent = HIS_STROKE;
    srckind = 1;
  }
  ELIF (entity, "STF")
  {
    hisent = HIS_STF;
    srckind = 1;
  }
  ELIF (entity, "SF")
  {
    hisent = HIS_SF;
    srckind = 1;
  }
  ELIF (entity, "SS")
  {
    hisent = HIS_SS;
    srckind = 1;
  }
  ELSE
  {
    PyErr_SetString (PyExc_ValueError, "Invalid entity");
    return NULL;
  }

  if (source == NULL) /* useful when entity is 'TIME' */
  {
    list_size = 1;
    list = new int;
    list [0] = 0;
    if (!parnum && !h5file)
    {
      PyErr_SetString (PyExc_ValueError, "No particle has been defined");
      return NULL;
    }
    kind = HIS_LIST;
  }
  else if (PyInt_Check (source))
  {
    list_size = 1;
    list = new int;
    list [0] = PyInt_AsLong (source);
    if (!h5file && srckind == 0 && (list[0] < 0 || list[0] >= parnum))
    {
      PyErr_SetString (PyExc_ValueError, "Particle index out of range");
      return NULL;
    }
    if (!h5file && srckind == 1 && (list[0] < 0 || list[0] >= sprnum))
    {
      PyErr_SetString (PyExc_ValueError, "Spring index out of range");
      return NULL;
    }
    kind = HIS_LIST;
  }
  else if (PyList_Check (source))
  {
    list_size = PyList_Size (source);
    list = new int[list_size];
    for (int j = 0; j < list_size; j ++)
    {
      list[j] = PyInt_AsLong (PyList_GetItem(source, j));
      if (!h5file && srckind == 0 && (list[j] < 0 || list[j] >= parnum))
      {
	PyErr_SetString (PyExc_ValueError, "Particle index out of range");
	return NULL;
      }
      if (!h5file && srckind == 1 && (list[j] < 0 || list[j] >= sprnum))
      {
	PyErr_SetString (PyExc_ValueError, "Spring index out of range");
	return NULL;
      }
    }
    kind = HIS_LIST;
  }
  else if (PyTuple_Check (source))
  {
    if (srckind == 1)
    {
      PyErr_SetString (PyExc_ValueError, "Invalid source kind for a spring entity");
      return NULL;
    }

    list_size = 0;
    list = NULL;

    if (PyTuple_Size(source) == 4)
    {
      s[0] = PyFloat_AsDouble(PyTuple_GetItem (source, 0));
      s[1] = PyFloat_AsDouble(PyTuple_GetItem (source, 1));
      s[2] = PyFloat_AsDouble(PyTuple_GetItem (source, 2));
      s[3] = PyFloat_AsDouble(PyTuple_GetItem (source, 3));
      kind = HIS_SPHERE;
    }
    else if (PyTuple_Size(source) == 6)
    {
      s[0] = PyFloat_AsDouble(PyTuple_GetItem (source, 0));
      s[1] = PyFloat_AsDouble(PyTuple_GetItem (source, 1));
      s[2] = PyFloat_AsDouble(PyTuple_GetItem (source, 2));
      s[3] = PyFloat_AsDouble(PyTuple_GetItem (source, 3));
      s[4] = PyFloat_AsDouble(PyTuple_GetItem (source, 4));
      s[5] = PyFloat_AsDouble(PyTuple_GetItem (source, 5));
      kind = HIS_BOX;
    }
    else
    {
      PyErr_SetString (PyExc_ValueError, "Invalid source tuple size");
      return NULL;
    }
  }

  history_buffer_grow (list_size);

  int i = hisnum ++;

  if (list) /* particle list source */
  {
    for (int j = 0; j < list_size; j ++)
    {
      hislst[hisidx[i]+j] = list[j];
    }
    hisidx[i+1] = hisidx[i]+list_size;
    delete list;
  }
  else /* sphere or box source */
  {
    hisidx[i+1] = hisidx[i];
  }

  parmec::hisent[i] = hisent;

  if (kind == HIS_LIST && list_size == 1 && point)
  {
    s[0] = PyFloat_AsDouble (PyTuple_GetItem (point, 0));
    s[1] = PyFloat_AsDouble (PyTuple_GetItem (point, 1));
    s[2] = PyFloat_AsDouble (PyTuple_GetItem (point, 2));

    kind = HIS_LIST|HIS_POINT;
  }

  hiskind[i] = kind;

  parmec::source[0][i] = s[0];
  parmec::source[1][i] = s[1];
  parmec::source[2][i] = s[2];
  parmec::source[3][i] = s[3];
  parmec::source[4][i] = s[4];
  parmec::source[5][i] = s[5];

  history[i] = PyList_New(0);

  if (h5file) parmec::h5file[i] = PyString_AsString(h5file);
  else parmec::h5file[i] = NULL;

  if (h5last == Py_True) output_h5history();

  return (PyObject*) history[i];
}

/* declare output entities */
static PyObject* OUTPUT (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("entities", "subset", "mode", "format");
  PyObject *entities, *subset, *mode, *format;

  subset = NULL;
  mode = NULL;
  format = NULL;
  entities = NULL;

  PARSEKEYS ("|OOOO", &entities, &subset, &mode, &format);

  TYPETEST (is_list (entities, kwl[0], 0) && is_list_or_number (subset, kwl[1], 0) &&
            is_string_or_list (mode, kwl[2]) && is_string_or_list (format, kwl[3]));

  int list_size = 0;

  if (subset)
  {
    if (PyList_Check (subset))
    {
      list_size = PyList_Size (subset);

      for (int j = 0; j < list_size; j ++)
      {
	int k = PyInt_AsLong (PyList_GetItem (subset, j));

	if (k < 0 || k >= parnum)
	{
	  PyErr_SetString (PyExc_ValueError, "Invalid particle number");
	  return NULL;
	}
      }
    }
    else
    {
      int k = PyInt_AsLong (subset);

      if (k < 0 || k >= parnum)
      {
	PyErr_SetString (PyExc_ValueError, "Invalid particle number");
	return NULL;
      }
 
      list_size = 1;
    }
  }

  output_buffer_grow (list_size);

  int i = outnum ++;

  outidx[i+1] = outidx[i];

  if (subset)
  {
    if (PyList_Check (subset))
    {
      for (int j = 0; j < list_size; j ++)
      {
	int k = PyInt_AsLong (PyList_GetItem (subset, j));

	outpart[outidx[i+1]++] = k;

	flags[k] &= ~OUTREST;
      }
    }
    else
    {
      int k = PyInt_AsLong (subset);

      outpart[outidx[i+1]++] = k;

      flags[k] &= ~OUTREST;
    }
  }

  if (mode)
  {
    if (PyString_Check (mode))
    {
      IFIS (mode, "SPH")
      {
	if (subset) outmode[i] = OUT_MODE_SPH;
	else outrest[1] = OUT_MODE_SPH;
      }
      ELIF (mode, "MESH")
      {
	if (subset) outmode[i] = OUT_MODE_MESH;
	else outrest[1] = OUT_MODE_MESH;
      }
      ELIF (mode, "RB")
      {
	if (subset) outmode[i] = OUT_MODE_RB;
	else outmode[1] = OUT_MODE_RB;
      }
      ELIF (mode, "CD")
      {
	if (subset) outmode[i] = OUT_MODE_CD;
	else outmode[1] = OUT_MODE_CD;
      }
      ELIF (mode, "SD")
      {
	if (subset) outmode[i] = OUT_MODE_SD;
	else outmode[i] = OUT_MODE_SD;
      }
      ELSE
      {
	PyErr_SetString (PyExc_ValueError, "Invalid mode");
	return NULL;
      }
    }
    else
    {
      if (subset) outmode[i] = 0;
      else outrest[1] = 0;

      for (int j = 0; j < PyList_Size (mode); j ++)
      {
	PyObject *item = PyList_GetItem (mode, j);

	if (PyString_Check (item))
	{
	  IFIS (item, "SPH")
	  {
	    if (subset) outmode[i] |= OUT_MODE_SPH;
	    else outrest[1] |= OUT_MODE_SPH;
	  }
	  ELIF (item, "MESH")
	  {
	    if (subset) outmode[i] |= OUT_MODE_MESH;
	    else outrest[1] |= OUT_MODE_MESH;
	  }
	  ELIF (item, "RB")
	  {
	    if (subset) outmode[i] |= OUT_MODE_RB;
	    else outmode[1] |= OUT_MODE_RB;
	  }
	  ELIF (item, "CD")
	  {
	    if (subset) outmode[i] |= OUT_MODE_CD;
	    else outmode[1] |= OUT_MODE_CD;
	  }
	  ELIF (item, "SD")
	  {
	    if (subset) outmode[i] |= OUT_MODE_SD;
	    else outmode[i] |= OUT_MODE_SD;
	  }
	  ELSE
	  {
	    PyErr_SetString (PyExc_ValueError, "Invalid mode");
	    return NULL;
	  }
	}
	else
	{
	  PyErr_SetString (PyExc_TypeError, "mode list item is not a string");
	  return NULL;
	}
      }
    }
  }
  else
  {
    if (subset) outmode[i] = OUT_MODE_SPH|OUT_MODE_MESH|OUT_MODE_RB|OUT_MODE_CD|OUT_MODE_SD;
    else outrest[1] = OUT_MODE_SPH|OUT_MODE_MESH|OUT_MODE_RB|OUT_MODE_CD|OUT_MODE_SD; 
  }

  if (entities)
  {
    if (subset) outent[i] = 0;
    else outrest[0] = 0;

    for (int j = 0; j < PyList_Size (entities); j ++)
    {
      PyObject *item = PyList_GetItem (entities, j);

      IFIS (item, "NUMBER")
      {
	if (subset) outent[i] |= OUT_NUMBER;
	else outrest[0] |= OUT_NUMBER;
      }
      ELIF (item, "COLOR")
      {
	if (subset) outent[i] |= OUT_COLOR;
	else outrest[0] |= OUT_COLOR;
      }
      ELIF (item, "DISPL")
      {
	if (subset) outent[i] |= OUT_DISPL;
	else outrest[0] |= OUT_DISPL;
      }
      ELIF (item, "LENGTH")
      {
	if (subset) outent[i] |= OUT_LENGTH;
	else outrest[0] |= OUT_LENGTH;
      }
      ELIF (item, "ORIENT")
      {
	if (subset) outent[i] |= OUT_ORIENT;
	else outrest[0] |= OUT_ORIENT;
      }
      ELIF (item, "ORIENT1")
      {
	if (subset) outent[i] |= OUT_ORIENT1;
	else outrest[0] |= OUT_ORIENT1;
      }
      ELIF (item, "ORIENT2")
      {
	if (subset) outent[i] |= OUT_ORIENT2;
	else outrest[0] |= OUT_ORIENT2;
      }
      ELIF (item, "ORIENT3")
      {
	if (subset) outent[i] |= OUT_ORIENT3;
	else outrest[0] |= OUT_ORIENT3;
      }
      ELIF (item, "LINVEL")
      {
	if (subset) outent[i] |= OUT_LINVEL;
	else outrest[0] |= OUT_LINVEL;
      }
      ELIF (item, "ANGVEL")
      {
	if (subset) outent[i] |= OUT_ANGVEL;
	else outrest[0] |= OUT_ANGVEL;
      }
      ELIF (item, "FORCE")
      {
	if (subset) outent[i] |= OUT_FORCE;
	else outrest[0] |= OUT_FORCE;
      }
      ELIF (item, "TORQUE")
      {
	if (subset) outent[i] |= OUT_TORQUE;
	else outrest[0] |= OUT_TORQUE;
      }
      ELIF (item, "F")
      {
	if (subset) outent[i] |= OUT_F;
	else outrest[0] |= OUT_F;
      }
      ELIF (item, "FN")
      {
	if (subset) outent[i] |= OUT_FN;
	else outrest[0] |= OUT_FN;
      }
      ELIF (item, "FT")
      {
	if (subset) outent[i] |= OUT_FT;
	else outrest[0] |= OUT_FT;
      }
      ELIF (item, "SF")
      {
	if (subset) outent[i] |= OUT_SF;
	else outrest[0] |= OUT_SF;
      }
      ELIF (item, "SS")
      {
	if (subset) outent[i] |= OUT_SS;
	else outrest[0] |= OUT_SS;
      }
      ELIF (item, "AREA")
      {
	if (subset) outent[i] |= OUT_AREA;
	else outrest[0] |= OUT_AREA;
      }
      ELIF (item, "PAIR")
      {
	if (subset) outent[i] |= OUT_PAIR;
	else outrest[0] |= OUT_PAIR;
      }
      ELSE
      {
	PyErr_SetString (PyExc_ValueError, "Invalid entity");
	return NULL;
      }
    }
  }
  else
  {
    if (subset) outent[i] = 0xffffffff;
    else outrest[0] = 0xffffffff;
  }

  if (format)
  {
    if (PyString_Check (format))
    {
      IFIS (format, "VTK")
      {
	parmec::outformat = OUT_FORMAT_VTK;
      }
      ELIF (format, "XDMF")
      {
	parmec::outformat = OUT_FORMAT_XDMF;
      }
      ELIF (format, "MED")
      {
	parmec::outformat = OUT_FORMAT_MED;
      }
      ELSE
      {
	PyErr_SetString (PyExc_ValueError, "Invalid format");
	return NULL;
      }
    }
    else
    {
      parmec::outformat = 0;

      for (int j = 0; j < PyList_Size (format); j ++)
      {
	PyObject *item = PyList_GetItem (format, j);

	IFIS (item, "VTK")
	{
	  parmec::outformat |= OUT_FORMAT_VTK;
	}
	ELIF (item, "XDMF")
	{
	  parmec::outformat |= OUT_FORMAT_XDMF;
	}
	ELIF (item, "MED")
	{
	  parmec::outformat |= OUT_FORMAT_MED;
	}
	ELSE
	{
	  PyErr_SetString (PyExc_ValueError, "Invalid format");
	  return NULL;
	}
      }
    }
  }

  Py_RETURN_NONE;
}

/* run DEM simulation */
static PyObject* DEM (PyObject *self, PyObject *args, PyObject *kwds)
{
  KEYWORDS ("duration", "step", "interval", "prefix", "adaptive");
  double duration, step, adaptive;
  PyObject *prefix, *interval;
  pointer_t dt_func[2];
  int dt_tms[2];
  REAL dt[2];
  char *pre;

  prefix = NULL;
  interval = NULL;
  adaptive = 0.0;

  PARSEKEYS ("dd|OOd", &duration, &step, &interval, &prefix, &adaptive);

  TYPETEST (is_positive (duration, kwl[0]) && is_positive (step, kwl[1]) &&
            is_string (prefix, kwl[3]) && is_ge_le (adaptive, 0.0, 1.0, kwl[4]));

  if (interval)
  {
    if (PyTuple_Check(interval))
    {
      if (PyTuple_Size(interval) != 2)
      {
	PyErr_SetString (PyExc_ValueError, "Invalid output interval");
	return NULL;
      }

      PyObject *dt0 = PyTuple_GetItem (interval, 0);

      if (PyCallable_Check (dt0))
      {
	dt[0] = 0.0;
	dt_func[0] = dt0;
	dt_tms[0] = -1;
      }
      else if (PyInt_Check (dt0))
      {
	dt[0] = 0.0;
	dt_func[0] = NULL;
	dt_tms[0] = PyInt_AsLong (dt0);

	if (dt_tms[0] < 0 || dt_tms[0] >= tmsnum)
	{
	  PyErr_SetString (PyExc_ValueError, "Invalid output interval TSERIES number");
	  return NULL;
	}
      }
      else
      {
	dt[0] = PyFloat_AsDouble (dt0);
	dt_func[0] = NULL;
	dt_tms[0] = -1;
      }

      PyObject *dt1 = PyTuple_GetItem (interval, 1);

      if (PyCallable_Check (dt1))
      {
	dt[1] = 0.0;
	dt_func[1] = dt1;
	dt_tms[1] = -1;
      }
      else if (PyInt_Check (dt1))
      {
	dt[1] = 0.0;
	dt_func[1] = NULL;
	dt_tms[1] = PyInt_AsLong (dt1);

	if (dt_tms[1] < 0 || dt_tms[1] >= tmsnum)
	{
	  PyErr_SetString (PyExc_ValueError, "Invalid output interval TSERIES number");
	  return NULL;
	}
      }
      else
      {
	dt[1] = PyFloat_AsDouble (dt1);
	dt_func[1] = NULL;
	dt_tms[1] = -1;
      }
    }
    else 
    {
      if (PyCallable_Check (interval))
      {
        dt[0] = dt[1] = 0.0;
	dt_func[0] = dt_func[1] = interval;
	dt_tms[0] = dt_tms[1] = -1;
      }
      else if (PyInt_Check (interval))
      {
        dt[0] = dt[1] = 0.0;
	dt_func[0] = dt_func[1] = NULL;
	dt_tms[0] = dt_tms[1] = PyInt_AsLong (interval);

	if (dt_tms[0] < 0 || dt_tms[0] >= tmsnum)
	{
	  PyErr_SetString (PyExc_ValueError, "Invalid output interval TSERIES number");
	  return NULL;
	}
      }
      else
      {
        dt[0] = dt[1] = PyFloat_AsDouble (interval);
	dt_func[0] = dt_func[1] = NULL;
	dt_tms[0] = dt_tms[1] = -1;
      }
    }

    if (dt[0] < 0.0 || dt[1] < 0.0)
    {
      PyErr_SetString (PyExc_ValueError, "Invalid, negative, output interval");
      return NULL;
    }
  }
  else
  {
    dt[0] = dt[1] = step;
    dt_func[0] = dt_func[1] = NULL;
    dt_tms[0] = dt_tms[1] = -1;
  }

  if (prefix)
  {
    pre = PyString_AsString (prefix);
  }
  else pre = NULL;

  duration = dem (duration, step, dt, dt_func, dt_tms, pre, 1, adaptive);

  return Py_BuildValue ("d", duration); /* PyFloat_FromDouble (dt) */
}

static PyMethodDef methods [] =
{
  {"ARGV", (PyCFunction)ARGV, METH_VARARGS|METH_KEYWORDS, "Command line arguments"},
  {"RESET", (PyCFunction)RESET, METH_NOARGS, "Reset simulation"},
  {"TSERIES", (PyCFunction)TSERIES, METH_VARARGS|METH_KEYWORDS, "Create time series"},
  {"MATERIAL", (PyCFunction)MATERIAL, METH_VARARGS|METH_KEYWORDS, "Create material"},
  {"SPHERE", (PyCFunction)SPHERE, METH_VARARGS|METH_KEYWORDS, "Create spherical particle"},
  {"MESH", (PyCFunction)MESH, METH_VARARGS|METH_KEYWORDS, "Create meshed particle"},
  {"ANALYTICAL", (PyCFunction)::ANALYTICAL, METH_VARARGS|METH_KEYWORDS, "Create analytical particle"},
  {"OBSTACLE", (PyCFunction)OBSTACLE, METH_VARARGS|METH_KEYWORDS, "Create obstacle"},
  {"SPRING", (PyCFunction)::SPRING, METH_VARARGS|METH_KEYWORDS, "Create translational spring"},
  {"UNSPRING", (PyCFunction)::UNSPRING, METH_VARARGS|METH_KEYWORDS, "Undo translational springs"},
  {"EQM", (PyCFunction)EQM, METH_VARARGS|METH_KEYWORDS, "Calculate equivalent point mass"},
  {"GRANULAR", (PyCFunction)GRANULAR, METH_VARARGS|METH_KEYWORDS, "Define surface pairing for the granular interaction model"},
  {"RESTRAIN", (PyCFunction)RESTRAIN, METH_VARARGS|METH_KEYWORDS, "Constrain particle motion"},
  {"PRESCRIBE", (PyCFunction)PRESCRIBE, METH_VARARGS|METH_KEYWORDS, "Prescribe particle motion"},
  {"VELOCITY", (PyCFunction)VELOCITY, METH_VARARGS|METH_KEYWORDS, "Set particle velocity"},
  {"GRAVITY", (PyCFunction)GRAVITY, METH_VARARGS|METH_KEYWORDS, "Set gravity"},
  {"DAMPING", (PyCFunction)DAMPING, METH_VARARGS|METH_KEYWORDS, "Set global damping"},
  {"CRITICAL", (PyCFunction)CRITICAL, METH_VARARGS|METH_KEYWORDS, "Estimate critical time step"},
  {"HISTORY", (PyCFunction)HISTORY, METH_VARARGS|METH_KEYWORDS, "Time history output"},
  {"OUTPUT", (PyCFunction)OUTPUT, METH_VARARGS|METH_KEYWORDS, "Declare output entities"},
  {"DEM", (PyCFunction)DEM, METH_VARARGS|METH_KEYWORDS, "Run DEM simulation"},
  {NULL, 0, 0, NULL}
};

namespace parmec
{ /* namespace */

/* interpret an input file (return 0 on success) */
int input (const char *path, char **argv, int argc)
{
  int error, len;
  char *line;

  parmec::argv = argv;
  parmec::argc = argc;

  len = strlen (path);
  ERRMEM (output_path = new char [len+1]);
  strcpy (output_path, path);
 
  if (output_path[len-3] != '.' ||
      output_path[len-2] != 'p' ||
      output_path[len-1] != 'y')
  {
    fprintf (stderr, "ERROR: input file does not have '.py' extension!\n");
    fprintf (stderr, "       the input path reads: %s\n", path);
    return 1;
  }
  else output_path[len-3] = '\0';

  Py_Initialize();

  if (!Py_InitModule3 ("parmec", methods, "parmec module")) return -1;

  PyRun_SimpleString ("from parmec import ARGV\n"
                      "from parmec import RESET\n"
                      "from parmec import TSERIES\n"
                      "from parmec import MATERIAL\n"
                      "from parmec import SPHERE\n"
                      "from parmec import MESH\n"
                      "from parmec import ANALYTICAL\n"
                      "from parmec import OBSTACLE\n"
                      "from parmec import SPRING\n"
                      "from parmec import UNSPRING\n"
                      "from parmec import EQM\n"
                      "from parmec import GRANULAR\n"
                      "from parmec import RESTRAIN\n"
                      "from parmec import PRESCRIBE\n"
                      "from parmec import VELOCITY\n"
                      "from parmec import GRAVITY\n"
                      "from parmec import DAMPING\n"
                      "from parmec import CRITICAL\n"
                      "from parmec import HISTORY\n"
                      "from parmec import OUTPUT\n"
                      "from parmec import DEM\n");

  ERRMEM (line = new char [128 + strlen (path)]);
  sprintf (line, "execfile ('%s')", path);

  error = PyRun_SimpleString (line); /* we do not run a file directly because FILE destriptors differe
					between WIN32 and UNIX while Python is often provided in binary form */
  delete line;

  return error;
}

/* update obstacles time histories from callbacks */
void obstaclev (int obsnum, REAL *obsang, REAL *obslin, pointer_t anghis[], pointer_t linhis[], REAL time)
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

/* prescribe particle acceleration */
void prescribe_acceleration (int prsnum, pointer_t tms[], int prspart[], pointer_t prslin[], int *tmslin[3], int linkind[],
  pointer_t prsang[], int *tmsang[3], int angkind[], REAL time, REAL mass[], REAL *inertia[9], REAL *force[3], REAL *torque[3])
{
  PyObject *result, *args;
  int i, j;

  for (i = 0; i < prsnum; i ++)
  {
    j = prspart[i];

    if ((prslin[i] || tmslin[0][i] >= 0) && linkind[i] == 1) /* prescribed acceleration --> set up force */
    {
      REAL acc[3];

      if (prslin[i])
      {
        args = Py_BuildValue ("(d)", time);

	result = PyObject_CallObject ((PyObject*)prslin[i], args);

	ASSERT (is_tuple (result, "Returned value", 3), "Prescribed linear acceleration callback did not return a (ax, ay, az) tuple");

	acc[0] = PyFloat_AsDouble(PyTuple_GetItem (result, 0));
	acc[1] = PyFloat_AsDouble(PyTuple_GetItem (result, 1));
	acc[2] = PyFloat_AsDouble(PyTuple_GetItem (result, 2));

        Py_DECREF (args);
	Py_DECREF (result);
      }
      else
      {
	acc[0] = TMS_Value ((TMS*)tms[tmslin[0][i]], time);
	acc[1] = TMS_Value ((TMS*)tms[tmslin[1][i]], time);
	acc[2] = TMS_Value ((TMS*)tms[tmslin[2][i]], time);
      }

      REAL ma = mass[j];

      force[0][j] = ma * acc[0];
      force[1][j] = ma * acc[1];
      force[2][j] = ma * acc[2];
    }
    else if ((prslin[i] || tmslin[0][i] >= 0) && linkind[i] == 0) /* prescribed velocity --> zero force */
    {
      force[0][j] = 0.0;
      force[1][j] = 0.0;
      force[2][j] = 0.0;
    }

    if ((prsang[i] || tmsang[0][i] >= 0) && angkind[i] == 1) /* prescribed acceleration --> set up torque */
    {
      REAL acc[3];

      if (prsang[i])
      {
        args = Py_BuildValue ("(d)", time);

	result = PyObject_CallObject ((PyObject*)prsang[i], args);

	ASSERT (is_tuple (result, "Returned value", 3), "Prescribed angular acceleration callback did not return a (ox, oy, oz) tuple");

	acc[0] = PyFloat_AsDouble(PyTuple_GetItem (result, 0));
	acc[1] = PyFloat_AsDouble(PyTuple_GetItem (result, 1));
	acc[2] = PyFloat_AsDouble(PyTuple_GetItem (result, 2));

        Py_DECREF (args);
        Py_DECREF (result);
      }
      else
      {
	acc[0] = TMS_Value ((TMS*)tms[tmsang[0][i]], time);
	acc[1] = TMS_Value ((TMS*)tms[tmsang[1][i]], time);
	acc[2] = TMS_Value ((TMS*)tms[tmsang[2][i]], time);
      }

      REAL in[9] = {inertia[0][j], inertia[1][j], inertia[2][j],
                    inertia[3][j], inertia[4][j], inertia[5][j],
		    inertia[6][j], inertia[7][j], inertia[8][j]};

      REAL to[3];

      NVMUL (in, acc, to);

      torque[0][j] = to[0];
      torque[1][j] = to[1];
      torque[2][j] = to[2];
    }
    else if ((prsang[i] || tmsang[0][i] >= 0) && angkind[i] == 0) /* prescribed velocity --> zero torque */
    {
      torque[0][j] = 0.0;
      torque[1][j] = 0.0;
      torque[2][j] = 0.0;
    }
  }
}

/* prescribe particle velocity */
void prescribe_velocity (int prsnum, pointer_t tms[], int prspart[], pointer_t prslin[], int *tmslin[3], int linkind[],
  pointer_t prsang[], int *tmsang[3], int angkind[], REAL time, REAL *rotation[9], REAL *linear[3], REAL *angular[6])
{
  PyObject *result, *args;
  int i, j;

  for (i = 0; i < prsnum; i ++)
  {
    j = prspart[i];

    if ((prslin[i] || tmslin[0][i] >= 0) && linkind[i] == 0)
    {
      if (prslin[i])
      {
	args = Py_BuildValue ("(d)", time);

	result = PyObject_CallObject ((PyObject*)prslin[i], args);

	ASSERT (is_tuple (result, "Returned value", 3), "Prescribed linear velocity callback did not return a (vx, vy, vz) tuple");
	
	linear[0][j] = PyFloat_AsDouble(PyTuple_GetItem (result, 0));
	linear[1][j] = PyFloat_AsDouble(PyTuple_GetItem (result, 1));
	linear[2][j] = PyFloat_AsDouble(PyTuple_GetItem (result, 2));

	Py_DECREF (args);
	Py_DECREF (result);
      }
      else
      {
	linear[0][j] = TMS_Value((TMS*)tms[tmslin[0][i]], time);
	linear[1][j] = TMS_Value((TMS*)tms[tmslin[1][i]], time);
	linear[2][j] = TMS_Value((TMS*)tms[tmslin[2][i]], time);
      }
    }

    if ((prsang[i] || tmsang[0][i] >= 0) && angkind[i] == 0)
    {
      REAL o[3];

      if (prsang[i])
      {
	args = Py_BuildValue ("(d)", time);

	result = PyObject_CallObject ((PyObject*)prsang[i], args);

	ASSERT (is_tuple (result, "Returned value", 3), "Prescribed angular velocity callback did not return a (ox, oy, oz) tuple");

        o[0] = PyFloat_AsDouble(PyTuple_GetItem (result, 0));
	o[1] = PyFloat_AsDouble(PyTuple_GetItem (result, 1));
	o[2] = PyFloat_AsDouble(PyTuple_GetItem (result, 2));

	Py_DECREF (args);
        Py_DECREF (result);
      }
      else
      {
	o[0] = TMS_Value((TMS*)tms[tmsang[0][i]], time);
	o[1] = TMS_Value((TMS*)tms[tmsang[1][i]], time);
	o[2] = TMS_Value((TMS*)tms[tmsang[2][i]], time);
      }

      REAL L[9] = {rotation[0][j], rotation[1][j], rotation[2][j],
                   rotation[3][j], rotation[4][j], rotation[5][j],
		   rotation[6][j], rotation[7][j], rotation[8][j]};

      REAL O[3];

      TVMUL (L, o, O);

      angular[0][j] = O[0];
      angular[1][j] = O[1];
      angular[2][j] = O[2];
      angular[3][j] = o[0];
      angular[4][j] = o[1];
      angular[5][j] = o[2];
    }
  }
}

/* read gravity and global damping */
void read_gravity_and_damping (REAL time, pointer_t *tms, pointer_t gravfunc[3], int gravtms[3],
  REAL gravity[3], pointer_t lindamp, int lindamptms[3], pointer_t angdamp, int angdamptms[3], REAL damping[6])
{
  for (int i = 0; i < 3; i ++)
  {
    if (gravfunc[i])
    {
      PyObject *result, *args;

      args = Py_BuildValue ("(d)", time);

      result = PyObject_CallObject ((PyObject*)gravfunc[i], args);
      ASSERT (PyNumber_Check (result), "Gravity callback component %d did not return a number", i);
      gravity[i] = PyFloat_AsDouble(result);
      Py_DECREF (result);

      Py_DECREF (args);
    }
    else if (gravtms[i] >= 0 && gravtms[i] < tmsnum)
    {
      gravity[i] = TMS_Value ((TMS*)tms[gravtms[i]], time);
    }
  }

  if (lindamp)
  {
    PyObject *result, *args;

    args = Py_BuildValue ("(d)", time);

    result = PyObject_CallObject ((PyObject*)lindamp, args);
    ASSERT (is_tuple (result, "Returned value", 3), "Prescribed linear damping did not return a (dvx, dvy, dvz) tuple");
    damping[0] = PyFloat_AsDouble(PyTuple_GetItem (result, 0));
    damping[1] = PyFloat_AsDouble(PyTuple_GetItem (result, 1));
    damping[2] = PyFloat_AsDouble(PyTuple_GetItem (result, 2));
    Py_DECREF (result);

    Py_DECREF (args);
  }
  else if (lindamptms[0] >= 0 && lindamptms[0] < tmsnum &&
           lindamptms[1] >= 0 && lindamptms[1] < tmsnum &&
           lindamptms[2] >= 0 && lindamptms[2] < tmsnum)
  {
    damping[0] = TMS_Value ((TMS*)tms[lindamptms[0]], time);
    damping[1] = TMS_Value ((TMS*)tms[lindamptms[1]], time);
    damping[2] = TMS_Value ((TMS*)tms[lindamptms[2]], time);
  }
  else
  {
    damping[0] = damping[1] = damping[2] = 0.0;
  }

  if (angdamp)
  {
    PyObject *result, *args;

    args = Py_BuildValue ("(d)", time);

    result = PyObject_CallObject ((PyObject*)angdamp, args);
    ASSERT (is_tuple (result, "Returned value", 3), "Prescribed linear damping did not return a (dox, doy, doz) tuple");
    damping[3] = PyFloat_AsDouble(PyTuple_GetItem (result, 0));
    damping[4] = PyFloat_AsDouble(PyTuple_GetItem (result, 1));
    damping[5] = PyFloat_AsDouble(PyTuple_GetItem (result, 2));
    Py_DECREF (result);

    Py_DECREF (args);
  }
  else if (angdamptms[0] >= 0 && angdamptms[0] < tmsnum &&
           angdamptms[1] >= 0 && angdamptms[1] < tmsnum &&
           angdamptms[2] >= 0 && angdamptms[2] < tmsnum)
  {
    damping[3] = TMS_Value ((TMS*)tms[angdamptms[0]], time);
    damping[4] = TMS_Value ((TMS*)tms[angdamptms[1]], time);
    damping[5] = TMS_Value ((TMS*)tms[angdamptms[2]], time);
  }
  else
  {
    damping[3] = damping[4] = damping[5] = 0.0;
  }
}

/* call interval callback */
REAL current_interval (pointer_t func, REAL time)
{
  PyObject *result, *args;
  REAL dt;

  args = Py_BuildValue ("(d)", time);

  result = PyObject_CallObject ((PyObject*)func, args);

  dt = PyFloat_AsDouble(result);

  Py_DECREF (result);

  return dt;
}
} /* namespace */

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
#include "parmec_ispc.h"
#include "partition_ispc.h"
#include "forces_ispc.h"
#include "dynamics_ispc.h"
#include "shapes_ispc.h"
#include "obstacles_ispc.h"

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

/* create sphere */
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

  ellrng[0][i] = j;
  ellrng[1][i] = j+1;

  double volume = (4./3.)*M_PI*rad*rad*rad;

  mass[i] = volume*mparam[DENSITY][material];

  rotation[0][i] = rotation[4][i] = rotation[8][i] = 1.0;
  rotation[1][i] = rotation[2][i] = rotation[3][i] =
  rotation[5][i] = rotation[6][i] = rotation[7][i] = 0.0;

  inertia[0][i] = inertia[4][i] = inertia[8][i] = 0.4*mass[i]*radii[0][j]*radii[0][j];
  inertia[1][i] = inertia[2][i] = inertia[3][i] =
  inertia[5][i] = inertia[6][i] = inertia[7][i] = 0.0;

  return PyInt_FromLong (i);
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
  iparam[SPRING][i] = spring;
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

    shapes (threads, ellnum, part, center, radii, orient, rotation, position);

    obstaclev (obsnum, obsang, obslin, anghis, linhis, 0.5*step);

    obstacles (threads, obsnum, trirng, obspnt, obsang, obslin, trinum, tri, step);
  }

  invert_inertia (threads, parnum, inertia, inverse, mass, invm);

  partitioning *tree = partitioning_create (threads, ellnum, center);

  /* time stepping */
  for (t0 = time = 0.0; time < duration; time += step)
  {
    if (partitioning_store (threads, tree, ellnum, ellcol, part, center, radii, orient) > 0)
    {
      partitioning_destroy (tree);

      tree = partitioning_create (threads, ellnum, center);

      ASSERT (partitioning_store (threads, tree, ellnum, ellcol, part, center, radii, orient) == 0, "Repartitioning failed");
    }

    condet (threads, tree, master, parnum, ellnum, ellcol, part,
            center, radii, orient, trinum, tricol, triobs, tri);

    forces (threads, master, slave, parnum, angular, linear, rotation, position, inertia, inverse,
            mass, invm, obspnt, obslin, obsang, parmat, mparam, pairnum, pairs, ikind, iparam, step);

    dynamics (threads, master, slave, parnum, angular, linear, rotation,
      position, inertia, inverse, mass, invm, force, torque, gravity, step);

    shapes (threads, ellnum, part, center, radii, orient, rotation, position);

    obstaclev (obsnum, obsang, obslin, anghis, linhis, time+step);

    obstacles (threads, obsnum, trirng, obspnt, obsang, obslin, trinum, tri, step);

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
  {"OBSTACLE", (PyCFunction)OBSTACLE, METH_VARARGS|METH_KEYWORDS, "Create obstacle"},
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

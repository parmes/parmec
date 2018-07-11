/*
The MIT License (MIT)

Copyright (c) 2018 EDF Energy

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

/* Contributors: Tomasz Koziara */

#include <map>
#include <set>

#include "macros.h"
#include "joints.h"
#include "solver.hpp"

namespace parmec {

typedef amgcl::backend::builtin<REAL> Backend;

#define DIRECT_SOLVE 1

#if DIRECT_SOLVE
typedef amgcl::solver::skyline_lu<REAL> Solver;
#else
typedef amgcl::solver::bicgstab<Backend> Solver;

typedef amgcl::amg<Backend,
  amgcl::coarsening::smoothed_aggregation,
  amgcl::relaxation::spai0> Precond;

typedef std::tuple<ptrdiff_t,
  std::vector<ptrdiff_t>,
  std::vector<ptrdiff_t>,
  std::vector<REAL>> System;

Precond *precond = NULL;

System *system = NULL;
#endif

Solver *solve = NULL;

/* reset joints matrix after change */
void reset_joints_matrix (int jnum, int *jpart[2], REAL *jpoint[3],
  REAL *position[6], REAL *rotation[9], REAL *inverse[9], REAL invm[],
  bool update_precond)
{
  /* matrix in CRS format */
  std::vector<ptrdiff_t> ptr;
  std::vector<ptrdiff_t> col;

  /* adjacency map of genralized inverse inertia */
  std::map<int,std::set<int>> partadj; /* particle to joint mapping */
  std::vector<std::map<int,int>> matadj(3*jnum); /* matrix adjacency structure and value mapping */

  for (int i = 0; i < jnum; i ++)
  {
    int part1 = jpart[0][i], part2 = jpart[1][i];

    partadj[part1].insert(i);
    if (part2 >= 0) partadj[part2].insert(i);
  }

  for (std::map<int,std::set<int>>::iterator it = partadj.begin(); it != partadj.end(); it ++)
  {
    std::set<int> &adj = it->second;
    for (std::set<int>::iterator jt = adj.begin(); jt != adj.end(); jt ++)
    {
      for (std::set<int>::iterator kt = adj.begin(); kt != adj.end(); kt ++)
      {
	for (int i0 = 0; i0 < 3; i0 ++)
	{
	  for (int j0 = 0; j0 < 3; j0 ++)
	  {
	    matadj[(*jt)*3+i0][(*kt)*3+j0] = 0; /* first create adjacency with zero index value mapping */
	  }
	}
      }
    }
  }

  /* next create 'ptr' and 'col' structures and set up value mapping */
  int j = 0;
  ptr.push_back(0);
  for (int i = 0; i < 3*jnum; i ++)
  {
    std::map<int,int> &adj = matadj[i];

    for (std::map<int,int>::iterator it = adj.begin(); it != adj.end(); it ++)
    {
      col.push_back(it->first);
      it->second = j ++;
    }

    ptr.push_back(ptr.back()+adj.size());
  }
  std::vector<REAL> val(j, 0.); /* create zeroed matrix values */

#if 0
  for (int i = 0; i < 3*jnum; i ++)
  {
    std::map<int,int> &adj = matadj[i];

    for (std::map<int,int>::iterator it = adj.begin(); it != adj.end(); it ++)
    {
      printf ("it->first = %d, it->second = %d\n", it->first, it->second);
    }
  }
#endif
	
  /* TODO: use OMP pragma to parallelise this loop */

  /* assemble joints matrix */
  for (int i = 0; i < jnum; i ++)
  {
    REAL Jiv[9], Rot[9], im, A[3], Ask[9], Hi[9], Hj[9], C[9], D[9], Wii[9], Wij[9];

    SET9(Wii, 0.); 

    for (int k = 0; k < 2; k ++)
    {
      int part = jpart[k][i];

      if (part >= 0)
      {
	Rot[0] = rotation[0][part];
	Rot[1] = rotation[1][part];
	Rot[2] = rotation[2][part];
	Rot[3] = rotation[3][part];
	Rot[4] = rotation[4][part];
	Rot[5] = rotation[5][part];
	Rot[6] = rotation[6][part];
	Rot[7] = rotation[7][part];
	Rot[8] = rotation[8][part];
	Jiv[0] = inverse[0][part];
	Jiv[1] = inverse[1][part];
	Jiv[2] = inverse[2][part];
	Jiv[3] = inverse[3][part];
	Jiv[4] = inverse[4][part];
	Jiv[5] = inverse[5][part];
	Jiv[6] = inverse[6][part];
	Jiv[7] = inverse[7][part];
	Jiv[8] = inverse[8][part];
	im = invm[part];
	A[0] = position[3][part] - jpoint[0][i];
	A[1] = position[4][part] - jpoint[1][i];
	A[2] = position[5][part] - jpoint[2][i];
	VECSKEW (A, Ask);
	NNMUL (Rot, Ask, Hi);
	NTMUL (Jiv, Hi, C);
	NNMUL (Hi, C, D);
	NNADD (Wii, D, Wii);
	Wii[0] += im;
	Wii[4] += im;
	Wii[8] += im;

	std::set<int> &adj = partadj[part];
	for (std::set<int>::iterator it = adj.begin(); it != adj.end(); it ++)
	{
	  int j = *it;
	  if (i != j) /* Wij */
	  {
	    A[0] = position[3][part] - jpoint[0][j];
	    A[1] = position[4][part] - jpoint[1][j];
	    A[2] = position[5][part] - jpoint[2][j];
	    VECSKEW (A, Ask);
	    NNMUL (Rot, Ask, Hj);
	    NTMUL (Jiv, Hj, C);
	    NNMUL (Hi, C, Wij);
	    Wij[0] += im;
	    Wij[4] += im;
	    Wij[8] += im;

	    if ((part == jpart[0][i] && part == jpart[1][j]) ||
	        (part == jpart[1][i] && part == jpart[0][j]))
	    {
	      SCALE9(Wij, -1.0);
	    }

	    /* accumulate Wij into (ptr, col, val) */
	    for (int i0 = 0; i0 < 3; i0 ++)
	    {
	      for (int j0 = 0; j0 < 3; j0 ++)
	      {
		int k0 = matadj[i*3+i0][j*3+j0];
		val[k0] += Wij[i0*3+j0];
	      }
	    }
	  }
	}
      }
    }

    /* accumulate Wii into (ptr, col, val) */
    for (int i0 = 0; i0 < 3; i0 ++)
    {
      for (int j0 = 0; j0 < 3; j0 ++)
      {
        int k0 = matadj[i*3+i0][i*3+j0];
        val[k0] += Wii[i0*3+j0];
      }
    }
  }

  /* problem size */
  int n = ptr.size() - 1;

  /* set up solver */
#if DIRECT_SOLVE
  auto A = amgcl::adapter::zero_copy(n, ptr.data(), col.data(), val.data());
  solve = new Solver(*A);
#else
  Precond::params pprm;
  Solver::params sprm;
  pprm.coarsening.aggr.block_size = 3;
  sprm.tol = REAL_EPS;
  sprm.maxiter = 100;

  if (system) delete system;

  system = new System(n, ptr, col, val);

  if (precond == NULL || update_precond)
  {
    if (precond) delete precond;

    precond = new Precond(*system, pprm);
  }

  if (solve) delete solve;

  solve = new Solver(n, sprm);
#endif
}

/* solve joints and update forces */
void solve_joints (int jnum, int *jpart[2], REAL *jpoint[3], REAL *jreac[3],
  REAL *position[6], REAL *rotation[9], REAL *inertia[9], REAL *inverse[9], REAL invm[],
  REAL *linear[3], REAL *angular[6], REAL *force[6], REAL *torque[6], REAL step)
{
  REAL half = 0.5*step;

  /* assemble joints-free local velocity vector */
  std::vector<REAL> rhs(3*jnum, 0.);
  REAL *B = &rhs[0];

  for (int i = 0; i < jnum; i ++, B += 3)
  {
    REAL Rot[9], Jiv[9], J[9], O[3], t[3], A[3], C[9], T[3], dRot[9], Hi[9], im;

    for (int k = 0; k < 2; k ++)
    {
      int part = jpart[k][i];

      if (part >= 0)
      {
	Rot[0] = rotation[0][part];
	Rot[1] = rotation[1][part];
	Rot[2] = rotation[2][part];
	Rot[3] = rotation[3][part];
	Rot[4] = rotation[4][part];
	Rot[5] = rotation[5][part];
	Rot[6] = rotation[6][part];
	Rot[7] = rotation[7][part];
	Rot[8] = rotation[8][part];
	Jiv[0] = inverse[0][part];
	Jiv[1] = inverse[1][part];
	Jiv[2] = inverse[2][part];
	Jiv[3] = inverse[3][part];
	Jiv[4] = inverse[4][part];
	Jiv[5] = inverse[5][part];
	Jiv[6] = inverse[6][part];
	Jiv[7] = inverse[7][part];
	Jiv[8] = inverse[8][part];
	J[0] = inertia[0][part];
	J[1] = inertia[1][part];
	J[2] = inertia[2][part];
	J[3] = inertia[3][part];
	J[4] = inertia[4][part];
	J[5] = inertia[5][part];
	J[6] = inertia[6][part];
	J[7] = inertia[7][part];
	J[8] = inertia[8][part];
	O[0] = angular[0][part];
	O[1] = angular[1][part];
	O[2] = angular[2][part];
	t[0] = torque[0][part];
	t[1] = torque[1][part];
	t[2] = torque[2][part];
	im = invm[part];
      
	A[0] = position[3][part] - jpoint[0][i];
	A[1] = position[4][part] - jpoint[1][i];
	A[2] = position[5][part] - jpoint[2][i];
	VECSKEW (A, C);
	NNMUL (Rot, C, Hi);

	expmap (-half*O[0], -half*O[1], -half*O[2], dRot[0], dRot[1],
	  dRot[2], dRot[3], dRot[4], dRot[5], dRot[6], dRot[7], dRot[8]);

	TVMUL (Rot, t, T);
	NVMUL (J, O, A);
	NVMUL (dRot, A, C);
	ADDMUL (C, half, T, C);
	NVMUL (Jiv, C, A); /* O(t+h/2) */

	NVMUL (J, A, C);
	PRODUCTSUB (A, C, T); /* T - O(t+h/2) x J O(t+h/2) */

	SCALE (T, step);
	NVADDMUL (O, Jiv, T, A); /* O(t+h) */
	NVMUL (Hi, A, C);

	if (part == jpart[0][i])
	{
	  B[0] -= C[0];
	  B[1] -= C[1];
	  B[2] -= C[2];
	}
	else
	{
	  B[0] += C[0];
	  B[1] += C[1];
	  B[2] += C[2];
	}

        C[0] = force[0][part];
        C[1] = force[1][part];
        C[2] = force[2][part];
	SCALE (C, step*im);
	C[0] += linear[0][part];
	C[1] += linear[1][part];
	C[2] += linear[2][part];

	if (part == jpart[0][i])
	{
	  B[0] -= C[0];
	  B[1] -= C[1];
	  B[2] -= C[2];
	}
	else
	{
	  B[0] += C[0];
	  B[1] += C[1];
	  B[2] += C[2];
	}
      }
    }
  }

  /* solve for joint forces */
  std::vector<REAL> x(rhs.size());
#if DIRECT_SOLVE
  (*solve)(rhs, x);
#else
  int iters;
  REAL error;
  REAL *R0 = &x[0];
  for (int i = 0; i < jnum; i ++, R0 += 3)
  {
    R0[0] = step*jreac[0][i];
    R0[1] = step*jreac[1][i];
    R0[2] = step*jreac[2][i];
  }
  std::tie(iters, error) = (*solve)(*system, *precond, rhs, x);
#endif

  /* accumulate joint forces into body force and torque vectors */
  REAL *R = &x[0];
  REAL invstep = 1./step;
  for (int i = 0; i < jnum; i ++, R += 3)
  {
    REAL Rot[9], A[3], a[3], t[3];

    SCALE (R, invstep); /* 0 = B + W hR --> x = hR --> R = x/h */

    jreac[0][i] = R[0]; /* save solution */
    jreac[1][i] = R[1];
    jreac[2][i] = R[2];

    for (int k = 0; k < 2; k ++)
    {
      int part = jpart[k][i];

      if (part >= 0)
      {
	Rot[0] = rotation[0][part];
	Rot[1] = rotation[1][part];
	Rot[2] = rotation[2][part];
	Rot[3] = rotation[3][part];
	Rot[4] = rotation[4][part];
	Rot[5] = rotation[5][part];
	Rot[6] = rotation[6][part];
	Rot[7] = rotation[7][part];
	Rot[8] = rotation[8][part];
	A[0] = position[3][part] - jpoint[0][i];
	A[1] = position[4][part] - jpoint[1][i];
	A[2] = position[5][part] - jpoint[2][i];
	NVMUL (Rot, A, a);
	PRODUCT (R, a, t);

	if (part == jpart[0][i])
	{
	  torque[0][part] += t[0];
	  torque[1][part] += t[1];
	  torque[2][part] += t[2];
	  force[0][part] += R[0];
	  force[1][part] += R[1];
	  force[2][part] += R[2];
	}
	else
	{
	  torque[0][part] -= t[0];
	  torque[1][part] -= t[1];
	  torque[2][part] -= t[2];
	  force[0][part] -= R[0];
	  force[1][part] -= R[1];
	  force[2][part] -= R[2];
	}
      }
    }
  }
}

} /* namespace */

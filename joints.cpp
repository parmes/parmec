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

#include <vector>
#include <tuple>
#include <map>
#include <set>

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "macros.h"
#include "joints.h"
#include "mem.h"

#if SUITESPARSE
#include "SuiteSparseQR.hpp"
#else
#include "amgcl.hpp"
#endif

namespace parmec {

#if SUITESPARSE
typedef SuiteSparseQR_factorization <REAL> Solver; /* QR factorization */
typedef SuiteSparse_long Long;
#else 
typedef amgcl::solver::skyline_lu<REAL> Solver; /* LU factorization */
typedef ptrdiff_t Long;
#endif

typedef std::tuple<Long,
  std::vector<Long>,
  std::vector<Long>,
  std::vector<REAL>> System;  /* System matrix */

std::map<int,std::set<int>> partadj; /* particle to joint mapping */
std::vector<int> jmap; /* joint index to independent joint set index mapping */
std::vector<std::vector<int>> isets; /* independent joint sets */
std::vector<std::vector<std::map<int,int>>> matadjs; /* matrix adjacency structures and value mappings */
std::vector<System*> systems; /* linear systems */

static void jmark (std::vector<std::set<int>> &jadj, std::set<int> &jall, int joint, std::vector<int> &iset)
{
  if (jall.find(joint) != jall.end())
  {
    jall.erase (joint);

    iset.push_back (joint);

    for (std::set<int>::iterator it = jadj[joint].begin(); it != jadj[joint].end(); it++)
    {
      jmark (jadj, jall, *it, iset);  
    }
  }
}

/* reset joints symbolic information after change */
void reset_joints_symbolic (int jnum, int *jpart[2])
{
  /* clear symbolic data */
  for (int i = 0; i < systems.size(); i ++) delete systems[i];
  systems.clear();
  partadj.clear();
  matadjs.clear();
  isets.clear();
  jmap.clear();

  /* create joint to particle adjacency mapping */
  for (int i = 0; i < jnum; i ++)
  {
    int part1 = jpart[0][i], part2 = jpart[1][i];

    partadj[part1].insert(i);
    if (part2 >= 0) partadj[part2].insert(i);
  }

  /* detect independent joint sets */
  std::vector<std::set<int>> jadj(jnum); /* joint adjacency */
  for (std::map<int,std::set<int>>::iterator it = partadj.begin(); it != partadj.end(); it ++)
  {
    std::set<int> &adj = it->second;
    for (std::set<int>::iterator jt = adj.begin(); jt != adj.end(); jt ++)
    {
      for (std::set<int>::iterator kt = adj.begin(); kt != adj.end(); kt ++)
      {
	if (*jt != *kt)
	{
          jadj[*jt].insert (*kt);
	}
      }
    }
  }

  /* create isets, jmap */
  jmap.resize(jnum);
  std::set<int> jall;
  for (int i = 0; i < jnum; i ++) jall.insert(i);
  while (jall.size() > 0)
  {
    std::vector<int> iset;

    int joint = *jall.begin();

    jmark (jadj, jall, joint, iset);

    for (int i = 0; i < iset.size(); i ++) jmap[iset[i]] = i;

    isets.push_back(iset);
  }

  for (int is = 0; is < isets.size(); is ++)
  {
    std::vector<Long> ptr; /* matrix in CRS format row pointers */
    std::vector<Long> col; /* and column indices */

    std::vector<int> &iset = isets[is];

    std::vector<std::map<int,int>> matadj(3*iset.size()); /* matrix adjacency structure and value mapping */

    for (int i0 = 0; i0 < iset.size(); i0 ++)
    {
      for (int j0 = 0; j0 < iset.size(); j0 ++)
      {
	for (int i1 = 0; i1 < 3; i1 ++)
	{
	  for (int j1 = 0; j1 < 3; j1 ++)
	  {
	    matadj[i0*3+i1][j0*3+j1] = 0; /* first create adjacency with zero index value mapping */
	  }
	}
      }
    }

    /* next create 'ptr' and 'col' structures and set up value mapping */
    int j = 0;
    ptr.push_back(0);
    for (int i = 0; i < 3*iset.size(); i ++)
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

    int n = ptr.size() - 1; /* problem size */

    matadjs.push_back(matadj);

    systems.push_back(new System (n, ptr, col, val));
  }
}

/* solve joints and update forces */
void solve_joints (int jnum, int *jpart[2], REAL *jpoint[3], REAL *jreac[3], int parnum,
  REAL *position[6], REAL *rotation[9], REAL *inertia[9], REAL *inverse[9], REAL mass[], REAL invm[],
  REAL damping[6], REAL *linear[3], REAL *angular[6], REAL *force[6], REAL *torque[6], REAL step0, REAL step1)
{
  REAL half = 0.5*step0;
  REAL step = 0.5*(step0+step1);

#if defined(_OPENMP)
  std::vector<omp_lock_t> lock(parnum);
  for (int i = 0; i < parnum; i ++) omp_init_lock (&lock[i]);

  #pragma omp parallel for
#endif
  for (int is = 0; is < isets.size(); is ++)
  {
    int n = std::get<0>(*systems[is]); /* system size */
    std::vector<ptrdiff_t> &ptr = std::get<1>(*systems[is]); /* row pointers */
    std::vector<ptrdiff_t> &col = std::get<2>(*systems[is]); /* column indices */
    std::vector<REAL> &val = std::get<3>(*systems[is]); /* matrix values */
    std::vector<std::map<int,int>> &matadj = matadjs[is]; /* matrix adjacency and value mapping */
    std::vector<int> &iset = isets[is]; /* independent joint set */

    /* zero matrix values */
    for (int i = 0; i < ptr[n]; i ++) val[i] = 0.0;
	  
    /* assemble joints matrix */
    for (std::vector<int>::iterator it = iset.begin(); it != iset.end(); it++)
    {
      REAL Jiv[9], Rot[9], im, A[3], Ask[9], Hi[9], Hj[9], C[9], D[9], Wii[9], Wij[9];

      int i0 = it-iset.begin();

      int i = *it;

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
	    int j0 = jmap[*it];
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
	      for (int i1 = 0; i1 < 3; i1 ++)
	      {
		for (int j1 = 0; j1 < 3; j1 ++)
		{
		  int k0 = matadj[i0*3+i1][j0*3+j1];
		  val[k0] += Wij[i1*3+j1];
		}
	      }
	    }
	  }
	}
      }

      /* accumulate Wii into (ptr, col, val) */
      for (int i1 = 0; i1 < 3; i1 ++)
      {
	for (int j1 = 0; j1 < 3; j1 ++)
	{
	  int k0 = matadj[i0*3+i1][i0*3+j1];
	  val[k0] += Wii[i1*3+j1];
	}
      }
    }

    /* create linear solver */
#if SUITESPARSE
    cholmod_common cc;
    cholmod_sparse A;

    A.nrow = n;
    A.ncol = n;
    A.p = ptr.data();
    A.i = col.data();
    A.nz = NULL;
    A.x = val.data();
    A.z = NULL;
    A.stype = 0;
    A.itype = CHOLMOD_LONG;
    A.xtype = CHOLMOD_REAL;
#if REAL_SIZE==4
    A.dtype = CHOLMOD_SINGLE;
#else
    A.dtype = CHOLMOD_DOUBLE;
#endif
    A.sorted = 1;
    A.packed = 1;

    cholmod_l_start (&cc);

    Solver *solver = SuiteSparseQR_factorize <REAL> (SPQR_ORDERING_DEFAULT, SPQR_DEFAULT_TOL, &A, &cc) ;
#else
    auto A = amgcl::adapter::zero_copy(n, ptr.data(), col.data(), val.data());

    Solver *solver = new Solver(*A);
#endif

    /* assemble joints-free local velocity vector */
    std::vector<REAL> rhs(3*iset.size(), 0.);
    REAL *B = &rhs[0];

    for (std::vector<int>::iterator it = iset.begin(); it != iset.end(); it++, B += 3)
    {
      REAL Rot[9], Jiv[9], J[9], O[3], o[3], v[3], f[3], t[3], A[3], V[3], C[9], T[3], dRot[9], Hi[9], im;

      int i = *it;

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
	  f[0] = force[0][part];
	  f[1] = force[1][part];
	  f[2] = force[2][part];
	  t[0] = torque[0][part];
	  t[1] = torque[1][part];
	  t[2] = torque[2][part];
	  v[0] = linear[0][part];
	  v[1] = linear[1][part];
	  v[2] = linear[2][part];
	  im = invm[part];
	
	  A[0] = position[3][part] - jpoint[0][i];
	  A[1] = position[4][part] - jpoint[1][i];
	  A[2] = position[5][part] - jpoint[2][i];
	  VECSKEW (A, C);
	  NNMUL (Rot, C, Hi);

	  NVMUL (Hi, O, V);
	  ACC (v, V); /* U(t) */

	  if (damping[0] != 0.0 || damping[1] != 0.0 || damping[2] != 0.0)
	  {
	    REAL ma = mass[i];

	    f[0] -= ma * damping[0] * v[0];
	    f[1] -= ma * damping[1] * v[1];
	    f[2] -= ma * damping[2] * v[2];
	  }
	  if (damping[3] != 0.0 || damping[4] != 0.0 || damping[5] != 0.0)
	  {
	    o[0] = damping[3]*angular[3][part];
	    o[1] = damping[4]*angular[4][part];
	    o[2] = damping[5]*angular[5][part];

	    TVMUL (Rot, o, C);
	    NVMUL (J, C, T);
	    NVMUL (Rot, T, C);

	    t[0] -= C[0];
	    t[1] -= C[1];
	    t[2] -= C[2];
	  }

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

	  NVMUL (Hi, A, C); /* rel. joint vel. U(t+h) = H O(t+h) + ... */
          ACC (v, C);
	  ADDMUL (C, step*im, f, C) /* ... U(t+h) += v(t+h) {== v(t) + step*force/mass} */

	  if (part == jpart[0][i])
	  {
	    B[0] -= V[0]+C[0]; /* U(t) + U(t+h) = 0 = V + C + W*step*R */
	    B[1] -= V[1]+C[1];
	    B[2] -= V[2]+C[2];
	  }
	  else
	  {
	    B[0] += V[0]+C[0];
	    B[1] += V[1]+C[1];
	    B[2] += V[2]+C[2];
	  }
	}
      }
    }

    /* solve for joint forces */
    std::vector<REAL> x(rhs.size());
#if SUITESPARSE
    cholmod_dense *X, *Y, Z;

    cholmod_l_start (&cc);

    Z.nrow = rhs.size();
    Z.ncol = 1;
    Z.nzmax = Z.nrow;
    Z.d = Z.nrow;
    Z.x = rhs.data();
    Z.z = NULL;
    Z.xtype = CHOLMOD_REAL;
#if REAL_SIZE==4
    Z.dtype = CHOLMOD_SINGLE;
#else
    Z.dtype = CHOLMOD_DOUBLE;
#endif

    Y = SuiteSparseQR_qmult (SPQR_QTX, solver, &Z, &cc);
    X = SuiteSparseQR_solve (SPQR_RETX_EQUALS_B, solver, Y, &cc);

    x.assign((REAL*)X->x, (REAL*)X->x + X->nrow);
   
    cholmod_l_free_dense (&Y, &cc);
    cholmod_l_free_dense (&X, &cc);

    SuiteSparseQR_free (&solver, &cc);
#else
    (*solver)(rhs, x);

    delete solver;
#endif

    /* accumulate joint forces into body force and torque vectors */
    REAL *R = &x[0];
    REAL invstep = 1./step;
    for (std::vector<int>::iterator it = iset.begin(); it != iset.end(); it++, R += 3)
    {
      REAL Rot[9], A[3], a[3], t[3];
      
      int i = *it;

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

#if defined(_OPENMP)
	  omp_set_lock (&lock[part]);
#endif

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
#if defined(_OPENMP)
	  omp_unset_lock (&lock[part]);
#endif
	}
      }
    }
  }

#if defined(_OPENMP)
  for (int i = 0; i < parnum; i ++) omp_destroy_lock (&lock[i]);
#endif
}

} /* namespace */

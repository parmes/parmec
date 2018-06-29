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

typedef amgcl::make_solver<
    amgcl::amg<
	Backend,
	amgcl::coarsening::smoothed_aggregation,
	amgcl::relaxation::spai0
	>,
    amgcl::solver::bicgstab<Backend>
    > joints_solver;

joints_solver *solve = NULL;

/* reset joints matrix after change */
void reset_joints_matrix (int jnum, int *jpart[2], REAL *jpoint[3],
  REAL *position[6], REAL *rotation[9], REAL *inverse[9], REAL invm[])
{
  /* matrix in CRS format */
  std::vector<int> ptr;
  std::vector<int> col;
  std::vector<REAL> val;

  /* adjacency map of genralized inverse inertia */
  std::map<int, std::set<int>> partadj;

  for (int i = 0; i < jnum; i ++)
  {
    int part1 = jpart[0][i], part2 = jpart[1][i];

    partadj[part1].insert(i);
    if (part2 >= 0) partadj[part2].insert(i);
  }

  /* TODO: use OMP pragma to parallelise this loop */

  /* assemble joints matrix */
  for (int i = 0; i < jnum; i ++)
  {
    REAL Jiv[9], Rot[9], im, A[3], Ask[9], Hi[9], Hj[9], tmp0[9], tmp1[9], Mii[9], Mij[9], mass;

    SET9(Mii, 0.0); 

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
	NTMUL (Jiv, Hi, tmp0);
	NNMUL (Hi, tmp0, tmp1);
	NNADD (Mii, tmp1, Mii);
	Mii[0] += im;
	Mii[4] += im;
	Mii[8] += im;

	std::set<int> &adj = partadj[part];
	for (std::set<int>::iterator it = adj.begin(); it != adj.end(); it ++)
	{
	  int j = *it;
	  if (i != j) /* Mij */
	  {
	    A[0] = position[3][part] - jpoint[0][j];
	    A[1] = position[4][part] - jpoint[1][j];
	    A[2] = position[5][part] - jpoint[2][j];
	    VECSKEW (A, Ask);
	    NNMUL (Rot, Ask, Hj);
	    NTMUL (Jiv, Hj, tmp0);
	    NNMUL (Hj, tmp0, Mij);
	    Mij[0] += im;
	    Mij[4] += im;
	    Mij[8] += im;

	    if ((part == jpart[0][i] && part == jpart[1][j]) ||
	        (part == jpart[1][i] && part == jpart[0][j]))
	    {
	      SCALE9(Mij, -1.0);
	    }

	    /* TODO: accumulate Mij into (ptr, col, val) */
	  }
	}
      }
    }

    /* TODO: accumulate Mii into (ptr, col, val) */
  }

  /* solver parameters */
  joints_solver::params prm;
  prm.solver.tol = REAL_EPS;
  prm.solver.maxiter = 100;

  /* problem size */
  int n = ptr.size() - 1;

  /* set up solver */
  if (solve) delete solve;

  solve = new joints_solver(std::tie(n, ptr, col, val), prm);
}

/* solve joints and update forces */
void solve_joints (int jnum, int *jpart[2], REAL *jpoint[3],
  REAL *position[6], REAL *rotation[9], REAL *inverse[9], REAL invm[],
  REAL *linear[3], REAL *angular[6], REAL *force[6], REAL *torque[6], REAL step)
{
  /* assemble joints-free local velocity vector */

  /* TODO */

  /* solve for joint forces */

  /* TODO */

  /* accumulate joint forces into body force and torque vectors */

  /* TODO */
}

} /* namespace */

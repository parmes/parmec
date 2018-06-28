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
  const std::vector<int> ptr;
  const std::vector<int> col;
  const std::vector<REAL> val;

  /* assemble joints matrix */

  /* TODO */

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

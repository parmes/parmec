/*
The MIT License (MIT)

Copyright (c) 2018 Tomasz Koziara

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

#include <vector>
#include <tuple>

#include <amgcl/make_solver.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/adapter/crs_tuple.hpp>

#ifndef __solver__
#define __solver__

void simple_solve (
        /* matrix in CRS format */
        const std::vector<int>    &ptr,
        const std::vector<int>    &col,
        const std::vector<REAL> &val,
        /* right hand side */
        const std::vector<REAL> &rhs,
        /* initial solution on input, solution on output */
        std::vector<REAL> &x,
	/* maximum iterations and solver tolerance */
	int maxiter,
	REAL tol
        )
{
  typedef amgcl::backend::builtin<REAL> Backend;

  typedef amgcl::make_solver<
      amgcl::amg<
	  Backend,
	  amgcl::coarsening::smoothed_aggregation,
	  amgcl::relaxation::spai0
	  >,
      amgcl::solver::bicgstab<Backend>
      > Solver;

  Solver::params prm;

  prm.solver.maxiter = maxiter;
  prm.solver.tol = tol;

  int n = ptr.size() - 1;

  /* stup the solver */
  Solver solve(std::tie(n, ptr, col, val), prm);

  /* solve the problem */
  int    iters;
  REAL error;
  std::tie(iters, error) = solve(rhs, x);
}

#endif

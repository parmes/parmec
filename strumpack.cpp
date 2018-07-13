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

#include <StrumpackSparseSolver.hpp>

#include "strumpack.h"

using namespace strumpack;

StrumpackSparseSolver<float,ptrdiff_t> *float_solver = NULL; 

StrumpackSparseSolver<double,ptrdiff_t> *double_solver = NULL; 

void StrumpackCreate (ptrdiff_t n, ptrdiff_t *ptr, ptrdiff_t *col, float *val)
{
  if (float_solver) delete float_solver;

  float_solver = new StrumpackSparseSolver<float,ptrdiff_t>(false);

  //float_solver->options().set_rel_tol(1E-5);
  //float_solver->options().set_gmres_restart(10);
  //float_solver->options().enable_HSS();
  //float_solver->options().set_Krylov_solver(KrylovSolver::DIRECT);
  float_solver->options().set_matching(MatchingJob::NONE);
  float_solver->set_csr_matrix(n, ptr, col, val);
  float_solver->reorder();
  float_solver->factor();
}

void StrumpackSolve (float *b, float *x)
{
  if (float_solver) float_solver->solve (b, x);
}

void StrumpackCreate (ptrdiff_t n, ptrdiff_t *ptr, ptrdiff_t *col, double *val)
{
  if (double_solver) delete float_solver;

  double_solver = new StrumpackSparseSolver<double,ptrdiff_t>(false);

  //double_solver->options().set_rel_tol(1E-10);
  //double_solver->options().set_gmres_restart(10);
  //double_solver->options().enable_HSS();
  //double_solver->options().set_Krylov_solver(KrylovSolver::DIRECT);
  double_solver->options().set_matching(MatchingJob::NONE);
  double_solver->set_csr_matrix(n, ptr, col, val);
  double_solver->reorder();
  double_solver->factor();
}

void StrumpackSolve (double *b, double *x)
{
  if (double_solver) double_solver->solve (b, x);
}

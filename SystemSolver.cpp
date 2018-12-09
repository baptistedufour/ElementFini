#ifndef _SYSTEMSOLVER_CPP

#include "SystemSolver.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>

using namespace std;
using namespace Eigen;

SystemSolver::SystemSolver(SystemBuilder* builder)
{
  _builder = builder;
  _nb_pts = builder->Get_nb_pts();
}

void SystemSolver::BuildSol()
{
  cout << "Creation du solveur lineaire" << endl;
  SparseMatrix<double, ColMajor> matTemp;
  matTemp = _builder->Get_matA();
  cout << "-------------------------------" << endl;
  _solverMethod.analyzePattern(matTemp);
  _solverMethod.factorize(matTemp);

  cout << "-------------------------------" << endl;
  VectorXd sourceTermTemp(_nb_pts);
  sourceTermTemp = MatrixXd(_builder->Get_sourceTerm());

  cout << "-------------------------------" << endl;
  VectorXd solTemp(_nb_pts);
  solTemp = _solverMethod.solve(sourceTermTemp);

  cout << "-------------------------------" << endl;
  _sol.resize(_nb_pts);
  _sol.setZero();
  _sol = solTemp.sparseView();
  cout << "-------------------------------" << endl;
}

// Sauvegarde la solution
void SystemSolver::SaveSol()
{
	string name_file = "Resultats/sol.dat";

  assert((_sol.size() == _nb_pts) && "The size of the solution vector is not the same than the number of 2 * nb_pts !");

	ofstream solution;
	solution.open(name_file, ios::out);
	solution.precision(7);

	solution.close();
}


const VectorXd SystemSolver::ExactSolution() const
{
  if ((_sigma_choice != constant)||(_source_fct_choice != constant))
  {
    std::cout << "The exact solution is not known !" << std::endl;
    abort();
  }
  VectorXd vec;
  vec.resize(10);
  vec.setZero();
  return vec;
}


double Error(Eigen::SparseVector<double> sol)
{
  return 0;
}

#define _SYSTEMSOLVER_CPP
#endif

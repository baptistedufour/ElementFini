#ifndef _SYSTEMSOLVER_CPP

#include "SystemSolver.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>

using namespace std;
using namespace Eigen;

SystemSolver::SystemSolver(SystemBuilder* builder) :
_source_fct_choice(builder->Get_source_fct_choice()),
_ur(builder->Get_right_BC()),_ul(builder->Get_left_BC()),
_sigma_choice(builder->Get_sigma_choice())
{
  _builder = builder;
  _nb_pts = builder->Get_nb_pts();
}

void SystemSolver::BuildSol()
{
  cout << "-------------------------------------------------" << endl;
  cout << "Résolution du systeme lineaire" << endl;
  SparseMatrix<double, ColMajor> matTemp;
  matTemp = _builder->Get_matA();

  _solverMethod.analyzePattern(matTemp);
  _solverMethod.factorize(matTemp);

  VectorXd sourceTermTemp(_nb_pts);
  sourceTermTemp = MatrixXd(_builder->Get_sourceTerm());

  VectorXd solTemp(_nb_pts);
  solTemp = _solverMethod.solve(sourceTermTemp);
  //cout << solTemp << endl;

  _sol.resize(_nb_pts);
  _sol.setZero();
  _sol = solTemp.sparseView();

  //cout << "Test A*solApprox = " << endl << matTemp*_sol << endl;
}

// Sauvegarde la solution
void SystemSolver::SaveSol()
{
  cout << "-------------------------------------------------" << endl;
  cout << "Sauvegarde de la solution " << endl;
	string name_file = "Resultats/sol.dat";
  string name_file2 = "Resultats/solbis.dat";

  assert((_sol.size() == 2*(_nb_pts+1)) && "The size of the solution vector is not the same than the number of 2 * nb_pts !");
  double dx = 1./(_nb_pts+1);

	ofstream solution;
	solution.open(name_file, ios::out);
	solution.precision(7);

  ofstream solution2;
  solution2.open(name_file2, ios::out);
  solution2.precision(7);
  double alpha = (1.+sqrt(3))/2., beta = (1.-sqrt(3))/2.;

  solution << 0 << " " << _ul << endl;
  solution2 << 0 << " " << _ul << endl;
  for(int i=0; i<_nb_pts+1;i++)
  {
      solution << i*dx + dx*(sqrt(3)-1)/(2*sqrt(3)) << " " << _sol.coeffRef(2*i) << endl;
      solution << i*dx + dx*(sqrt(3)+1)/(2*sqrt(3)) << " " << _sol.coeffRef(2*i+1) << endl;

      solution2 << i*dx + 1e-6 << " " << alpha*_sol.coeffRef(2*i) + beta*_sol.coeffRef(2*i+1) << endl;
      solution2 << (i+1)*dx - 1e-6 << " " << beta*_sol.coeffRef(2*i) + alpha*_sol.coeffRef(2*i+1) << endl;
  }
  solution2 << 1 << " " << _ur << endl;
  solution2.close();

  solution << 1 << " " << _ur << endl;
	solution.close();
  cout << "-------------------------------------------------" << endl;
}


const VectorXd SystemSolver::ExactSolution() const
{
  if ((_sigma_choice != "constant")||(_source_fct_choice != "constant"))
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

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
}

// Sauvegarde la solution
void SystemSolver::SaveSol()
{
  cout << "-------------------------------------------------" << endl;
  cout << "Sauvegarde de la solution " << endl;
	string name_file = "Resultats/sol.dat";

  assert((_sol.size() == 2*_nb_pts) && "The size of the solution vector is not the same than the number of 2 * nb_pts !");
  double dx = 1./(_nb_pts+1);

	ofstream solution;
	solution.open(name_file, ios::out);
	solution.precision(7);
  for(int i=0; i<_nb_pts+1;i++)
  {
    if ((i!=0)&&(i!=_nb_pts))
    {
      solution << i*dx + dx*(sqrt(3)-1)/(2*sqrt(3)) << " " << _sol.coeffRef(2*i-1) << endl;
      solution << i*dx + dx*(sqrt(3)+1)/(2*sqrt(3)) << " " << _sol.coeffRef(2*i) << endl;
    }
    else if (i==0)
    {
      solution << dx/2. << " " << _sol.coeffRef(0) << endl;
    }
    else
    {
      solution << 1-dx/2. << " " << _sol.coeffRef(2*_nb_pts-1) << endl;
    }
  }
	solution.close();
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

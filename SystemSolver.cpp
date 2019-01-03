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
  _sol.resize(2*(_nb_pts+1));
  _sol.setZero();
  _solEx.resize(2*(_nb_pts+1));
  _solEx.setZero();
  _errorL = 0;
  _errorH = 0;

  _dx = 1./(_nb_pts+1);
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

  _sol = solTemp.sparseView();
}

// Sauvegarde la solution
void SystemSolver::SaveSol()
{
  cout << "-------------------------------------------------" << endl;
  cout << "Sauvegarde de la solution " << endl;
	string name_file = "Resultats/sol.dat";
  string name_file2 = "Resultats/solbis.dat";

  bool checkSolEx = false;
  if (((_builder->Get_source_fct_choice() == "constant")||(_builder->Get_source_fct_choice() == "line"))
    &&(_builder->Get_sigma_choice() == "constant"))
  {
    checkSolEx = true;
    BuildSolEx();
  }

  assert((_sol.size() == 2*(_nb_pts+1)) && "The size of the solution vector is not the same than the number of 2 * nb_pts !");

	ofstream solution;
	solution.open(name_file, ios::out);
	solution.precision(7);

  ofstream solution2;
  solution2.open(name_file2, ios::out);
  solution2.precision(7);
  double alpha = (1.+sqrt(3))/2., beta = (1.-sqrt(3))/2.;

  solution << 0 << " " << _ul;
  if (checkSolEx) solution << " " << _ul;
  solution << endl;
  solution2 << 0 << " " << _ul;
  if (checkSolEx) solution2 << " " << _ul;
  solution2 << endl;
  for(int i=0; i<_nb_pts+1;i++)
  {
      solution << i*_dx + _dx*(sqrt(3)-1)/(2*sqrt(3)) << " " << _sol.coeffRef(2*i);
      if (checkSolEx) solution << " " << _solEx.coeffRef(2*i);
      solution << endl;
      solution << i*_dx + _dx*(sqrt(3)+1)/(2*sqrt(3)) << " " << _sol.coeffRef(2*i+1);
      if (checkSolEx) solution << " " << _solEx.coeffRef(2*i+1);
      solution << endl;

      solution2 << i*_dx + 1e-6 << " " << alpha*_sol.coeffRef(2*i) + beta*_sol.coeffRef(2*i+1);
      if (checkSolEx) solution2 << " " << alpha*_solEx.coeffRef(2*i) + beta*_solEx.coeffRef(2*i+1);
      solution2 << endl;
      solution2 << (i+1)*_dx - 1e-6 << " " << beta*_sol.coeffRef(2*i) + alpha*_sol.coeffRef(2*i+1);
      if (checkSolEx) solution2 << " " << beta*_solEx.coeffRef(2*i) + alpha*_solEx.coeffRef(2*i+1);
      solution2 << endl;
  }
  solution2 << 1 << " " << _ur;
  if (checkSolEx) solution2 << " " << _ur;
  solution2 << endl;
  solution2.close();

  solution << 1 << " " << _ur;
  if (checkSolEx) solution << " " << _ur;
  solution << endl;
	solution.close();
  cout << "-------------------------------------------------" << endl;
}

void SystemSolver::BuildSolEx()
{
  double xk, xk1;
  double a,b,c,d;
  double sigma = _builder->Get_sigma0();
  if(_builder->Get_source_fct_choice() == "line")
  {
    a = _builder->Get_param_d();
    b = _builder->Get_param_e();
    c = (_ur-_ul)+a/(6.*sigma)+b/(2.*sigma);
    d = _ul;
    a = -a/(6.*sigma);
    b = -b/(2.*sigma);

    for(int i = 0; i < _nb_pts+1; i++)
    {
      xk  = i*_dx + _dx*(sqrt(3)-1.)/(2*sqrt(3));
      xk1 = i*_dx + _dx*(sqrt(3)+1.)/(2*sqrt(3));

      _solEx.coeffRef(2*i)   = a*xk*xk*xk + b*xk*xk + c*xk + d;
      _solEx.coeffRef(2*i+1) = a*xk1*xk1*xk1 + b*xk1*xk1 + c*xk1 + d;
    }
  }
  else if (_builder->Get_source_fct_choice() == "line")
  {
    a = _builder->Get_param_d();
    b = (_ur-_ul)+b/(2.*sigma);
    c = _ul;
    a = -a/(2.*sigma);

    for(int i = 0; i < _nb_pts+1; i++)
    {
      xk  = i*_dx + _dx*(sqrt(3)-1)/(2*sqrt(3));
      xk1 = i*_dx + _dx*(sqrt(3)+1)/(2*sqrt(3));

      _solEx.coeffRef(2*i)   = a*xk*xk + b*xk + c;
      _solEx.coeffRef(2*i+1) = a*xk1*xk1 + b*xk1 + c;
    }
  }
}

void SystemSolver::ErrorLinf()
{
  double error;

  VectorXd exact_sol  = MatrixXd(_solEx);
  VectorXd approx_sol = MatrixXd(_sol);

  error = ((exact_sol-approx_sol).array().abs()).maxCoeff();

  cout << " -- Norme Linf : " << error << endl;
}

void SystemSolver::ErrorL2()
{
  _errorL = 0;
  for(int i=0; i<_nb_pts+1;i++)
  {
    _errorL += _errorL;
  }
  cout << " -- Norme L2   : indisponible" << endl;
}

void SystemSolver::ErrorH1()
{
  cout << " -- Norme H1   : indisponible" << endl;
}

#define _SYSTEMSOLVER_CPP
#endif

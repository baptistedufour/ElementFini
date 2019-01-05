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
  if (((_builder->Get_source_fct_choice() == "constant")||(_builder->Get_source_fct_choice() == "line")||(_builder->Get_source_fct_choice() == "sinus"))
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
  else if (_builder->Get_source_fct_choice() == "constant")
  {
    a = _builder->Get_param_d();
    b = (_ur-_ul)+a/(2.*sigma);
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
  else if (_builder->Get_source_fct_choice() == "sinus")
  {
    for(int i = 0; i < _nb_pts+1; i++)
    {
      xk  = i*_dx + _dx*(sqrt(3)-1)/(2*sqrt(3));
      xk1 = i*_dx + _dx*(sqrt(3)+1)/(2*sqrt(3));

      _solEx.coeffRef(2*i)   = sin(PI * xk);
      _solEx.coeffRef(2*i+1) = sin(PI * xk1);
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

void SystemSolver::ErrorL2_sin()
{
  _errorL = 0;
  double Ck, alphk, alphk1, xk, xk1;
  SparseVector<double> sigVec = _builder->Get_sigma();


  for(int i=0; i<_nb_pts+1;i++)
  {
    xk     =  i*_dx;
    xk1    =  (i+1)*_dx;
    alphk  =  _dx/2.*(1.-1./sqrt(3)) + xk;
    alphk1 =  _dx/2.*(1.+1./sqrt(3)) + xk;
    Ck= (_sol.coeffRef(2*i+1) - _sol.coeffRef(2*i))/(alphk1 - alphk);

    _errorL += pow(Ck,2)/3. * (pow(xk1,3)-pow(xk,3));
    _errorL += (_sol.coeffRef(2*i) - Ck*alphk)*( Ck*(pow(xk1,2)-pow(xk,2)) + 2./(PI*sigVec.coeffRef(i))*(cos(PI*xk1)-cos(PI*xk)) );
    _errorL += 2.*Ck/(PI*sigVec.coeffRef(i))*( (xk1*cos(PI*xk1)-xk*cos(PI*xk)) - (sin(PI*xk1)-sin(PI*xk))/PI );
    _errorL += (1./2. *(xk1-xk) - 1./(4.*PI) *(sin(2*PI*xk1)-sin(2*PI*xk)))/pow(sigVec.coeffRef(i),2);
    _errorL += pow(Ck*alphk - _sol.coeffRef(2*i),2)*(xk1-xk);

  }
  cout << " -- Norme L2   : " << _errorL << endl;
}


void SystemSolver::ErrorL2()
{
  _errorL = 0;
  double Ck, alphk, alphk1, xk, xk1;
  double a, b, c, d;
  double sigma = _builder->Get_sigma0();
  if(_builder->Get_source_fct_choice() == "line")
  {
    a = _builder->Get_param_d();
    b = _builder->Get_param_e();
    c = (_ur-_ul)+a/(6.*sigma)+b/(2.*sigma);
    d = _ul;
    a = -a/(6.*sigma);
    b = -b/(2.*sigma);
    //cout << " -- poly  line: " << a << "   "<< b << "   "<< c << "   "<< d << "   "  << endl;

  }
  else if (_builder->Get_source_fct_choice() == "constant")
  {
    a = _builder->Get_param_d();
    b = -a/(2.*sigma);
    c = (_ur-_ul)+a/(2.*sigma);
    d = _ul;
    a = 0.;
    //cout << " -- poly  cst: " << a << "   "<< b << "   "<< c << "   "<< d << "   "  << endl;

  }
  for(int i=0; i<_nb_pts+1;i++)
  {
    xk     =  i*_dx;
    xk1    =  (i+1)*_dx;
    alphk  =  _dx/2.*(1.-1./sqrt(3)) + xk;
    alphk1 =  _dx/2.*(1.+1./sqrt(3)) + xk;
    Ck= (_sol.coeffRef(2*i+1) - _sol.coeffRef(2*i))/(alphk1 - alphk);



    _errorL += ( a*a )                                                                                                           *(pow(xk1,7)-pow(xk,7))/7.;
    _errorL += ( 2.*a*b )                                                                                                        *(pow(xk1,6)-pow(xk,6))/6.;
    _errorL += ( b*b + 2.*a*c - 2.*Ck*a )                                                                                        *(pow(xk1,5)-pow(xk,5))/5.;
    _errorL += 2.*(a*d + b*c - Ck*b + (Ck*alphk-_sol.coeffRef(2*i))*a )                                                          *(pow(xk1,4)-pow(xk,4))/4.;
    _errorL += (c*c + 2.*b*d - 2.*Ck*c + Ck*Ck + 2.*(Ck*alphk-_sol.coeffRef(2*i))*b )                                            *(pow(xk1,3)-pow(xk,3))/3.;
    _errorL += 2.*( c*d - Ck*d + Ck*_sol.coeffRef(2*i) - Ck*Ck*alphk + (Ck*alphk-_sol.coeffRef(2*i))*c )                         *(pow(xk1,2)-pow(xk,2))/2.;
    _errorL += ( d*d + Ck*Ck*alphk*alphk - 2.*Ck*_sol.coeffRef(2*i)*alphk + _sol.coeffRef(2*i)*_sol.coeffRef(2*i) + 2.*(Ck*alphk-_sol.coeffRef(2*i))*d)*_dx;


  }
  cout << " -- Norme L2   : " << sqrt(abs(_errorL)) << endl;
}

void SystemSolver::ErrorH1()
{
  cout << " -- Norme H1   : indisponible" << endl;
}

#define _SYSTEMSOLVER_CPP
#endif

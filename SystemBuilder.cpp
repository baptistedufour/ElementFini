#ifndef _SYSTEMBUILDER_CPP

#include "SystemBuilder.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <chrono>

using namespace Eigen;
using namespace std;

SystemBuilder::SystemBuilder(DataFile* data_file) :
_nb_pts(data_file->Get_Nmesh()), _src_choice(data_file->Get_source_fct_choice()),
BC_right(data_file->Get_right_BC_choice()),_ur(data_file->Get_right_BC()),
BC_left(data_file->Get_left_BC_choice()),_ul(data_file->Get_left_BC()),
_gamma0(data_file->Get_param_gamma0()), _d(data_file->Get_param_d()),
_Pk_choice(data_file->Get_Pk_choice()), _a(data_file->Get_param_a()),
_sigma_choice(data_file->Get_sigma_choice())
{
  //Parametre
  _dx = 1./(_nb_pts+1.);
  _Asize = 2*(_nb_pts+1);

  //sigma
  if(_sigma_choice == "line")
  {
    _b = data_file->Get_param_b();
  }
  else if (_sigma_choice == "creneau")
  {
    _b = data_file->Get_param_b();
    _c = data_file->Get_param_c();
  }

  _sigma.resize(_Asize);
  double xk, xk1;
  for(int i = 0; i < _nb_pts+1; i++)
  {
    if(_sigma_choice == "constant")
    {
      _sigma.coeffRef(2*i) = _a;
      _sigma.coeffRef(2*i+1) = _a;
    }
    else
    {
      xk  = i*_dx + _dx*(sqrt(3)-1)/(2*sqrt(3));
      xk1 = i*_dx + _dx*(sqrt(3)+1)/(2*sqrt(3));
      if(_sigma_choice == "creneau")
      {
        if ((xk<=_b)&&(xk>=_a))
        {
          _sigma.coeffRef(2*i) = _c;
        }
        else
        {
          _sigma.coeffRef(2*i) = 1.;
        }

        if ((xk1<=_b)&&(xk>=_a))
        {
          _sigma.coeffRef(2*i+1) = _c;
        }
        else
        {
          _sigma.coeffRef(2*i+1) = 1.;
        }
      }
      else if(_sigma_choice == "line")
      {
        _sigma.coeffRef(2*i)=xk;
        _sigma.coeffRef(2*i+1)=xk1;
      }
      else
      {
        _sigma.coeffRef(2*i) = 1.;
        _sigma.coeffRef(2*i+1) = 1.;
      }
    }
  }

  _matAm.resize(_Asize,_Asize);
  _matAm.setZero();
  _matAs.resize(_Asize,_Asize);
  _matAs.setZero();
  _matA.resize(_Asize,_Asize);
  _matA.setZero();


  //Terme source
  if(_src_choice == "line")
  {
    _e = data_file->Get_param_e();
  }
  else if (_src_choice == "creneau")
  {
    _e = data_file->Get_param_e();
    _f = data_file->Get_param_f();
  }
}


void SystemBuilder::Build_matA()
{
  cout << "Création de A (Taille " << _Asize << ") : " << endl;

  double alpha_m = 1+sqrt(3)/2.;
  double beta_m  = 1-sqrt(3)/2.;
  double alpha_s = (3+sqrt(3))/(2.*_dx);
  double beta_s  = (3-sqrt(3))/(2.*_dx);
  double gamma_s = 3/(2.*_dx);
  double mu =_gamma0/_dx;
  double sigma = _sigma.coeffRef(0);

  for (int i = 0; i < _nb_pts+1; i++)
  {
     cout.flush();
     cout << "Progression : " << (double)i/((double)(_nb_pts))*100 << "% \r";

     if( (i > 0) && (i < _nb_pts) )
     {
       _matAm.coeffRef(2*i,2*i-2) += 0.5;
       _matAm.coeffRef(2*i,2*i-1) += -alpha_m;
       _matAm.coeffRef(2*i,2*i)   += 2;
       _matAm.coeffRef(2*i,2*i+1) += -1;
       _matAm.coeffRef(2*i,2*i+2) += 0.5;
       _matAm.coeffRef(2*i,2*i+3) += -beta_m;

       _matAm.coeffRef(2*i+1,2*i-2) += -beta_m;
       _matAm.coeffRef(2*i+1,2*i-1) += 0.5;
       _matAm.coeffRef(2*i+1,2*i)   += -1;
       _matAm.coeffRef(2*i+1,2*i+1) += 2;
       _matAm.coeffRef(2*i+1,2*i+2) += -alpha_m;
       _matAm.coeffRef(2*i+1,2*i+3) += 0.5;


       _matAs.coeffRef(2*i,2*i-2) += -gamma_s;
       _matAs.coeffRef(2*i,2*i-1) += alpha_s;
       _matAs.coeffRef(2*i,2*i+2) += -gamma_s;
       _matAs.coeffRef(2*i,2*i+3) += beta_s;

       _matAs.coeffRef(2*i+1,2*i-2) += beta_s;
       _matAs.coeffRef(2*i+1,2*i-1) += -gamma_s;
       _matAs.coeffRef(2*i+1,2*i+2) += alpha_s;
       _matAs.coeffRef(2*i+1,2*i+3) += -gamma_s;
     }

     else if (i==0)
     {
       _matAm.coeffRef(2*i,2*i)   += 2;
       _matAm.coeffRef(2*i,2*i+1) += -1;
       _matAm.coeffRef(2*i,2*i+2) += 0.5;
       _matAm.coeffRef(2*i,2*i+3) += -beta_m;

       _matAm.coeffRef(2*i+1,2*i)   += -1;
       _matAm.coeffRef(2*i+1,2*i+1) += 2;
       _matAm.coeffRef(2*i+1,2*i+2) += -alpha_m;
       _matAm.coeffRef(2*i+1,2*i+3) += 0.5;


       _matAs.coeffRef(2*i,2*i)   += -alpha_s;
       _matAs.coeffRef(2*i,2*i+1) += gamma_s;
       _matAs.coeffRef(2*i,2*i+2) += -gamma_s;
       _matAs.coeffRef(2*i,2*i+3) += beta_s;

       _matAs.coeffRef(2*i+1,2*i)   += gamma_s;
       _matAs.coeffRef(2*i+1,2*i+1) += -beta_s;
       _matAs.coeffRef(2*i+1,2*i+2) += alpha_s;
       _matAs.coeffRef(2*i+1,2*i+3) += -gamma_s;
     }

     else
     {
       _matAm.coeffRef(2*i,2*i-2) += 0.5;
       _matAm.coeffRef(2*i,2*i-1) += -alpha_m;
       _matAm.coeffRef(2*i,2*i)   += 2;
       _matAm.coeffRef(2*i,2*i+1) += -1;

       _matAm.coeffRef(2*i+1,2*i-2) += -beta_m;
       _matAm.coeffRef(2*i+1,2*i-1) += 0.5;
       _matAm.coeffRef(2*i+1,2*i)   += -1;
       _matAm.coeffRef(2*i+1,2*i+1) += 2;


       _matAs.coeffRef(2*i,2*i-2) += -gamma_s;
       _matAs.coeffRef(2*i,2*i-1) += alpha_s;
       _matAs.coeffRef(2*i,2*i)   += -beta_s;
       _matAs.coeffRef(2*i,2*i+1) += gamma_s;

       _matAs.coeffRef(2*i+1,2*i-2) += beta_s;
       _matAs.coeffRef(2*i+1,2*i-1) += -gamma_s;
       _matAs.coeffRef(2*i+1,2*i)   += gamma_s;
       _matAs.coeffRef(2*i+1,2*i+1) += -alpha_s;
     }
  }

  _matAm = mu*_matAm;
  //_matAs = sigma*matA;
  for(int i=0; i < _Asize; i++)
  {
    for(int j=0; j < _Asize; j++)
    {
      _matAs.coeffRef(i,j) = _matAs.coeffRef(i,j)*_sigma.coeffRef(i);
    }
  }
  _matA = _matAm + _matAs;

  Matrix<double, Dynamic, Dynamic> Matrix;
  Matrix = MatrixXd(_matA);
  //cout << Matrix << endl;
  //cout << endl << " --- " << endl << endl;

  SparseMatrix<double, ColMajor> A(_matA);

  cout << "Factorisation LU                 " << endl;
  _solverMethod.analyzePattern(A);
  _solverMethod.factorize(A);
}

void SystemBuilder::Build_sourceTerm()
{
  cout << "Création de b (Taille " << _Asize << ") : " << endl;
  _sourceTerm.resize(_Asize);
  _sourceTerm.setZero();

  double sigma = _sigma.coeffRef(0);
  double mu =_gamma0/_dx;
  double xk, xk1;
  for (int i = 0; i < _nb_pts+1; i++)
  {
    cout.flush();
    cout << "Progression : " << (double)i/((double)(_nb_pts))*100 << "% \r";

    if(_src_choice == "constant")
    {
      _sourceTerm.coeffRef(2*i)   = _d;
      _sourceTerm.coeffRef(2*i+1) = _d;

      _sourceTerm.coeffRef(2*i)   *= _dx*0.5;
      _sourceTerm.coeffRef(2*i+1) *= _dx*0.5;
    }
    else if( _src_choice == "line")
    {
      xk  = i;
      xk1 = i+1;
      _sourceTerm.coeffRef(2*i)   = _d*( _dx*sqrt(3)/3.+(1+sqrt(3))*xk/2.+(1-sqrt(3))*xk1/2.) + _e;
      _sourceTerm.coeffRef(2*i+1) = _d*(-_dx*sqrt(3)/3.+(1-sqrt(3))*xk/2.+(1+sqrt(3))*xk1/2.) + _e;

      _sourceTerm.coeffRef(2*i)   *= _dx*0.5;
      _sourceTerm.coeffRef(2*i+1) *= _dx*0.5;
    }
    else if( _src_choice == "creneau")
    {
      xk  = i*_dx + _dx*(sqrt(3)-1)/(2*sqrt(3));
      xk1 = i*_dx + _dx*(sqrt(3)+1)/(2*sqrt(3));

      if((xk>=_d)&&(xk<=_e))
        _sourceTerm.coeffRef(2*i) = _f;

      if((xk1>=_d)&&(xk1<=_e))
        _sourceTerm.coeffRef(2*i+1) = _f;

      _sourceTerm.coeffRef(2*i)   *= _dx*0.5;
      _sourceTerm.coeffRef(2*i+1) *= _dx*0.5;
    }
    else
    {
      _sourceTerm.coeffRef(2*i)   = 0;
      _sourceTerm.coeffRef(2*i+1) = 0;
    }
  }

  //Condition aux bords
  if(BC_left == "dirichlet")
  {
    _sourceTerm.coeffRef(0) += -_sigma.coeffRef(0)*sqrt(3)*_ul/_dx;
    _sourceTerm.coeffRef(1) += _sigma.coeffRef(1)*sqrt(3)*_ul/_dx;

    _sourceTerm.coeffRef(0) += mu*(1.+sqrt(3))*_ul/2.;
    _sourceTerm.coeffRef(1) += mu*(1.-sqrt(3))*_ul/2.;
  }
  else
  {

    _sourceTerm.coeffRef(0) -= mu*(_ul*_ul);
    _sourceTerm.coeffRef(1) -= mu*(_ul*_ul);
  }

  if(BC_right == "dirichlet")
  {
    _sourceTerm.coeffRef(_Asize-2) += _sigma.coeffRef(_Asize-2)*sqrt(3)*_ur/_dx;
    _sourceTerm.coeffRef(_Asize-1) += -_sigma.coeffRef(_Asize-1)*sqrt(3)*_ur/_dx;

    _sourceTerm.coeffRef(_Asize-2) += mu*(1-sqrt(3))*_ur/2.;
    _sourceTerm.coeffRef(_Asize-1) += mu*(1+sqrt(3))*_ur/2.;
  }
  else
  {
    _sourceTerm.coeffRef(_Asize-2) -= mu*(_ur*_ur);
    _sourceTerm.coeffRef(_Asize-1) -= mu*(_ur*_ur);
  }
  //cout << _sourceTerm << endl;
}


#define _SYSTEMBUILDER_CPP
#endif

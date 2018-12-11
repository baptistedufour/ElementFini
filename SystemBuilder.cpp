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
  _dx = 1./(_nb_pts+1);
  _Asize = 2*_nb_pts;

  //Conditions aux bords (neumann, dirichlet u_g)
  if(BC_right == "dirichlet")
  {
    _aR = 2; _bR = -1;
  }
  else
  {
    _aR = 1; _bR = -_dx/2.;
  }

  if(BC_left == "dirichlet")
  {
    _aL = 2; _bL = -1;
  }
  else
  {
    _aL = 1; _bL = -_dx/2.;
  }

  //Terme source
  if( _src_choice == "line" )
  {
    _b = data_file->Get_param_b();
  }
  else if( _src_choice == "curve" )
  {
    _b = data_file->Get_param_b();
    _c = data_file->Get_param_b();
  }

  //system("mkdir -p ./Results");
}


void SystemBuilder::Build_matA()
{
  cout << "Création de A (Taille " << _Asize << ") : " << endl;
  _matA.resize(_Asize,_Asize);
  _matA.setZero();

  double alpha = sqrt(3)/2.+1;
  double beta = 1-sqrt(3)/2.;
  double mu =_gamma0/_dx;

  for (int i = 0; i < _nb_pts+1; i++)
  {
     cout.flush();
     cout << "Progression : " << (double)i/((double)(_nb_pts))*100 << "% \r";

     if( (i > 1) && (i < _nb_pts-1) )
     {
       _matA.coeffRef(2*i-1,2*i-3) = 0.5;
       _matA.coeffRef(2*i-1,2*i-2) = -alpha;
       _matA.coeffRef(2*i-1,2*i-1) = 2;
       _matA.coeffRef(2*i-1,2*i)   = -1;
       _matA.coeffRef(2*i-1,2*i+1) = 0.5;
       _matA.coeffRef(2*i-1,2*i+2) = -beta;

       _matA.coeffRef(2*i,2*i-3) = -beta;
       _matA.coeffRef(2*i,2*i-2) = 0.5;
       _matA.coeffRef(2*i,2*i-1) = -1;
       _matA.coeffRef(2*i,2*i)   = 2;
       _matA.coeffRef(2*i,2*i+1) = -alpha;
       _matA.coeffRef(2*i,2*i+2) = 0.5;
     }

     else if (i==1)
     {
       _matA.coeffRef(2*i-1,2*i-1) = 2;
       _matA.coeffRef(2*i-1,2*i)   = -1;
       _matA.coeffRef(2*i-1,2*i+1) = 0.5;
       _matA.coeffRef(2*i-1,2*i+2) = -beta;

       _matA.coeffRef(2*i,2*i-1) = -1;
       _matA.coeffRef(2*i,2*i)   = 2;
       _matA.coeffRef(2*i,2*i+1) = -alpha;
       _matA.coeffRef(2*i,2*i+2) = 0.5;

       if(BC_left == "dirichlet")
       {
         _matA.coeffRef(2*i-1,2*i-2) = -_aL*(1+sqrt(3))/2.;
         _matA.coeffRef(2*i  ,2*i-2) = -_aL*(1-sqrt(3))/2.;
       }
       else
       {
         _matA.coeffRef(2*i-1,2*i-2) = 0;
         _matA.coeffRef(2*i  ,2*i-2) = 0;
       }
     }

     else if (i==_nb_pts-1)
     {
       _matA.coeffRef(2*i-1,2*i-3) = 0.5;
       _matA.coeffRef(2*i-1,2*i-2) = -alpha;
       _matA.coeffRef(2*i-1,2*i-1) = 2;
       _matA.coeffRef(2*i-1,2*i)   = -1;

       _matA.coeffRef(2*i,2*i-3) = -beta;
       _matA.coeffRef(2*i,2*i-2) = 0.5;
       _matA.coeffRef(2*i,2*i-1) = -1;
       _matA.coeffRef(2*i,2*i)   = 2;

       if(BC_right == "dirichlet")
       {
         _matA.coeffRef(2*i-1,2*i+1) = -_aR*(1-sqrt(3))/2.;
         _matA.coeffRef(2*i  ,2*i+1) = -_aR*(1+sqrt(3))/2.;
       }
       else
       {
         _matA.coeffRef(2*i-1,2*i+1) = 0;
         _matA.coeffRef(2*i  ,2*i+1) = 0;
       }
     }

     else if (i==0)
     {
       _matA.coeffRef(2*i,2*i)   = 2*_aL;
       _matA.coeffRef(2*i,2*i+1) = -(1+sqrt(3));
       _matA.coeffRef(2*i,2*i+2) = -(1-sqrt(3));
     }

     else
     {
       _matA.coeffRef(2*i-1,2*i-3) = -(1-sqrt(3));
       _matA.coeffRef(2*i-1,2*i-2) = -(1+sqrt(3));
       _matA.coeffRef(2*i-1,2*i-1) = 2*_aR;
     }
  }

  //Affichage
  Matrix<double, Dynamic, Dynamic> Matrix;
  Matrix = MatrixXd(_matA);
  cout << Matrix << endl;

  _matA = mu*_matA;
  SparseMatrix<double, ColMajor> A(_matA);

  cout << "Utilisation de la librairie Eigen" << endl;
  _solverMethod.analyzePattern(A);
  _solverMethod.factorize(A);
}

void SystemBuilder::Build_sourceTerm()
{
  cout << "Création de b (Taille " << _Asize << ") " << endl;
  _sourceTerm.resize(_Asize);
  _sourceTerm.setZero();


  double mu =_gamma0/_dx;
  for (int i = 0; i < _Asize; i++)
  {
    cout.flush();
    cout << "Progression : " << (double)i/((double)(_Asize))*100 << "% \r";

    if(_src_choice == "constant")
    {
      _sourceTerm.coeffRef(i)=_d;
    }
    else if( _src_choice == "line")
    {
      _sourceTerm.coeffRef(i)=_d*i*_dx + _e;
    }
    else
    {
      _sourceTerm.coeffRef(i)=_d*i*_dx*i*_dx + _e*i*_dx + _f;
    }
  }

  cout << "ul = " << _ul << "  ; ur = " << _ur <<  "          " << endl;
  cout << "bL = " << _bL << " ; bR = " << _bR <<  "          " << endl;

  //Condition aux bords
  if(BC_left == "dirichlet")
  {
    _sourceTerm.coeffRef(0) = _sourceTerm.coeffRef(0)-mu*2*_bL*_ul;
    _sourceTerm.coeffRef(1) = _sourceTerm.coeffRef(1)+mu*(1.+sqrt(3))*_bL*_ul/2.;
    _sourceTerm.coeffRef(2) = _sourceTerm.coeffRef(2)+mu*(1.-sqrt(3))*_bL*_ul/2.;
  }
  else
  {
    _sourceTerm.coeffRef(0) = _sourceTerm.coeffRef(0)-mu*2*_bL*_ul;
    _sourceTerm.coeffRef(1) = _sourceTerm.coeffRef(1)+mu*(1.-sqrt(3))*_ul/2.;
    _sourceTerm.coeffRef(2) = _sourceTerm.coeffRef(2)+mu*(1.+sqrt(3))*_ul/2.;
  }

  if(BC_right == "dirichlet")
  {
    _sourceTerm.coeffRef(_Asize-3) = _sourceTerm.coeffRef(_Asize-3)+mu*(1-sqrt(3))*_bR*_ur/2.;
    _sourceTerm.coeffRef(_Asize-2) = _sourceTerm.coeffRef(_Asize-2)+mu*(1+sqrt(3))*_bR*_ur/2.;
    _sourceTerm.coeffRef(_Asize-1) = _sourceTerm.coeffRef(_Asize-1)-mu*2*_bR*_ur;
  }
  else
  {
    _sourceTerm.coeffRef(_Asize-3) = _sourceTerm.coeffRef(_Asize-3)+mu*(1-sqrt(3))*_bR*_ur/2.;
    _sourceTerm.coeffRef(_Asize-2) = _sourceTerm.coeffRef(_Asize-2)+mu*(1+sqrt(3))*_bR*_ur/2.;
    _sourceTerm.coeffRef(_Asize-1) = _sourceTerm.coeffRef(_Asize-1)-mu*2*_bR*_ur;
  }
  cout << _sourceTerm << endl;
}


#define _SYSTEMBUILDER_CPP
#endif

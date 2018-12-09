#ifndef _SYSTEMBUILDER_CPP

#include "SystemBuilder.h"
#include <iostream>
#include <fstream>
#include <math.h>
#include <chrono>

using namespace Eigen;
using namespace std;

SystemBuilder::SystemBuilder(DataFile* data_file) :
_ur(data_file->Get_right_BC()), _ul(data_file->Get_left_BC()),
_a(data_file->Get_param_a()), _d(data_file->Get_param_d()),
_Pk_choice(data_file->Get_Pk_choice())
{
  if(_Pk_choice==1)
  {
    _nb_pts = 2*data_file->Get_Nmesh();
  }
  else
  {
    _nb_pts = 3*data_file->Get_Nmesh()+1;
  }
  
  //system("mkdir -p ./Results");
}


void SystemBuilder::Build_matA()
{
  _matA.resize(_nb_pts,_nb_pts);
  _matA.setZero();

  cout << "Création de A : " << endl;
  for (int i = 0; i < _nb_pts; i++)
  {
    cout.flush();
    cout << "Progression : " << (double)i/((double)(_nb_pts))*100 << "% \r";
    for (int j = 0 ; j < _nb_pts ; i++)
    {
      if(i==j)
        _matA.coeffRef(i,j)=1;
    }
  }

  SparseMatrix<double, ColMajor> A(_matA);
  BoundaryCondition_matA(A);

  cout << "Creation du solveur lineaire" << endl;
  _solverMethod.analyzePattern(A);
  _solverMethod.factorize(A);
}

void SystemBuilder::Build_sourceTerm()
{
  _sourceTerm.resize(_nb_pts);
  _sourceTerm.setZero();
  for (int i = 0; i < _nb_pts; i++)
  {
    _sourceTerm.coeffRef(i)=0;
  }
}

void SystemBuilder::BoundaryCondition_matA(SparseMatrix<double, RowMajor> A)
{
  cout << "-------------------------------------------------" << endl;
  cout << "Modification de A en fonction des conditions aux bords : " << endl;

  for (int i = 0; i < 2*_nb_pts; i++)
  {
    cout.flush();
    cout << "Progression : " << (double)i/(double)_nb_pts*100 << "% \r";
  }

  cout << "-------------------------------------------------" << endl;
}

void SystemBuilder::BoundaryCondition_sourceTerm(SparseVector<double>& B)
{
  cout << "-------------------------------------------------" << endl;
  cout << "Modification de b en fonction des conditions aux bords : " << endl;
  for (int i = 0; i < _nb_pts; i++)
  {
    if(i==0)
      _sourceTerm.coeffRef(i)=_sourceTerm.coeffRef(i)+1;
  }
}


#define _SYSTEMBUILDER_CPP
#endif

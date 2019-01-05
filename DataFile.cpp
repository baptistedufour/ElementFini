#ifndef _DATA_FILE_CPP

#include "DataFile.h"
#include <fstream>
#include <iostream>
#include <cmath>

using namespace std;

DataFile::DataFile(std::string file_name)
: _if_results(false), _if_gamma0_choice(false),
_file_name(file_name),  _if_Nmesh_choice(false), _if_sigma_choice(false),
_if_left_boundary_condition_choice(false),_if_right_boundary_condition_choice(false),
_if_Pk_choice(false),_if_solver_choice(false),_if_norm_L2_choice(false),_if_norm_H1_choice(false)
{}

  void DataFile::ReadDataFile()
  {
    ifstream data_file(_file_name.data());
    if (!data_file.is_open())
    {
      cout << "Unable to open file " << _file_name << endl;
      abort();
    }
    else
    {
      cout << "-------------------------------------------------" << endl;
      cout << "Reading data file " << _file_name << endl;
    }

    string file_line;

    while (!data_file.eof())
    {
      getline(data_file, file_line);
      if (file_line.find("Nmesh") != std::string::npos)
      {
        data_file >> _Nmesh_choice; _if_Nmesh_choice = true;
      }

      getline(data_file, file_line);
      if (file_line.find("gamma0") != std::string::npos)
      {
        data_file >> _gamma0; _if_gamma0_choice = true;
      }

      if (file_line.find("sigma") != std::string::npos)
      {
        data_file >> _sigma_choice; _if_sigma_choice = true;
        if (_sigma_choice == "constant")
        {
          data_file >> _a;
        }
        else if ((_sigma_choice == "line"))
        {
          data_file >> _a >> _b;
        }
        else if (_sigma_choice == "creneau")
        {
          data_file >> _a >> _b >> _c;
        }
        else
        {
          cout << "Only constant, line and creneau sigma are implemented." << endl;
          abort();
        }
      }

      if (file_line.find("source_fct") != std::string::npos)
      {
        data_file >> _source_fct_choice; _if_source_fct_choice = true;
        if (_source_fct_choice == "constant")
        {
          data_file >> _d;
        }
        else if (_source_fct_choice == "line")
        {
          data_file >> _d >> _e;
        }
        else if (_source_fct_choice == "creneau")
        {
          data_file >> _d >> _e >> _f;
        }
        else if (_source_fct_choice == "sinus")
        {
          //data_file >> _d >> _e >> _f;
        }
        else
        {
          cout << "Only constant, line, sinus and creneau source function are implemented." << endl;
          abort();
        }
      }

      if (file_line.find("left_boundary_condition") != std::string::npos)
      {
        data_file >> _left_boundary_condition_choice >> _ul;  _if_left_boundary_condition_choice=true;
        if ((_left_boundary_condition_choice != "neumann") && (_left_boundary_condition_choice != "dirichlet"))
        {
          cout << "Only Neumann and Dirichlet boundary condition are implemented." << endl;
          abort();
        }
      }

      if (file_line.find("right_boundary_condition") != std::string::npos)
      {
        data_file >> _right_boundary_condition_choice >> _ur; _if_right_boundary_condition_choice=true;
        if ((_right_boundary_condition_choice != "neumann") && (_right_boundary_condition_choice != "dirichlet"))
        {
          cout << "Only Neumann and Dirichlet boundary condition are implemented." << endl;
          abort();
        }
      }

      if (file_line.find("Pk_element") != std::string::npos)
      {
        data_file >> _Pk_choice; _if_Pk_choice = true;
        if ((_Pk_choice != 1) && (_Pk_choice != 2))
        {
          cout << "Only P1 and P2 numerical flows are implemented." << endl;
          abort();
        }
      }

      if (file_line.find("solver") != std::string::npos)
      {
        data_file >> _solver_choice; _if_solver_choice = true;
        if ((_solver_choice != "cholesky") && (_solver_choice != "LU"))
        {
          cout << "Only Cholesky and LU solver are implemented." << endl;
          abort();
        }
      }

      if (file_line.find("norm_L2") != std::string::npos)
      {
        data_file >> _norm_L2_choice; _if_norm_L2_choice = true;
      }

      if (file_line.find("norm_H1") != std::string::npos)
      {
        data_file >> _norm_H1_choice; _if_norm_H1_choice = true;
      }

      if (file_line.find("results") != std::string::npos)
      {
        data_file >> _results; _if_results = true;
      }
    }

    if (!_if_Nmesh_choice)
    {
      cout << "-------------------------------------------------" << endl;
      cout << "Beware - The default size Nmesh = 100 is used." << endl;
      _Nmesh_choice = 100;
    }
    if (!_if_sigma_choice)
    {
      cout << "-------------------------------------------------" << endl;
      cout << "Beware - The default sigma = 1 is used." << endl;
      _sigma_choice = "constant"; _a=1;
    }
    if (!_if_source_fct_choice)
    {
      cout << "-------------------------------------------------" << endl;
      cout << "Beware - The default source function 0 is used." << endl;
      _sigma_choice = "constant"; _d=0;
    }
    if (!_if_left_boundary_condition_choice)
    {
      cout << "-------------------------------------------------" << endl;
      cout << "Beware - The default left boundary condition is homogeous Dirichlet (0)." << endl;
      _left_boundary_condition_choice = "dirichlet"; _ul = 0;
    }
    if (!_if_right_boundary_condition_choice)
    {
      cout << "-------------------------------------------------" << endl;
      cout << "Beware - The default right boundary condition is homogeous Dirichlet (0)." << endl;
      _right_boundary_condition_choice = "dirichlet"; _ur = 0;
    }
    if (!_if_Pk_choice)
    {
      cout << "-------------------------------------------------" << endl;
      cout << "Beware - The default value (1) is used for the degree of Pk elements." << endl;
      _Pk_choice = 1;
    }
    if (!_if_gamma0_choice)
    {
      cout << "-------------------------------------------------" << endl;
      cout << "Beware - The default value (1000) is used for gamma0." << endl;
      _gamma0 = 1000;
    }
    if (!_if_norm_L2_choice)
    {
      cout << "-------------------------------------------------" << endl;
      cout << "Do not forget to say if you want the L2-norm of the solution in the data file." << endl;
      abort();
    }
    if (!_if_norm_H1_choice)
    {
      cout << "-------------------------------------------------" << endl;
      cout << "Do not forget to say if you want the H1-norm of the solution in the data file." << endl;
      abort();
    }
    if (!_if_results)
    {
      cout << "-------------------------------------------------" << endl;
      cout << "Beware - The default results folder name (Results) is used." << endl;
      _results = "Results/";
    }
  }

  #define _DATA_FILE_CPP
  #endif

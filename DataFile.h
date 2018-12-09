#ifndef _DATA_FILE_H

#include <string>
#include <vector>
#include <iostream>
// Définition de la classe

class DataFile {
private:
  std::string _file_name;
  double _Nmesh_choice,_Pk_choice;
  double _a, _b, _c;
  double _d, _e, _f;
  double _ul, _ur;

  std::string _sigma_choice;
  std::string _source_fct_choice;
  std::string _left_boundary_condition_choice;
  std::string _right_boundary_condition_choice;
  std::string _solver_choice;
  std::string _norm_L2_choice;
  std::string _results;
  std::string _gamma0_choice;

  bool _if_Nmesh_choice;
  bool _if_sigma_choice;
  bool _if_source_fct_choice;
  bool _if_left_boundary_condition_choice;
  bool _if_right_boundary_condition_choice;
  bool _if_Pk_choice;
  bool _if_solver_choice;;
  bool _if_norm_L2_choice;
  bool _if_results;
  bool _if_gamma0_choice;

public: // Méthodes et opérateurs de la classe
  DataFile(std::string file_name);
  void ReadDataFile();

  std::string Get_file_name() const {return _file_name;}
  std::string Get_sigma_choice() const {return _sigma_choice;};
  std::string Get_source_fct_choice() const {return _source_fct_choice;};
  std::string Get_right_BC_choice() const {return _right_boundary_condition_choice;};
  std::string Get_left_BC_choice() const {return _left_boundary_condition_choice;};
  std::string Get_solver_choice() const {return _solver_choice;};
  std::string Get_norm_L2_choice() const {return _norm_L2_choice;};
  std::string Get_results() const {return _results;};
  std::string Get_gamma0_choice() const {return _gamma0_choice;};

  double Get_Nmesh() const {return _Nmesh_choice;};
  double Get_Pk_choice() const { return _Pk_choice;}
  double Get_right_BC() const { return _ur;}
  double Get_left_BC() const { return _ul;}
  double Get_param_a() const { return _a;}
  double Get_param_b() const { return _b;}
  double Get_param_c() const { return _c;}
  double Get_param_d() const { return _d;}
  double Get_param_e() const { return _e;}
  double Get_param_f() const { return _f;}
};

#define _DATA_FILE_H
#endif

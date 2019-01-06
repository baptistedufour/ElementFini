#ifndef _SYSTEMBUILDER_H

#include "Sparse"
#include "Dense"
#include <SparseLU>
#include <SparseCholesky>
#include <fstream>
#include "DataFile.h"

#define PI 3.141592653

class SystemBuilder
{
  private:
    Eigen::SparseMatrix<double, Eigen::ColMajor> _matAm, _matAs, _matA;
    Eigen::SparseVector<double> _sourceTerm;
		Eigen::SparseVector<double> _sigma;


    double _nb_pts, _Asize;
    double  _dx, _Pk_choice;
    double _a, _b, _c;
    double _d, _e, _f;
    double _ur, _ul;
    double _gamma0;

    double _aL, _bL, _aR, _bR;
    std::string BC_left;
    std::string BC_right;
    std::string _src_choice;
    std::string _sigma_choice;

  public:
    SystemBuilder(DataFile* data_file);
    void Build_matA();
    void Build_sourceTerm();

    const Eigen::SparseMatrix<double, Eigen::ColMajor> & Get_matA() const {return _matA;};
    const Eigen::SparseVector<double> & Get_sourceTerm() const {return _sourceTerm;};
    const double & Get_nb_pts() const {return _nb_pts;};
    double Get_right_BC() const { return _ur;};
    double Get_left_BC() const { return _ul;};


    double Get_param_d() const { return _d;}
    double Get_param_e() const { return _e;}
    double Get_param_f() const { return _f;}
    double Get_sigma0();
    Eigen::SparseVector<double> Get_sigma() const { return _sigma;};

    std::string Get_sigma_choice() const {return _sigma_choice;};
    std::string Get_source_fct_choice() const {return _src_choice;};

};

#define _SYSTEMBUILDER_H
#endif

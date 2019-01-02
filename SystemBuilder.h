#ifndef _SYSTEMBUILDER_H

#include "Sparse"
#include "Dense"
#include <SparseLU>
#include <SparseCholesky>
#include <fstream>
#include "DataFile.h"

class SystemBuilder
{
  private:
    Eigen::SparseMatrix<double, Eigen::ColMajor> _matAm, _matAs, _matA;
    Eigen::SparseVector<double> _sourceTerm;
		std::vector<double> _sigma;


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

    Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>, Eigen::COLAMDOrdering<int> > _solverMethod;
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > _solverMethod2;

  public:
    SystemBuilder(DataFile* data_file);
    void Build_matA();
    void Build_sourceTerm();

    const Eigen::SparseMatrix<double, Eigen::ColMajor> & Get_matA() const {return _matA;};
    const Eigen::SparseVector<double> & Get_sourceTerm() const {return _sourceTerm;};
    const double & Get_nb_pts() const {return _nb_pts;};
    double Get_right_BC() const { return _ur;};
    double Get_left_BC() const { return _ul;};

    std::string Get_sigma_choice() const {return _sigma_choice;};
    std::string Get_source_fct_choice() const {return _src_choice;};

};

#define _SYSTEMBUILDER_H
#endif

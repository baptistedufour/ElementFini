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
    Eigen::SparseMatrix<double, Eigen::ColMajor> _matA;
    Eigen::SparseVector<double> _sourceTerm;
		std::vector<double> _sigma;


    double _Pk_choice;
    double _nb_pts;
    double _a, _b, _c;
    double _d, _e, _f;
    double _ur, _ul;

    Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>, Eigen::COLAMDOrdering<int> > _solverMethod;
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > _solverMethod2;

  public:
    SystemBuilder(DataFile* data_file);
    void Build_matA();
    void Build_sourceTerm();

    void BoundaryCondition_matA(Eigen::SparseMatrix<double, Eigen::RowMajor> A);
    void BoundaryCondition_sourceTerm(Eigen::SparseVector<double>& B);

    const Eigen::SparseMatrix<double, Eigen::ColMajor> & Get_matA() const {return _matA;};
    const Eigen::SparseVector<double> & Get_sourceTerm() const {return _sourceTerm;};
    const double & Get_nb_pts() const {return _nb_pts;};
};

#define _SYSTEMBUILDER_H
#endif

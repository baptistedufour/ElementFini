#ifndef _SYSTEMSOLVER_H

#include "Sparse"
#include "Dense"
#include <SparseLU>
#include <SparseCholesky>
#include <fstream>
#include "SystemBuilder.h"

class SystemSolver {
  private:
    SystemBuilder* _builder;
    std::string _source_fct_choice;
    std::string _sigma_choice;
    enum {constant,line,curve};

    Eigen::SparseVector<double>  _sol;
    Eigen::SparseVector<double>  _solEx;

    double _nb_pts, _ur, _ul, _dx;
    double _errorL, _errorH;

    Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>, Eigen::COLAMDOrdering<int> > _solverMethod;
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > _solverMethod2;

	public: // Méthodes et opérateurs de la classe
    SystemSolver(SystemBuilder* builder);

    void BuildSol();
    void BuildSolEx();
    void SaveSol();

    Eigen::SparseVector<double> & Get_Sol() {return _sol;};
    Eigen::SparseVector<double> & Get_ExactSol() {return _solEx;};

    void ErrorLinf();
    void ErrorL2();
    void ErrorH1();
};

#define _SYSTEMSOLVER_H
#endif

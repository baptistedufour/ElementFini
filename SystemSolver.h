#ifndef _SYSTEMSOLVER_H

#include "Sparse"
#include "Dense"
#include <SparseLU>
#include <SparseCholesky>
#include <fstream>
#include "SystemBuilder.h"

#define PI 3.141592653

class SystemSolver {
  private:
    SystemBuilder* _builder;
    std::string _source_fct_choice;
    std::string _sigma_choice;
    enum {constant,line,curve};

    Eigen::SparseVector<double> _sigma;
    Eigen::SparseVector<double>  _sol;
    Eigen::SparseVector<double>  _solEx;

    double _nb_pts, _ur, _ul, _dx;
    double _errorL, _errorH;

    Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>, Eigen::COLAMDOrdering<int> > _solverMethod;
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > _solverMethod2;
    Eigen::ConjugateGradient< Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper > _solverMethod3;


	public: // Méthodes et opérateurs de la classe
    SystemSolver(SystemBuilder* builder);

    void BuildSol();
    void BuildSolEx();
    void SaveSol();

    Eigen::SparseVector<double> & Get_Sol() {return _sol;};
    Eigen::SparseVector<double> & Get_ExactSol() {return _solEx;};

    void ErrorLinf();
    void ErrorL2_poly();
    void ErrorL2_sin();
    void ErrorH1_poly();
    void ErrorH1_sin();
};

#define _SYSTEMSOLVER_H
#endif

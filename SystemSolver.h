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
    Eigen::VectorXd  _sol2;

    double _nb_pts;

    Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>, Eigen::COLAMDOrdering<int> > _solverMethod;
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double> > _solverMethod2;

	public: // Méthodes et opérateurs de la classe
    SystemSolver(SystemBuilder* builder);

    void BuildSol();
    void SaveSol();

    Eigen::VectorXd & Get_Sol() {return _sol2;};

    const Eigen::VectorXd ExactSolution() const;
    double Error(Eigen::SparseVector<double> sol);


};

#define _SYSTEMSOLVER_H
#endif

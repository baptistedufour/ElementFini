#include <iostream>
#include <fstream>
#include <chrono>
#include "SystemSolver.h"

using namespace std;
using namespace Eigen;

int main(int argc, char** argv)
{

  if (argc < 2)
  {
    cout << "Please, enter the name of your data file." << endl;
    abort();
  }
  const string data_file_name = argv[1];

  // --------------------------- Fichier de données ----------------------------
  DataFile* data_file = new DataFile(data_file_name);
  data_file->ReadDataFile();
  // ---------------------------------------------------------------------------

  cout << "-------------------------------------------------" << endl;
  cout << "Search u such that : " << endl;
  cout << "dx(sigma*dx(u)) = f,  B" << endl;
  cout << "-------------------------------------------------" << endl;

  // -------------------- Construction du système linéaire  --------------------
  auto start = chrono::high_resolution_clock::now();
  SystemBuilder* builder = new SystemBuilder(data_file);
  builder->Build_matA();
  builder->Build_sourceTerm();
  auto finish = chrono::high_resolution_clock::now();
  double t = chrono::duration_cast<chrono::milliseconds>(finish-start).count();
  // ---------------------------------------------------------------------------

  // --------------------- Resolution du système linéaire  ---------------------
  SystemSolver* solver = NULL;
  solver = new SystemSolver(builder);

  start = chrono::high_resolution_clock::now();
  solver->BuildSol();
  solver->SaveSol();
  finish = chrono::high_resolution_clock::now();
  t = chrono::duration_cast<chrono::milliseconds>(finish-start).count();

  cout << "Time to compute : "<< t*0.001 << " seconds" << endl;
  // ---------------------------------------------------------------------------

  //------------------- Si on connait la solution exacte -----------------------

  if (((data_file->Get_source_fct_choice() == "constant")||(data_file->Get_source_fct_choice() == "line") || (data_file->Get_source_fct_choice() == "sinus"))
    &&(data_file->Get_sigma_choice() == "constant"))
  {
    cout << "Error : " << endl;
    solver->ErrorLinf();
  if(((data_file->Get_source_fct_choice() == "constant")&&(data_file->Get_norm_L2_choice() == "yes"))||((data_file->Get_source_fct_choice() == "line")&&(data_file->Get_norm_L2_choice() == "yes")))
      solver->ErrorL2();
  else if ((data_file->Get_norm_L2_choice() == "yes")||(data_file->Get_source_fct_choice() == "sinus"))
      solver->ErrorL2_sin();
    if((data_file->Get_source_fct_choice() == "constant")&&(data_file->Get_norm_H1_choice() == "yes"))
      solver->ErrorH1();
  }
  // ---------------------------------------------------------------------------

  cout << "-------------------------------------------------" << endl;
  delete solver;
  delete builder;
  delete data_file;

  return 0;
}

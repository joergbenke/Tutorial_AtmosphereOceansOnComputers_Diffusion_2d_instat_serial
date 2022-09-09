#ifndef DIFFUSION_FDM_1D_H
#define DIFFUSION_FDM_1D_H

#include <vector>
#include <string>

class Diffusion_FDM_1d
{

 public:

  // (Standard-)Constructor
  // Old: AdvectionDiffusionReaction_FDM_Jacobi( double epsilon, unsigned long int n, unsigned long int iter ); //: error(epsilon), n(n), n_iter(iter){};
  Diffusion_FDM_1d( ); 
  
  // Create and init fields and boundary conditions
  double* createAndInitFields( double *x, double value );

  // Do the iteration in the diffusin case
  unsigned long int diffusion_ftcs( double *x_new, double *x_old, double *rhs );

  // Write data to disc
  int writeData( double *x );

  // Read namelist
  int readData( );

  // Set parameters from namelist
  void setParameters( std::vector<std::string> lines );

  // Caclculations of FLOPs
  double calculateFlops( unsigned long int );

  
private:

  unsigned long int n;       // Number of mesh nodes in one direction
  unsigned long int n_iter;  // Number of iterations of the iteration process
  double error;              // error value, when the calculation stopps

  double Axx; // Diffusioncoefficient xx direction
  double Ayy; // Diffusioncoefficient yy direction
  double Azz; // Diffusioncoefficient zz direction
  double c; // Reactioncoefficient
  
};

#endif

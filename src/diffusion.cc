#include <iostream>
#include <fstream>
#include <cmath>
#include <new>
#include <vector>

#include "Diffusion_FDM_2d.h"


//#define DEBUG false

using namespace std;

double calc_elapsed_time( clock_t start, clock_t end )
{

  return ( ( double ) end - start ) / CLOCKS_PER_SEC;

}



int main( int argc, char **argv )
{

  clock_t start, end;
  double elapsed_time{0.0};

  double *x_new{nullptr}, *x_old{nullptr}, *rhs{nullptr};

  //
  // Instantiate and initialize new object

  Diffusion_FDM_2d adr{}; 

  //
  // Read namelist

  adr.readData( ); 
  
  
  //
  // Create and init fields

  start = clock();

  x_new = adr.createAndInitFields( x_new, 0.0 ); // Create and init new mesh with 0 
  x_old = adr.createAndInitFields( x_old, 0.0 ); // Create and init old mesh with 0
  rhs   = adr.createAndInitFields( rhs, 0.90 );  // Create and init right hand side

  end = clock();
  elapsed_time = calc_elapsed_time( start, end ); // Calculate elapsed time for creating ind initializing meshes
  cout << "Elapsed time initialization: " << elapsed_time << endl;


  //
  // Jacobian calculations
  
  start = clock();
  
  unsigned long int n_iterations = adr.diffusion_ftcs( x_new, x_old, rhs ); // Calculate the values in the mesh points with jacobi 
  double flops = adr.calculateFlops( ); // Calculate the number of FLOPs of jacobi discretization

  end = clock();
  elapsed_time = calc_elapsed_time( start, end ); // Calculate elapsed time for jacobi discrtization 
  cout << "Elapsed time jacobi: " << elapsed_time << endl;
  flops = flops / elapsed_time; // Calculate the number of FLOPS of jacobi discretization
  cout << "Number of FLOPS: " << flops << endl; 
  
  //
  // Output of the results

  start = clock();

  adr.writeData( x_new );

  end = clock();
  elapsed_time = calc_elapsed_time( start, end );
  cout << "Elapsed time IO: " << elapsed_time << endl;

  //
  // Delete arrays

  delete x_new;
  delete x_old;
  delete rhs;

  
  return EXIT_SUCCESS; 

}

// #include <unistd.h> 
// usleep(10000000);
// clock_t: Arithmetic (until C11)Real (since C11) type capable of representing the processor time used by a process. It has implementation-defined range and precision. 

//using std::cout;
//using std::endl;

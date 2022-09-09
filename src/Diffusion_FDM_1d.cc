#include <cmath>
#include <iostream>
#include <fstream>
#include <new>
#include <string>
#include <vector>

#include "Diffusion_FDM_1d.h"

using namespace std;

// Standardconstructor
Diffusion_FDM_1d::Diffusion_FDM_1d( )  
{

  this->error = 10.0;  
  this->n = 10;
  this->n_iter = 100;
  
}
  
// Create and init fields and boundary conditions
double* Diffusion_FDM_1d::createAndInitFields( double *x, double value )
{

  //
  // Create mesh
  
  try
    {
	
      x = new double[ n * n ];
	
    }
  catch( bad_alloc &badAlloc )
    {
	
      cerr << "bad_alloc caught, not enough memory: " << badAlloc.what() << endl;
      exit( 1 );
	
    }
    
  //
  // Create and init problem
  for( unsigned long int i = 0; i < n * n; ++i )
    x[ i ] = value;
    
  //
  // Initialize boundary conditions

  // Upper boundary
  for( unsigned long int i = 0; i < n; ++i )
    x[ i ] = 1.0;
 
  // Lower boundary
  for( unsigned long int i = 0; i < n; ++i )
    x[ ( n - 1 ) * n + i ] = 1.0;

  // Left boundary
  for( unsigned long int i = 0; i < n; ++i )
    x[ n * i ] = 1.0;

  // Right Boundary
  for( unsigned long int i = 0; i < n; ++i )
    x[ n * (i + 1) - 1 ] = 1.0;

  return x;
    
}
  

//
// Do the iteration in the diffusin case
unsigned long int Diffusion_FDM_1d::diffusion_ftcs( double *x_new, double *x_old, double *rhs )
{

  unsigned long int iter{0};
  double err_it{10.0};
  register double sum{0.0};
  double h = static_cast<double>( 1 ) / ( n - 1 ); 
  double h_squared = pow( h, 2 );
  double scaling = static_cast<double> (1) / ( 2 * ( Axx + Ayy ) - c * h_squared ); 
  
  //
  // Jacobi iteration
  
  while( ( err_it > error ) && ( iter < n_iter ) )
    {

#if DEBUG
      cout << "Iteration: " << iter << " " << "Error: " << err_it << endl ; 
#endif
      
      //
      // Calculation

      for( unsigned long int i = 1; i < n - 1; ++i )  
	for( unsigned long int j = 1; j < n - 1; ++j )  
	  {

	    //sum = 0.25 * ( x_old[ i*n+(j+1) ] + x_old[ i*n+(j-1) ] + x_old[ (i+1)*n+j ] + x_old[ (i-1)*n+j ] + h_squared*rhs[ i*n+j ]); 
	    x_new[ i*n+j ] = scaling * ( Axx*x_old[ i*n+(j+1) ] + Axx*x_old[ i*n+(j-1) ] + Ayy*x_old[ (i+1)*n+j ] + Ayy*x_old[ (i-1)*n+j ] + h_squared*rhs[ i*n+j ]); 

	  }

      sum = 0.0;
      for( unsigned long int i = 1; i < n-1; ++i )  
	for( unsigned long int j = 1; j < n-1; ++j )  
	  sum += ( x_new[ i*n+j ] - x_old[ i*n+j ] ) * ( x_new[ i*n+j ] - x_old[ i*n+j ] );
	  
      err_it = sqrt( sum );
      sum = 0.0;
      ++iter;

      for( unsigned long int i = 1; i < n-1; ++i )  
	for( unsigned long int j = 1; j < n-1; ++j )  
	  x_old[ i*n+j ] = x_new[ i*n+j ]; 

    }

#if DEBUG
  cout << "Iteration: " << iter << " " << "Error: " << err_it << endl ; 
#endif

  return iter;

}


// Write data to disk
int Diffusion_FDM_1d::writeData( double *x )
{

  ofstream jacobi_output_file( "jacobi_output_1.txt" );

  //
  // Write to the file

  if ( jacobi_output_file.is_open() )
    {
      
      for( unsigned long int i = 0; i < n; ++i )
	{
	  
	  for( unsigned long int j = 0; j < n; ++j )
	    jacobi_output_file << x[ i * n + j ] << " "; 
	  
	  jacobi_output_file << endl;
	  
	}

    } 

  //
  // Close the file

  jacobi_output_file.close();
  
  return EXIT_SUCCESS;
  
}

// Read the paramters for the numerical discretization from the namelist file
int Diffusion_FDM_1d::readData( )
{

  size_t pos = 0;

  string filename{"jacobi_namelist.txt"};
  string line{};
  string del=":";
  string token;
  vector<string> lines;
    
  ifstream jacobi_namelist_file( filename );
  
  // Write to the file
  if ( !jacobi_namelist_file.is_open() )
    {
      
      cerr << "Could not open the file - '" << filename << "'" << endl;
      return EXIT_FAILURE;

    }

  if ( jacobi_namelist_file.is_open() )
    {
      
      while( getline( jacobi_namelist_file, line ) )
	{
	  
	  while( (pos = line.find( del ) ) != string::npos )
	    {
	      
	      token = line.substr( 0, pos );
	      line.erase( 0, pos + del.length() );

	    }
	  
	  lines.push_back( line ); // Add new vector element

	}
      
      for (const auto &i : lines )
        cout << i << endl;

      // Assign new values to numerical model
      setParameters( lines );
      
      // Close the file
      jacobi_namelist_file.close();
  
    }

  return EXIT_SUCCESS;

}

// Set parameters 
void Diffusion_FDM_1d::setParameters( std::vector<std::string> lines )
{

  n = stoi( lines.at( 0 ) );
  n_iter = stoi( lines.at( 1 ) );
  error = stod( lines.at( 2 ) );
  
  Axx = stod( lines.at( 3 ) );
  Ayy = stod( lines.at( 4 ) );
  Azz = stod( lines.at( 5 ) );
  c = stod( lines.at( 6 ) );
  
  cout << "n = " << n << "\n n_iter = " << n_iter << "\n error = " << error << endl;
  cout << "Axx = " << Axx << "\n Ayy = " << Ayy << "\n Azz = " << Azz << "\nc = " << c << endl;

  return;
  
}

// Calculate FLOPs
double Diffusion_FDM_1d::calculateFlops( unsigned long int n_iter )
{
  
  return  n_iter * ( ( 6+4 ) * (n-2) * (n-2) + 1 );
  
}

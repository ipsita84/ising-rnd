// 2nd Renyi entropy for classical 2d Ising model in zero magnetic field 
//Metropolis algorithm employed
//Parameters that can be changed for different runs:
//J, axis1, axis2, N_mc

#include <iostream>
#include <fstream>
#include <ctime>
#include <boost/multi_array.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include <gsl/gsl_integration.h>

// gen is a variable name
// Its data-type is boost::random::mt19937
boost::random::mt19937 gen;
using namespace std;

typedef
 boost::multi_array < int, 2 > array_2d;
// typedef keyword allows you to create an alias fo a data type

// Define global scopes to use them across all functions
double J = 1.0;
const unsigned int axis1 = 10, axis2 = 10;
//axis1 should be of even no. of sites
// above assigns length along each dimension of the 2d configuration
//No.of Monte Carlo updates we want
unsigned int N_mc = 100;

//Function templates
double en_avg(double beta);
double modi_en(double beta) ;
int roll_coin(int a, int b);
double random_real(int a, int b);
double energy_tot(array_2d sitespin);
double nn_energy(array_2d sitespin, unsigned int row, unsigned int col);
double f (double x, void * params);
double g (double x, void * params);

int main()
{	

	double kT(0), T_min(0), T_max(0), del_T(0);

	cout << "Enter minimum T (in units of k)" << endl;
	cin >> T_min;

	cout << "Enter maximum T (in units of k)" << endl;
	cin >> T_max;

	cout << "Enter increment of T (in units of k) at each step" << endl;
	cin >> del_T;

	
	double mut_info(0); //mutual information I_2

	ofstream fout("mutual-info.dat"); // Opens a file for output
	
	gsl_integration_cquad_workspace *w
	 = gsl_integration_cquad_workspace_alloc (1000);

	for (kT = T_min; kT < T_max + del_T; kT += del_T)
	{	double beta = 1.0/kT ;
	 	mut_info = 0 ;
	 	
	 	double term1(0), term2(0), term3(0), abs_error(0);
	 	gsl_function F;
  		F.function = &f;
		F.params = NULL;
//Function: int gsl_integration_qags (const gsl_function * f,double a,double b, double epsabs, double epsrel,size_t limit,gsl_integration_workspace * workspace,double * result, double *abserr)  int
//gsl_integration_cquad (const gsl_function * f, double a, double b, double epsabs, double epsrel, gsl_integration_cquad_workspace * workspace, double * result, double * abserr, size_t * nevals)		
  		gsl_integration_cquad (&F, 0, beta, 1e-6, 1e-4, w, &term2, &abs_error, NULL);
  		gsl_integration_cquad (&F,beta, 2.0*beta, 1e-6, 1e-4, w, &term3, &abs_error, NULL);
  		
  		F.function = &g;
  		gsl_integration_cquad (&F, 0, beta, 1e-6, 1e-4, w, &term1, &abs_error, NULL);
		
	        
	        mut_info =2.0*term1 +3.0*term2 + term3;
	  
		fout << kT / J << '\t' << mut_info << endl;
	
	}

	

	fout.close();
	return 0;

}

//performing numerical integration using gsl
double f (double beta, void * params) 
{
  double f = -en_avg(beta) ;
  return f;
}

double g (double beta, void * params) 
{
  double g = modi_en(beta) ;
  return g;
}

//function to calculate avg energy for replica spin config at temp 1/beta
//logic: for a[n1][n2], a[n1] is n1 copies of 1d array of length n2

double modi_en(double beta)
{
	unsigned int sys_size = axis1 * axis2;
	unsigned int row, col, label;
	double r(0), acc_ratio(0) ;
	
	//define replica 1 spin configuration array
	array_2d sitespin1(boost::extents[axis1][axis2]);
	
	//define replica 2 spin configuration array
	array_2d sitespin2(boost::extents[axis1][axis2]);
	
	//For subsystem A,both replicas have same spin configuration
	for (unsigned int i = 0; i < axis1/2; ++i)
		for (unsigned int j = 0; j < axis2; ++j)
			{	sitespin1[i][j] = 2 * roll_coin(0, 1) - 1;
				sitespin2[i][j] = sitespin1[i][j];
			}
			
	//For subsystem B, the two replicas have independent spin configurations
	for (unsigned int i = axis1/2 ; i < axis1; ++i)
		for (unsigned int j = 0; j < axis2; ++j)
			{	sitespin1[i][j] = 2 * roll_coin(0, 1) - 1;
				sitespin2[i][j] = 2 * roll_coin(0, 1) - 1;
			}
			

	double energy = energy_tot(sitespin1);
	energy += energy_tot(sitespin2);
		

	double en_sum(0);

	for (unsigned int i = 1; i <= N_mc; ++i)
	{
		for (unsigned int j = 1; j <= 3*sys_size/2; ++j)
		{	//Choose a random spin site for the entire 2 replica system
			double energy_diff(0);
			label = roll_coin(1,2*sys_size);
			
			//if the random spin site is located in layer 1
			if (label <= sys_size)
			{	if (label % axis2 == 0)
				{	row = (label / axis2) - 1;
					col = axis2 -1 ; 
				}
				else
				{	col = label % axis2 - 1;
					row = (label-col-1)/axis2;
				} 
			
				energy_diff=-2.0*nn_energy(sitespin1, row, col);
				if (row < axis1/2)
					energy_diff +=-2.0*nn_energy(sitespin2, row, col);

				//Generate a random no. r such that 0 < r < 1

				r = random_real(0, 1);
				acc_ratio = exp(-1.0 * energy_diff *beta);

				//Spin flipped if r <= acceptance ratio
				if (r <= acc_ratio)
				{
					sitespin1[row][col] *= -1;
					if (row < axis1/2) 
						sitespin2[row][col] *=-1;
							//this line on addition creates trouble
					energy += energy_diff;
				}
			}
			
			//if the random spin site is located in layer 2
			if (label > sys_size)
			{	label -= sys_size;
				if (label % axis2 == 0)
				{	row = (label / axis2) - 1;
					col = axis2 -1 ; 
				}
				else
				{	col = label % axis2 - 1;
					row = (label-col-1)/axis2;
				} 
			
			
				energy_diff=-2.0*nn_energy(sitespin2, row, col);
				if (row < axis1/2)
					energy_diff +=-2.0*nn_energy(sitespin1, row, col);

					 
				//Generate a random no. r such that 0 < r < 1

				r = random_real(0, 1);
				acc_ratio = exp(-1.0 * energy_diff *beta);

				//Spin flipped if r <= acceptance ratio
				if (r <= acc_ratio)
				{
					sitespin2[row][col] *= -1;
					if (row < axis1/2) 
						sitespin1[row][col] *=-1;
						
					energy += energy_diff;
				}
			}

		}

		en_sum += energy;
	}
	double avg_en =en_sum / N_mc ;
	return avg_en ;
}



//function to calculate avg energy for a spin config at temp 1/beta
double en_avg(double beta)
{
	unsigned int sys_size = axis1 * axis2;
	array_2d sitespin(boost::extents[axis1][axis2]);

	// stores the spin configuration of the system
	//initial state chosen by random no. generator above
	for (unsigned int i = 0; i < axis1; ++i)
		for (unsigned int j = 0; j < axis2; ++j)
			sitespin[i][j] = 2 * roll_coin(0, 1) - 1;

	double energy = energy_tot(sitespin);
		

	double en_sum(0);

	// spins labelled by label = row*axis2 + col + 1
	//where 0 <= row < axis1 and 0 <= col < axis2

	for (unsigned int i = 1; i <= N_mc; ++i)
	{
		for (unsigned int j = 1; j <= sys_size; ++j)
		{	//Now choose a random spin site with site no.=label

			unsigned int label, row, col ;
			label = roll_coin(1, sys_size);
			if (label % axis2 == 0)
			{	row = (label / axis2) - 1;
				col = axis2 -1 ; 
			}
			else
			{	col = label % axis2 - 1;
				row = (label-col-1)/axis2;
			}

			double energy_diff =-2 * nn_energy(sitespin, row, col);

			//Generate a random no. r such that 0 < r < 1

			double r = random_real(0, 1);
			double acc_ratio = exp(-1.0 * energy_diff *beta);

			//Spin flipped if r <= acceptance ratio
			if (r <= acc_ratio)
			{
				sitespin[row][col] *= -1;
				energy += energy_diff;
			}
		}

		en_sum += energy;
	}
	double avg_en = en_sum / N_mc ;
	return avg_en ;
}

//function to generate random integer
// between 2 integers a & b, including a & b
int roll_coin(int a, int b)
{
	boost::random::uniform_int_distribution <> dist(a, b);

	return dist(gen);

}

//function to generate random real no.
// between 2 integers a & b, including a & excluding b

double random_real(int a, int b)
{
	boost::random::uniform_real_distribution <> dist(a, b);
	// uniform_real_distribution: continuous uniform distribution 
	//on some range [min, max) of real number
	return dist(gen);

}

//function to calculate total energy
//for a given spin configuration
//with periodic boundary conditions

double energy_tot(array_2d sitespin)
{
	double energy = 0;

	for (unsigned int i = 0; i < axis1 - 1; ++i)
	{
		for (unsigned int j = 0; j < axis2 - 1; ++j)
		{
			energy -= J * sitespin[i][j] * sitespin[i + 1][j];

			energy -= J * sitespin[i][j] * sitespin[i][j + 1];

		}
	}

	//periodic boundary conditions
	for (unsigned int j = 0; j < axis2; ++j)
		energy -= J * sitespin[axis1-1][j] * sitespin[0][j];

	for (unsigned int i = 0; i < axis1; ++i)
		energy -= J * sitespin[i][axis2-1] * sitespin[i][0];

	return energy;
}


//Calculating interaction energy change for spin 
//at random site->(row,col) with its nearest neighbours
double nn_energy(array_2d sitespin, unsigned int row, unsigned int col)
{
	double nn_en = 0;

	if ( row > 0 && row < axis1 - 1)
	{
		nn_en -= J * sitespin[row][col] * sitespin[row-1][col];
		nn_en -= J * sitespin[row][col] * sitespin[row+1][col];
	}

	if (col > 0 && col < axis2-1)
	{
		nn_en -= J * sitespin[row][col] * sitespin[row][col-1];
		nn_en -= J * sitespin[row][col] * sitespin[row][col+1];
	}

	if (row == 0)
	{
		nn_en -= J * sitespin[0][col] * sitespin[axis1-1][col];
		nn_en -= J * sitespin[0][col] * sitespin[1][col];

	}

	if (row == axis1-1)
	{
		nn_en -= J * sitespin[axis1-1][col] * sitespin[axis1-2][col];
		nn_en -= J * sitespin[axis1-1][col] * sitespin[0][col];

	}

	if (col == 0)
	{
		nn_en -= J * sitespin[row][0] * sitespin[row][axis2-1];
		nn_en -= J * sitespin[row][0] * sitespin[row][1];

	}

	if (col == axis2-1)
	{
		nn_en -= J * sitespin[row][axis2-1] * sitespin[row][axis2-2];
		nn_en -= J * sitespin[row][axis2-1] * sitespin[row][0];

	}
	return nn_en;
}

//  g++ -Wall -O3 -lgsl -lgslcblas -lm ising-2d-Renyi2.cc -o testo

// Considering modified 2d ising model using a 2d array in zero magnetic field
//Metropolis algorithm employed
// Modification: J has random sign at each site but same magnitude

//fix the Js at the start of the simulation
//keep the simulation running always with the same J
//-> called one "realization of disorder"

#include <iostream>
#include <fstream>
#include <ctime>
#include <vector>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>

// gen is a variable name
// Its data-type is boost::random::mt19937
boost::random::mt19937 gen;
using namespace std;

// Define global scopes 'J' and 'size' to use it across all functions
double J = 0;
const unsigned int size = 10;

//function to generate random integer
// between 2 integers a & b, including a & b
int roll_coin(int a, int b)
{
	boost::random::uniform_int_distribution <> dist(a, b);

	return dist(gen);

}

//function to calculate total energy
//for a given spin configuration
//with periodic boundary conditions

double energy_tot(int sitespin[][size], int J_x[][size], int J_y[][size])
{
	double energy = 0;

	for (unsigned int i = 0; i < size-1; ++i)
	{
		for (unsigned int j = 0; j < size-1; ++j)
		{
			energy -= J*J_x[i][j]*sitespin[i][j]*sitespin[i+1][j];
			    
			energy -= J*J_y[i][j]*sitespin[i][j]*sitespin[i][j+1];
			    
		}
	}
	
	//periodic boundary conditions
	for (unsigned int i=0 ; i < size ; ++i)
	{
		energy -= J*J_x[size-1][i]*sitespin[size-1][i] * sitespin[0][i];
		energy -= J*J_y[i][size-1]*sitespin[i][size-1] * sitespin[i][0];
	}
			
	return energy ;
}


//function to calculate magnetization
//for a given spin configuration
	
double mag_tot(int sitespin[][size])
{
	int mag = 0;

	for (unsigned int i = 0; i < size ; i++)
	{
		for (unsigned int j = 0; j < size ; j++) mag +=sitespin[i][j];
	
	}
			
	return mag ;
}


//function to generate random real no.
// between 2 integers a & b, including a & excluding b

double random_real(int a, int b)
{
	boost::random::uniform_real_distribution <> dist(a, b);
	// uniform_real_distribution: continuous uniform distribution 
	//on some range [min, max) of real numbers
	

	return dist(gen);

}


int main()
{
	int sitespin[size][size];
	// stores the spin configuration of the system
	//initial state chosen by random no. generator above
	
	
	for (unsigned int i = 0; i < size; ++i)
		for (unsigned int j = 0; j < size; ++j)
			//sitespin[i][j] = pow(-1,roll_coin(0,1));
			sitespin[i][j] = 2*roll_coin(0,1) - 1;

	double kT(0), beta(0);

	cout << "Enter T (in units of k)" << endl;
	cin >> kT;
	
	beta = 1.0/kT;
	
	
	
	cout << "Enter magnitude of coupling constant J > 0 " << endl;
	cin >> J;
	
	//Assign random sign to each NN bond
	//store in an array
	
	int J_x[size][size] , J_y[size][size] ;
	for (unsigned int i = 0; i < size; ++i)
	{
		for (unsigned int j = 0; j < size; ++j)
		{
			J_x[i][j] = 2*roll_coin(0,1) - 1; //pow(-1,roll_coin(0,1));
			J_y[i][j] = 2*roll_coin(0,1) - 1; //pow(-1,roll_coin(0,1));
		}
	}
	
	
	
	
	
	double energy = energy_tot(sitespin, J_x, J_y);
	
	//Input from terminal number of equilibration steps you wnat

	unsigned int maxstep(0);

	cout << "Enter no. of equilibration steps" << endl;
	cin >> maxstep;
	
	double en_sum(0), mag_sum(0);
	
	
	ofstream fout("2d_rnd_J.dat"); // Opens a file for output

	for (unsigned int i = 1; i <= maxstep; i++)
	{

		//Now choose a random spin site, say, i

		unsigned int rnd_row , rnd_col;
		rnd_row = roll_coin(1, size)-1;
		rnd_col = roll_coin(1, size)-1;

		//Calculating new energy for flipping spin at 
		//random site->(rnd_row,rnd_col)
		sitespin[rnd_row] [rnd_col] *= -1;

		double new_energy = energy_tot(sitespin, J_x, J_y);


		//Calculating change in energy for the above spin flip

		double energy_diff = new_energy - energy;
		
		

		//Generate a random no. r such that 0 < r < 1

		double r = random_real(0, 1);
		double acc_ratio = exp(-1.0 * energy_diff * beta);

		if (r > acc_ratio) sitespin[rnd_row] [rnd_col] *= -1;
		//Spin not flipped if r > acceptance ratio
			
		
		//Given the energy of ising system at a selection of times 
		// during the simulation, we can average them to find 
		//the estimators of internal energy
		//dividing it by no. of sites gives internal energy per site
		
		en_sum += energy_tot(sitespin, J_x, J_y ) / (size*size) ;
		mag_sum += mag_tot (sitespin ) * 1.0 / (size*size) ;
		
		if ((i % 1000) == 0)
			fout << i * 1.0 / (size*size) << '\t' << en_sum / i 
			<< '\t' << mag_sum / i << endl;

	}
	
	
	fout.close();
	
	
	return 0;
	
	}
	
	
	
	
//  g++ -Wall rnd-sign-ising-2d.cc -o testo
// ./testo


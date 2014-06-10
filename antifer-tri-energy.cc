//g++ -Wall -O3 antifer-tri-energy.cc -o testo

// Avg energy per spin vs beta for 2d antiferromagnetic triangular ising model
//in zero magnetic field 
//--http://www.boost.org/doc/libs/1_39_0/libs/multi_array/doc/user.html


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
// above assigns length along each dimension of the 2d configuration
unsigned int N_mc = 10000; //No.of Monte Carlo updates we want


//Function templates
int roll_coin(int a, int b);
double random_real(int a, int b);
double energy_tot(array_2d sitespin);
double nn_energy(array_2d sitespin, unsigned int row, unsigned int col);
double avg_en (double kT);
double f (double x, void * params);

int main()
{	

	double beta_min(0), beta_max(0), del_beta(0);
	

	cout << "Enter minimum beta" << endl;
	cin >> beta_min;

	cout << "Enter maximum beta" << endl;
	cin >> beta_max;

	cout << "Enter increment of beta at each step" << endl;
	cin >> del_beta;
	

	ofstream fout("triang_en_vs_beta.dat");	// Opens a file for output
	

	for (double beta = beta_min; beta < beta_max + del_beta; beta += del_beta)
	{	
		double en_per_spin = avg_en(beta)/(axis1 * axis2);
	  
		fout << beta << '\t' << en_per_spin << endl;
	
	}
	

	fout.close();
	
	return 0;
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
	double energy(0);
	int u1(0), v1(0), v2(0), sum(0);
	

	for (unsigned int row = 0; row < axis1 ; ++row)
	{	u1 = row+1;
		if (row == axis1-1) u1=0;
		for (unsigned int col = 0; col < axis2 ; ++col)
		{	v1 = col+1 ; v2 = col-1 ;
			if (col == axis2-1) v1 = 0;
			if (col == 0) v2 = axis2-1;
			sum = sitespin[row][v1]+sitespin[u1][col]+sitespin[u1][v2];
			energy +=0.25*J * sitespin[row][col] * sum ;
                        // 0.25 factor is due to (1/2)*(1/2) from s_i s_j
		}
	}

	//periodic boundary conditions employed
	
	return energy;
}



//Calculating interaction energy change for spin 
//at random site->(row,col) with its eight nearest neighbours
double nn_energy(array_2d sitespin, unsigned int row, unsigned int col)
{
	double nn_en(0);
	int u1(0), u2(0), v1(0), v2(0), sum(0);
	u1 = row + 1;
	u2 = row - 1;
	v1 = col + 1;
	v2 = col - 1;
	if (row == axis1-1) u1 = 0;
	if (row == 0) u2 = axis1-1;
	if (col == axis2-1) v1 = 0;
	if (col == 0) v2 = axis2-1;
	sum =sitespin[u1][col] + sitespin[u2][col] + sitespin[u1][v2]; 
	sum += sitespin[u2][v1] + sitespin[row][v1] + sitespin[row][v2];
	nn_en = 0.25*J * sitespin[row][col] * sum ;
	return nn_en;
}

//calculating avg total energy at inverse temp beta
double avg_en (double beta)
{	// Create a 2d array that is axis1 * axis2
		array_2d sitespin(boost::extents[axis1][axis2]);

		//stores the spin configuration of the system
		//initial state chosen by random no. generator above
		for (unsigned int i = 0; i < axis1; ++i)
			for (unsigned int j = 0; j < axis2; ++j)
				sitespin[i][j] = 2 * roll_coin(0, 1) - 1;

		double energy = energy_tot(sitespin);

		double en_sum(0);

		unsigned int sys_size = axis1 * axis2;

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
				double acc_ratio = exp(-1.0*energy_diff*beta);

				//Spin flipped if r <= acceptance ratio
				if (r <= acc_ratio)
				{
					sitespin[row][col] *= -1;
					energy += energy_diff;
				}
			}
			
			en_sum += energy;

		}

		return en_sum / N_mc;
}




// Considering 2d ising model using ./boost multi array 
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
const unsigned int axis1 = 15, axis2 = 15;
// above assigns length along each dimension of the 2d configuration

//Function templates
int roll_coin(int a, int b);
double random_real(int a, int b);
double energy_tot(array_2d sitespin);
double mag_tot(array_2d sitespin);
double nn_energy(array_2d sitespin, unsigned int rnd_row, unsigned int rnd_col);

int main()
{	//No.of Monte Carlo updates we want
	unsigned int N_mc = 10000;

	double kT(0), T_min(0), T_max(0), del_T(0);

	cout << "Enter minimum T (in units of k)" << endl;
	cin >> T_min;

	cout << "Enter maximum T (in units of k)" << endl;
	cin >> T_max;

	cout << "Enter increment of T (in units of k) at each step" << endl;
	cin >> del_T;

	ofstream fout("2d_T_dep.dat");	// Opens a file for output

	for (kT = T_min; kT < T_max + del_T; kT += del_T)

	{	// Create a 2d array that is axis1 * axis2
		array_2d sitespin(boost::extents[axis1][axis2]);

		// stores the spin configuration of the system
		//initial state chosen by random no. generator above
		for (unsigned int i = 0; i < axis1; ++i)
			for (unsigned int j = 0; j < axis2; ++j)
				sitespin[i][j] = 2 * roll_coin(0, 1) - 1;

		double energy = energy_tot(sitespin);
		double mag = mag_tot(sitespin) * 1.0;

		double en_sum(0), mag_sum(0), en_sq_sum(0), mag_sq_sum(0);
		double abs_mag_sum(0), sp_heat(0), susc(0);

		unsigned int sys_size = axis1 * axis2;

		for (unsigned int i = 1; i <= N_mc; ++i)
		{
			for (unsigned int j = 1; j <= sys_size; ++j)
			{	//Now choose a random spin site, say, i

				unsigned int rnd_row, rnd_col;
				rnd_row = roll_coin(1, axis1) - 1;
				rnd_col = roll_coin(1, axis2) - 1;

				double
				    energy_diff =
				    -2 * nn_energy(sitespin, rnd_row, rnd_col);

				//Generate a random no. r such that 0 < r < 1

				double r = random_real(0, 1);
				double acc_ratio = exp(-1.0 * energy_diff / kT);

				//Spin flipped if r <= acceptance ratio
				if (r <= acc_ratio)
				{
					sitespin[rnd_row][rnd_col] *= -1;
					energy += energy_diff;
					mag += 2.0 * sitespin[rnd_row][rnd_col];
				}
			}

			//Constitute a Markov chain assuming it is given at 
			// every N steps where N = system size
			// Find sums and averages at the end of ech N steps
			//instead of every step, where subsequest configuration 
			//is not very different as it involves single spin flip

			//Given the energy of ising system at a selection of times 
			// during the simulation, we can average them to find 
			//the estimators of internal energy
			//dividing it by no. of sites gives internal energy per site

			en_sum += energy;
			en_sq_sum += energy * energy;
			mag_sum += mag;
			mag_sq_sum += mag * mag;
			abs_mag_sum += abs(mag);

			//heat capacity per spin = (<E^2> - <E>^2 )/(system size*k*T^2)
			sp_heat = en_sq_sum / i - en_sum * en_sum / (i * i);

			//susceptibility = (<M^2> - <|M|>^2 )/(k*T)
			susc =
			    mag_sq_sum / i -
			    abs_mag_sum * abs_mag_sum / (i * i);

		}

		fout << kT / J
		    << '\t' << sp_heat / (sys_size * kT * kT)
		    << '\t' << susc / kT
		    << '\t' << en_sum / N_mc
		    << '\t' << en_sq_sum / N_mc
		    << '\t' << mag_sum / N_mc << '\t' << mag_sq_sum /
		    N_mc << endl;
	}

	fout.close();
	cout << "Output is phys quantitites vs kT / J " << endl;
	cout << "Data file columns are as follows:" << endl;
	cout << "kT/J--heat capacity per spin--susceptibility--total avg energy"
	    << "--total avg energy square--total avg magnetization"
	    << "--total avg magnetization square" << endl;

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
		energy -= J * sitespin[axis1 - 1][j] * sitespin[0][j];

	for (unsigned int i = 0; i < axis1; ++i)
		energy -= J * sitespin[i][axis2 - 1] * sitespin[i][0];

	return energy;
}

//function to calculate magnetization
//for a given spin configuration

double mag_tot(array_2d sitespin)
{
	int mag = 0;

	for (unsigned int i = 0; i < axis1; i++)
		for (unsigned int j = 0; j < axis2; j++)
			mag += sitespin[i][j];

	return mag;
}

//Calculating interaction energy change for spin 
//at random site->(rnd_row,rnd_col) with its nearest neighbours
double nn_energy(array_2d sitespin, unsigned int rnd_row, unsigned int rnd_col)
{
	double nn_en = 0;

	if (rnd_row > 0 && rnd_row < axis1 - 1)
	{
		nn_en -=
		    J * sitespin[rnd_row][rnd_col] * sitespin[rnd_row -
							      1][rnd_col];
		nn_en -=
		    J * sitespin[rnd_row][rnd_col] * sitespin[rnd_row +
							      1][rnd_col];
	}

	if (rnd_col > 0 && rnd_col < axis2 - 1)
	{
		nn_en -=
		    J * sitespin[rnd_row][rnd_col] * sitespin[rnd_row][rnd_col -
								       1];
		nn_en -=
		    J * sitespin[rnd_row][rnd_col] * sitespin[rnd_row][rnd_col +
								       1];
	}

	if (rnd_row == 0)
	{
		nn_en -=
		    J * sitespin[0][rnd_col] * sitespin[axis1 - 1][rnd_col];
		nn_en -= J * sitespin[0][rnd_col] * sitespin[1][rnd_col];

	}

	if (rnd_row == axis1 - 1)
	{
		nn_en -=
		    J * sitespin[axis1 - 1][rnd_col] * sitespin[axis1 -
								2][rnd_col];
		nn_en -=
		    J * sitespin[axis1 - 1][rnd_col] * sitespin[0][rnd_col];

	}

	if (rnd_col == 0)
	{
		nn_en -=
		    J * sitespin[rnd_row][0] * sitespin[rnd_row][axis2 - 1];
		nn_en -= J * sitespin[rnd_row][0] * sitespin[rnd_row][1];

	}

	if (rnd_col == axis2 - 1)
	{
		nn_en -=
		    J * sitespin[rnd_row][axis2 - 1] * sitespin[rnd_row][axis2 -
									 2];
		nn_en -=
		    J * sitespin[rnd_row][axis2 - 1] * sitespin[rnd_row][0];

	}
	return nn_en;
}

//  g++ -Wall -O3 ising-2d-T-dep.cc -o testo

// Considering 2d ising model using ./boost multi array 
//in zero magnetic field 
//--http://www.boost.org/doc/libs/1_39_0/libs/multi_array/doc/user.html
//Metropolis algorithm employed

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

typedef boost::multi_array<int, 2> array_2d;
// typedef keyword allows you to create an alias fo a data type


// Define global scopes to use them across all functions
double J = 0;
const unsigned int axis1 = 10, axis2=10;
// above assigns length along each dimension of the 2d configuration

//Function templates
int roll_coin(int a, int b);
double random_real(int a, int b);
double energy_tot(array_2d sitespin);
double mag_tot(array_2d sitespin);
double nn_energy(array_2d sitespin,unsigned int rnd_row, unsigned int rnd_col);



int main()
{

	// Create a 2d array that is axis1 * axis2
  	array_2d sitespin(boost::extents[axis1][axis2]);

	
	// stores the spin configuration of the system
	//initial state chosen by random no. generator above
	for (unsigned int i = 0; i < axis1; ++i)
		for (unsigned int j = 0; j < axis2; ++j)
			sitespin[i][j] = 2*roll_coin(0,1) - 1;
	
	

	double kT(0), beta(0);

	cout << "Enter T (in units of k)" << endl;
	cin >> kT;
	
	beta = 1.0/kT;
	
	
	//Calculating initial energy for this configuration above
	cout << "Enter value of coupling constant J > 0 " << endl;
	cin >> J;
	
	double energy = energy_tot(sitespin);
	double mag = mag_tot (sitespin ) * 1.0;
	
	//Input from terminal no.of Monte Carlo updates we want
	unsigned int N_mc ;

	cout << "Enter no. of Monte Carlo updates N_mc" << endl;
	cout << "total steps = N_mc *  system size" << endl;
	cin >> N_mc;
	
	double en_sum(0), mag_sum(0), en_sq_sum(0), mag_sq_sum(0);
	double abs_mag_sum(0);
	
	
	ofstream fout("2d.dat"); // Opens a file for output

	unsigned int sys_size = axis1*axis2 ;
	
	for (unsigned int i = 1; i <= N_mc; ++i)
		{for (unsigned int j = 1; j <= sys_size; ++j)
			{//Now choose a random spin site, say, i

			unsigned int rnd_row , rnd_col;
			rnd_row = roll_coin(1, axis1)-1;
			rnd_col = roll_coin(1, axis2)-1;
		
			double energy_diff=-2*nn_energy(sitespin,rnd_row,rnd_col);
		

			//Generate a random no. r such that 0 < r < 1

			double r = random_real(0, 1);
			double acc_ratio = exp(-1.0 * energy_diff * beta);

			//Spin flipped if r <= acceptance ratio
			if (r <= acc_ratio) 
				{	sitespin[rnd_row] [rnd_col] *= -1;
				energy += energy_diff;
				mag +=2.0 * sitespin[rnd_row] [rnd_col];
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
		
		en_sum += energy ;
		en_sq_sum += energy * energy ;
		mag_sum += mag ;
		abs_mag_sum +=abs(mag);
		
		//heat capacity per spin = (<E^2> - <E>^2 )/(system size*k*T^2)
		double sp_heat = en_sq_sum/i - en_sum*en_sum/(i*i);
		
		//susceptibility = (<M^2> - <|M|>^2 )/(k*T)
		double susc = mag_sq_sum/i - abs_mag_sum*abs_mag_sum/(i*i);
		
		
		fout << i 
		<< '\t' << en_sum / i 
		<< '\t' << en_sq_sum / i 
		<< '\t' << sp_heat/(sys_size*kT*kT)
		<< '\t' << mag_sum / i 
		<< '\t' << susc / kT
		<< endl;
		

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
	double energy = 0;

	for (unsigned int i = 0; i < axis1-1; ++i)
	{
		for (unsigned int j = 0; j < axis2-1; ++j)
		{
			energy -= J*sitespin[i][j]*sitespin[i+1][j];
			    
			energy -= J*sitespin[i][j]*sitespin[i][j+1];
			    
		}
	}
	
	//periodic boundary conditions
	for (unsigned int j=0 ; j < axis2 ; ++j)
		energy -= J * sitespin[axis1-1][j] * sitespin[0][j];
	
	for (unsigned int i=0 ; i < axis1 ; ++i)
		energy -= J * sitespin[i][axis2-1] * sitespin[i][0];
			
	return energy ;
}



//function to calculate magnetization
//for a given spin configuration
	
double mag_tot(array_2d sitespin)
{
	int mag = 0;

	for (unsigned int i = 0; i < axis1 ; i++)
		for (unsigned int j = 0; j < axis2 ; j++) 
			mag +=sitespin[i][j];
	
	return mag ;
}



//Calculating interaction energy change for spin 
//at random site->(rnd_row,rnd_col) with its nearest neighbours
double nn_energy(array_2d sitespin,unsigned int rnd_row, unsigned int rnd_col)
{	double nn_en=0;
		
	if (rnd_row>0 && rnd_row < axis1-1)
	{	nn_en -=J*sitespin[rnd_row][rnd_col]
			*sitespin[rnd_row-1][rnd_col];
		nn_en -=J*sitespin[rnd_row][rnd_col]
			*sitespin[rnd_row+1][rnd_col];
		}
	
	if (rnd_col>0 && rnd_col < axis2-1)
	{	nn_en -=J*sitespin[rnd_row][rnd_col]
			*sitespin[rnd_row][rnd_col-1];
		nn_en -=J*sitespin[rnd_row][rnd_col]
			*sitespin[rnd_row][rnd_col+1];
	}
		
		
		
	if (rnd_row ==0)
	{	nn_en -=J*sitespin[0][rnd_col]
			*sitespin[axis1-1][rnd_col];
		nn_en -=J*sitespin[0][rnd_col]
			*sitespin[1][rnd_col];
	
	}
		
	if (rnd_row ==axis1-1)
	{
		nn_en -=J*sitespin[axis1-1][rnd_col]
			*sitespin[axis1-2][rnd_col];
		nn_en -=J*sitespin[axis1-1][rnd_col]
			*sitespin[0][rnd_col];
	
	}
				
	if (rnd_col ==0)
	{
		nn_en -=J*sitespin[rnd_row][0]
			*sitespin[rnd_row][axis2-1];
		nn_en -=J*sitespin[rnd_row][0]
			*sitespin[rnd_row][1];
	
		}
		
	if (rnd_col ==axis2-1)
	{
		nn_en -=J*sitespin[rnd_row][axis2-1]
			*sitespin[rnd_row][axis2-2];
		nn_en -=J*sitespin[rnd_row][axis2-1]
			*sitespin[rnd_row][0];
	
	}
	return nn_en;
}

//  g++ -Wall -O3 ising-2d-boost-multiarr.cc -o testo



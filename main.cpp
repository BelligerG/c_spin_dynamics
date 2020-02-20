#include "spin_dynamics_plusplus.h"
#include <iostream>
#include <Eigen/Dense>
#include <fstream>

/*TODO
 * Changing spin to 1 causes an error -- need to fix this
 * Tidy the code to make it general;
 * 	Use the interaction partners for hyperfine
 * 	Set the zeeman interaction as a whole
 * 	Use internal hamiltonian variable
 */

/*Steps
 * Exchange
 */

/*CHECKED CALCULATIONS
 * 3 SPINS Zeeman matches
 * 3 SPINS Dipolar matches
 * 6 SPINS Hyperfine 3e- & 3n+ matches
 */


int main()
{
	SpinDynamics SPD;

	SPD.SetNumberOfElectrons(3);
	//These are all the spins (nuclei too)
	float spins [3] = {0.5, 0.5, 0.5};

	//Takes in the spins and the total number of spins
	SPD.spinsToSpinOperators(spins, sizeof(spins)/sizeof(spins[0]));

	//Currently saving the spin_operators outside of the object
	std::vector<std::vector<Eigen::MatrixXcd>> spin_ops = SPD.GetSpinOperators();

	SPD.SetMagneticField({0, 0, 1.0});
	//std::vector<std::pair<int, int>> interaction_partners = {{0, 3}, {1, 4}, {2, 5}};







	//Change this to be anisotropic?
	//Eigen::MatrixXcd h_hyperfine = SPD.hyperfine(spin_ops[0], spin_ops[3], 10);
	//Eigen::MatrixXcd h_hyperfine = SPD.hyperfine(spin_ops[interaction_partners[0].first], spin_ops[interaction_partners[0].second], 10.3172);
	//h_hyperfine += SPD.hyperfine(spin_ops[interaction_partners[1].first], spin_ops[interaction_partners[1].second], 10.3172);
	//h_hyperfine += SPD.hyperfine(spin_ops[interaction_partners[2].first], spin_ops[interaction_partners[2].second], 10.3172);

	//std::cout << zeeman+h_hyperfine << std::endl;

	//Must have a coordinate for each electron
	std::vector<Eigen::Vector3d> coordinates;
	coordinates.push_back({0, 0, -4.5});
	coordinates.push_back({0, 0, 4.5});
	coordinates.push_back({0, 0, 9.0});

	//Make sure can calculate the individual dipolar contribution
	Eigen::MatrixXcd h_dipolar = SPD.calculateDipolar(coordinates);

	//Save all the hamiltonian sections as we go along, then save the total hamiltonian and reset the sections/sum up the hamiltonian as we go?
	Eigen::MatrixXcd hamiltonian = h_dipolar;//+h_hyperfine;//+zeeman;


	std::vector<Eigen::Vector3d> distances = SPD.calculateDistances(coordinates);


	//Set these in the functions? But need to check how this is affected by different combinations of electrons
	double kSc = 0.01;

	SPD.SetBeta(1.4);
	SPD.SetRadicalRadius(4.5);
	SPD.SetkS0(0.2);
	Eigen::MatrixXcd K1 = SPD.calculateK1(coordinates);

	std::vector<std::vector<double>> yields;

	omp_lock_t lck;
	omp_init_lock(&lck);

	#pragma omp parallel for
	for(int B=0; B<10; B+=2){
		SpinDynamics SPD_loop  = SPD;
		double z_field = B;
		SPD_loop.SetMagneticField({0, 0, z_field});

		Eigen::MatrixXcd zeeman = SPD_loop.calculateZeeman();

		Eigen::MatrixXcd hamiltonian_loop = hamiltonian+zeeman;
		double yield = SPD_loop.singletYield(hamiltonian_loop, K1, kSc);

		omp_set_lock(&lck);
			yields.push_back({0, 0, z_field, yield});
		omp_unset_lock(&lck);
		
	}
	omp_destroy_lock(&lck);

	int number_of_yields = yields.size();
	std::ofstream outFile("results.txt");
	outFile << "Magnetic Field x, y, z, Yield\n";
	for (int i=0; i<number_of_yields; i++){
		outFile << yields[i][0] << "," << yields[i][1] << "," << yields[i][2] << "," << yields[i][3] << "\n";
	}
	
	//double yield = SPD.singletYield(hamiltonian, K1, kSc);
	//std::cout << yield << std::endl;


	return 0;
}


#include "spin_dynamics_plusplus.h"
#include <iostream>
#include <Eigen/Dense>


/*TODO
 * Link it to github, just in case!
 * Update code so that we can iterate over the vectors rather than having to explicitly mention the tuples
 * Calculate the Hamiltonian sections
 * Initial density operator
 * Calculate the singlet yield
 */

/*Steps
 * Exchange
 * Yields
 */

//Dipolar and exchange distances need to be calculated for all radical pairs

int main()
{
	SpinDynamics SPD;

	//This needs to be dynamic eventually could this be done in spinsToSpinOperators?
	SPD.SetSizeOfMatrix(4);
	std::cout << SPD.GetSizeOfMatrix() << std::endl;

	float spins [2] = {0.5, 0.5};

	SPD.spinsToSpinOperators(spins);
	std::vector<std::vector<Eigen::MatrixXcd>> spin_ops = SPD.GetSpinOperators();
        std::cout << spin_ops[0][2] << std::endl;



	//Calculating the zeeman interaction
	double magnetic_field [3] = {0, 0, 1.0};
	Eigen::MatrixXcd zeeman = SPD.zeeman(spin_ops[0], magnetic_field);
	zeeman += SPD.zeeman(spin_ops[1], magnetic_field);

	std::cout << zeeman << std::endl;

	Eigen::MatrixXcd h_hyperfine = SPD.hyperfine(spin_ops[0], spin_ops[1], 10.3172);

	std::cout << h_hyperfine << std::endl;

	Eigen::MatrixXcd h_dipole = SPD.dipolar(spin_ops[0], spin_ops[1], Eigen::Vector3d(1, 2, 3));

	std::cout << h_dipole << std::endl;


	std::vector<Eigen::Vector3d> coordinates;
	coordinates.push_back({0, 0, 0});
	coordinates.push_back({0, 0, 4.5});
	coordinates.push_back({0, 0, 9});
	coordinates.push_back({0, 0, 20});

	std::vector<Eigen::Vector3d> distances = SPD.calculateDistances(coordinates);

	int size = distances.size();
	for (int i=0;i<size;i++){
		std::cout << distances[i] << std::endl;
	}

	return 0;
}


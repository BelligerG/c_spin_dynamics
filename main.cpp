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
//Dipolar done

int main()
{
	SpinDynamics SPD;

	//This needs to be dynamic eventually could this be done in spinsToSpinOperators?
	//SPD.SetSizeOfMatrix(4);
	SPD.SetNumberOfElectrons(2);

	float spins [3] = {0.5, 0.5, 0.5};

	SPD.spinsToSpinOperators(spins, sizeof(spins)/sizeof(spins[0]));

	//int size_of_matrix = SPD.GetSizeOfMatrix();

	std::vector<std::vector<Eigen::MatrixXcd>> spin_ops = SPD.GetSpinOperators();
        //std::cout << spin_ops[0][2] << std::endl;

	double magnetic_field [3] = {0, 0, 1.0};
	Eigen::MatrixXcd zeeman = SPD.zeeman(spin_ops[0], magnetic_field);
	zeeman += SPD.zeeman(spin_ops[1], magnetic_field);
	Eigen::MatrixXcd h_hyperfine = SPD.hyperfine(spin_ops[0], spin_ops[2], 10.3172);

	std::vector<Eigen::Vector3d> coordinates;
	coordinates.push_back({0, 0, 0});
	coordinates.push_back({0, 0, 4.5});

	Eigen::MatrixXcd h_dipolar = SPD.calculateDipolar(coordinates);

	return 0;
}


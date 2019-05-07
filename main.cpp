#include "spin_dynamics_plusplus.h"
#include <iostream>
#include <Eigen/Dense>


/*TODO
 * Link it to github, just in case!
 * Update code so that we can iterate over the vectors rather than having to explicitly mention the tuples
 */

/*Steps
 * Exchange
 */

int main()
{
	SpinDynamics SPD;


	SPD.SetNumberOfElectrons(2);
	//These are all the spins (nuclei too)
	float spins [3] = {0.5, 0.5, 0.5};

	//Takes in the spins and the total number of spins
	SPD.spinsToSpinOperators(spins, sizeof(spins)/sizeof(spins[0]));

	//Currently saving the spin_operators outside of the object
	std::vector<std::vector<Eigen::MatrixXcd>> spin_ops = SPD.GetSpinOperators();

	std::vector<double> field = {0, 0, 1.0};
	SPD.SetMagneticField(field);





	//Make a total zeeman that takes in the number of electrons and the spins and runs the individual zeeman functions
	Eigen::MatrixXcd zeeman = SPD.zeeman(spin_ops[0]);
	zeeman += SPD.zeeman(spin_ops[1]);




	//Change this to be anisotropic?
	Eigen::MatrixXcd h_hyperfine = SPD.hyperfine(spin_ops[0], spin_ops[2], 10.3172);

	//Must have a coordinate for each electron
	std::vector<Eigen::Vector3d> coordinates;
	coordinates.push_back({0, 0, -4.5});
	coordinates.push_back({0, 0, 4.5});

	//Make sure can calculate the individual dipolar contribution
	Eigen::MatrixXcd h_dipolar = SPD.calculateDipolar(coordinates);

	//Save all the hamiltonian sections as we go along, then save the total hamiltonian and reset the sections/sum up the hamiltonian as we go?
	Eigen::MatrixXcd hamiltonian = zeeman+h_hyperfine+h_dipolar;



	//Set these in the functions? But need to check how this is affected by different combinations of electrons
	double kS0 = 0.2;
	double kSc = 0.01;

	//Pass the electron indices here rather than the spin ops, then save as a function in the class
	Eigen::MatrixXcd singlet_projector = SPD.singletProjector(spin_ops[0], spin_ops[1]);
	std::cout << singlet_projector << std::endl;
	


	double eScale = SPD.expScaling(1.4, 4.5, 9.0);
	
	Eigen::MatrixXcd K1 = 0.5*kS0*eScale*singlet_projector;

	double yield = SPD.singletYield(hamiltonian, K1, kSc);
	std::cout << yield << std::endl;

	return 0;
}


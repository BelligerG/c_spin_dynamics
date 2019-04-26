#ifndef SPIN_DYNAMICS_PLUSPLUS_H
#define SPIN_DYNAMICS_PLUSPLUS_H

#include <iostream>
#include <Eigen/Dense>
#include <complex>
#include <vector>
#include <math.h>
#include <tuple>

class SpinDynamics
{
public:

	SpinDynamics();
	//Getters and Setters
	void SetSizeOfMatrix(int);
	int GetSizeOfMatrix();
	void SetSpinOperators(std::vector<std::vector<Eigen::MatrixXcd>>);
	std::vector<std::vector<Eigen::MatrixXcd>> GetSpinOperators();
	void SetNumberOfElectrons(int);
	int GetNumberOfElectrons();

	//Setting up spin operators
	void spinsToSpinOperators(float [], int);

	//Hamiltonian sections
	Eigen::MatrixXcd zeeman(std::vector<Eigen::MatrixXcd>, double [3]);
	Eigen::MatrixXcd hyperfine(std::vector<Eigen::MatrixXcd>, std::vector<Eigen::MatrixXcd>, double);
	Eigen::MatrixXcd calculateDipolar(std::vector<Eigen::Vector3d>);

private:
	//Internal variables
	int size_of_matrix;
	double magnetic_field [3];
	int number_of_electrons;
	std::vector<std::vector<Eigen::MatrixXcd>> spin_operators;

	//Functions used for the internal workings
	Eigen::MatrixXcd kroneckerProductComplex(Eigen::MatrixXcd a, Eigen::MatrixXcd b);
	std::vector<float> generateMlValues(float spin);
	std::vector<Eigen::MatrixXcd> spinOperatorCorrectSpace(int pos, int numberOfSpins, double spin);
	std::vector<Eigen::MatrixXcd> deriveSpinOperator(float spin);
	Eigen::MatrixXcd calculateCombinationMatrix(Eigen::MatrixXcd, std::vector<Eigen::MatrixXcd>, std::vector<Eigen::MatrixXcd>);
	Eigen::MatrixXcd kroneckerProductComplexSlow(Eigen::MatrixXcd, Eigen::MatrixXcd);
	std::vector<Eigen::Vector3d> calculateDistances(std::vector<Eigen::Vector3d>);
	
	//Individual hamiltonian matrices
	Eigen::MatrixXcd dipolar(std::vector<Eigen::MatrixXcd>, std::vector<Eigen::MatrixXcd>, Eigen::Vector3d r);
};

#endif

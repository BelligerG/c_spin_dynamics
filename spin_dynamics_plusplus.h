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
	void SetSizeOfMatrix(int);
	int GetSizeOfMatrix();
	void SetSpinOperators(std::vector<std::vector<Eigen::MatrixXcd>>);
	std::vector<std::vector<Eigen::MatrixXcd>> GetSpinOperators();

	//Eigen::MatrixXcd zeeman(Eigen::MatrixXcd [3], double [3]);
	Eigen::MatrixXcd zeeman(std::vector<Eigen::MatrixXcd>, double [3]);
	void spinsToSpinOperators(float []);

private:
	int size_of_matrix;
	double magnetic_field [3];

	std::vector<std::vector<Eigen::MatrixXcd>> spin_operators;
	Eigen::MatrixXcd kroneckerProductComplex(Eigen::MatrixXcd a, Eigen::MatrixXcd b);
	std::vector<float> generateMlValues(float spin);
	std::vector<Eigen::MatrixXcd> spinOperatorCorrectSpace(int pos, int numberOfSpins, double spin);
	std::vector<Eigen::MatrixXcd> deriveSpinOperator(float spin);
};

#endif

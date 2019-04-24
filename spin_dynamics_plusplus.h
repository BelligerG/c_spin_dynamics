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

	void spinsToSpinOperators(float []);

	Eigen::MatrixXcd zeeman(std::vector<Eigen::MatrixXcd>, double [3]);
	Eigen::MatrixXcd hyperfine(std::vector<Eigen::MatrixXcd>, std::vector<Eigen::MatrixXcd>, double);
	Eigen::MatrixXcd dipolar(std::vector<Eigen::MatrixXcd>, std::vector<Eigen::MatrixXcd>, Eigen::Vector3d r);

	std::vector<Eigen::Vector3d> calculateDistances(std::vector<Eigen::Vector3d>);

private:
	int size_of_matrix;
	double magnetic_field [3];

	std::vector<std::vector<Eigen::MatrixXcd>> spin_operators;
	Eigen::MatrixXcd kroneckerProductComplex(Eigen::MatrixXcd a, Eigen::MatrixXcd b);
	std::vector<float> generateMlValues(float spin);
	std::vector<Eigen::MatrixXcd> spinOperatorCorrectSpace(int pos, int numberOfSpins, double spin);
	std::vector<Eigen::MatrixXcd> deriveSpinOperator(float spin);

	Eigen::MatrixXcd calculateCombinationMatrix(Eigen::MatrixXcd, std::vector<Eigen::MatrixXcd>, std::vector<Eigen::MatrixXcd>);
	Eigen::MatrixXcd kroneckerProductComplexSlow(Eigen::MatrixXcd, Eigen::MatrixXcd);
};

#endif

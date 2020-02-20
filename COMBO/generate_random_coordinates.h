#ifndef GENERATE_RANDOM_COORDINATES_H
#define GENERATE_RANDOM_COORDINATES_H

#include <iostream>
#include <vector>
#include <random>
#include "dc.h"

class GenerateCoordinates
{
public:

	GenerateCoordinates();
	
	void SetRMax(double);
	void SetRGridStep(double);
	void SetPhiGridStep(double);
	void SetXMax(double);
	void SetYMax(double);
	void SetZMax(double);
	void SetXGridStep(double);
	void SetYGridStep(double);

	void GeneratePossibleRValues();
	std::vector<double> GetPossibleRValues();

	void GeneratePossiblePhiValues();
	std::vector<double> GetPossiblePhiValues();

	void GeneratePossibleXValues();
	void GeneratePossibleYValues();

	std::vector<double> RandomPolarCoordinates();
	std::vector<double> RandomCartesianCoordinates();
	std::vector<double> RandomCoordinateAtR(std::vector<double>, double);
	std::vector<double> RandomCoordinateNotInr(std::vector<std::vector<double>>, double);

	std::vector<double> RandomContinuousCartesianCoordinates();
	std::vector<double> RandomContinuousCoordinateAtR(std::vector<double>, double);
	std::vector<double> RandomContinuousCoordinateNotInr(std::vector<std::vector<double>>, double);
	//Thread safe versions - the final argument is the thread number
	std::vector<double> RandomContinuousCoordinateAtR(std::vector<double>, double, int);
	std::vector<double> RandomContinuousCartesianCoordinates(int);
	std::vector<double> RandomContinuousCoordinateNotInr(std::vector<std::vector<double>>, double, int);

	std::vector<double> RandomContinuousCoordinateAtR3D(std::vector<double>, double, int);
	std::vector<double> RandomContinuousCartesianCoordinates3D(int);
	std::vector<double> RandomContinuousCoordinateNotInr3D(std::vector<std::vector<double>>, double, int);

	GenerateCoordinates(const GenerateCoordinates& gcd);

private:
	double r_max;
	double r_grid_step;
	std::vector<double> possible_r_values;
	double phi_grid_step;
	std::vector<double> possible_phi_values;

	double x_max;
	double y_max;
	double z_max;
	double x_grid_step;
	double y_grid_step;
	std::vector<double> possible_x_values;
	std::vector<double> possible_y_values;

	int RandomInt(int);
	double RandomDouble(int);
	double RandomDouble(int, int);

	std::vector<double> ARange(double, double, double);
};

#endif

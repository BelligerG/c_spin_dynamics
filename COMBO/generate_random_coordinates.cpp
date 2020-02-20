#include "generate_random_coordinates.h"

GenerateCoordinates::GenerateCoordinates(){}

void GenerateCoordinates::SetRMax( double rmax ){ r_max = rmax; }

void GenerateCoordinates::SetRGridStep( double gridstep ){ r_grid_step = gridstep; }

void GenerateCoordinates::SetPhiGridStep( double gridstep ){ phi_grid_step = gridstep; }

void GenerateCoordinates::SetXMax( double xmax ){ x_max = xmax; }
void GenerateCoordinates::SetYMax( double ymax ){ y_max = ymax; }
void GenerateCoordinates::SetZMax( double zmax ){ z_max = zmax; }
void GenerateCoordinates::SetXGridStep( double xgridstep ){ x_grid_step = xgridstep; }
void GenerateCoordinates::SetYGridStep( double ygridstep ){ y_grid_step = ygridstep; }

void GenerateCoordinates::GeneratePossibleRValues(){ possible_r_values = ARange(0, r_max, r_grid_step);  }
std::vector<double> GenerateCoordinates::GetPossibleRValues(){ return possible_r_values; }

void GenerateCoordinates::GeneratePossibleXValues(){ possible_x_values = ARange(0, x_max, x_grid_step);  }
void GenerateCoordinates::GeneratePossibleYValues(){ possible_y_values = ARange(0, y_max, y_grid_step);  }


void GenerateCoordinates::GeneratePossiblePhiValues(){ possible_phi_values = ARange(0, 2*M_PI, phi_grid_step); }
std::vector<double> GenerateCoordinates::GetPossiblePhiValues(){ return possible_phi_values; }


std::vector<double> GenerateCoordinates::ARange(double min, double max, double step){
	std::vector<double> values;
	for (double val=min; val<=max; val+=step){
		values.push_back(val);
	}

	return values;
}

int GenerateCoordinates::RandomInt(int max_value){
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> dis(0, max_value);

	return dis(gen);
}

double GenerateCoordinates::RandomDouble(int max_value){
	std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dis(0, max_value);

	return dis(gen);
}

//thread safe version of the above function using code from: Dynamic Creation of Pseudorandom Number Generators - Matsumoto & Nishimura
double GenerateCoordinates::RandomDouble(int max_value, int thread){
	std::random_device rd;

	mt_struct *mts;

	mts = get_mt_parameter_id_st(32, 521, thread, rd());
	sgenrand_mt(rd(), mts);
	uint32_t num = genrand_mt(mts);

	free_mt_struct(mts);

	return ((double)num/UINT32_MAX)*max_value;
}



std::vector<double> GenerateCoordinates::RandomCartesianCoordinates(){
	std::vector<double> random_coordinates;

	int number_of_x_values = possible_x_values.size();
	int number_of_y_values = possible_y_values.size();

	int pos_x = RandomInt(number_of_x_values-1);
	int pos_y = RandomInt(number_of_y_values-1);

	double x = possible_x_values[pos_x];
	double y = possible_y_values[pos_y];

	random_coordinates.push_back(x);
	random_coordinates.push_back(y);

	return random_coordinates;

}

std::vector<double> GenerateCoordinates::RandomContinuousCartesianCoordinates(){
	std::vector<double> random_coordinates;
	
	double x = RandomDouble(x_max);
	double y = RandomDouble(y_max);

	random_coordinates.push_back(x);
	random_coordinates.push_back(y);

	return random_coordinates;
}

std::vector<double> GenerateCoordinates::RandomContinuousCartesianCoordinates(int thread){
	std::vector<double> random_coordinates;
	
	double x = RandomDouble(x_max, thread);
	double y = RandomDouble(y_max, thread);

	random_coordinates.push_back(x);
	random_coordinates.push_back(y);

	return random_coordinates;
}

std::vector<double> GenerateCoordinates::RandomContinuousCartesianCoordinates3D(int thread){
	std::vector<double> random_coordinates;
	
	double x = RandomDouble(x_max, thread);
	double y = RandomDouble(y_max, thread);
	double z = RandomDouble(z_max, thread);

	random_coordinates.push_back(x);
	random_coordinates.push_back(y);
	random_coordinates.push_back(z);

	return random_coordinates;
}


std::vector<double> GenerateCoordinates::RandomPolarCoordinates(){
	std::vector<double> random_coordinates;

	int number_of_r_values = possible_r_values.size();
	int number_of_phi_values = possible_phi_values.size();

	int pos_r = RandomInt(number_of_r_values-1);
	int pos_phi = RandomInt(number_of_phi_values-1);

	double r = possible_r_values[pos_r];
	double phi = possible_phi_values[pos_phi];

	random_coordinates.push_back(r);
	random_coordinates.push_back(phi);

	return random_coordinates;

}

GenerateCoordinates::GenerateCoordinates(const GenerateCoordinates& gcd){
	r_max = gcd.r_max;
	r_grid_step = gcd.r_grid_step;
	phi_grid_step = gcd.phi_grid_step;

	possible_phi_values = gcd.possible_phi_values;
	possible_r_values = gcd.possible_r_values;

	x_max = gcd.x_max;
	y_max = gcd.y_max;
	z_max = gcd.z_max;
	x_grid_step = gcd.x_grid_step;
	y_grid_step = gcd.y_grid_step;

	possible_x_values = gcd.possible_x_values;
	possible_y_values = gcd.possible_y_values;

}

std::vector<double> GenerateCoordinates::RandomContinuousCoordinateAtR(std::vector<double> coord, double R){

	double x = -1;
	double y = -1;
	while (x < 0 || x > x_max || y < 0 || y > y_max){
		double phi = RandomDouble(M_PI*2);

		x = R*std::cos(phi)+coord[0];
		y = R*std::sin(phi)+coord[1];
	}

	std::vector<double> new_coords;
	new_coords.push_back(x);
	new_coords.push_back(y);

	return new_coords;
}

//Slightly different handling is required of theta and phi to generate truley random variables on a sphere.
//See: http://mathworld.wolfram.com/SpherePointPicking.html for more details
std::vector<double> GenerateCoordinates::RandomContinuousCoordinateAtR3D(std::vector<double> coord, double R, int thread){

	double x = -1;
	double y = -1;
	double z = -1;
	while (x < 0 || x > x_max || y < 0 || y > y_max || z < 0 || z > z_max){

		//THIS NEEDS TO BE UPDATED WITH THE CORRECT SPP CODE

		double u = RandomDouble(1, thread)*2.0-1.0;
		double theta = RandomDouble(1, thread)*2.0*M_PI;
		double norm = std::sqrt(std::abs(1-std::pow(u, 2)));
		x = R*norm*std::cos(theta)+coord[0];
		y = R*norm*std::sin(theta)+coord[1];
		z = R*u+coord[2];

	}

	std::vector<double> new_coords;
	new_coords.push_back(x);
	new_coords.push_back(y);
	new_coords.push_back(z);

	return new_coords;
}


std::vector<double> GenerateCoordinates::RandomContinuousCoordinateAtR(std::vector<double> coord, double R, int thread){

	double x = -1;
	double y = -1;
	while (x < 0 || x > x_max || y < 0 || y > y_max){
		double phi = RandomDouble(2*M_PI, thread);

		x = R*std::cos(phi)+coord[0];
		y = R*std::sin(phi)+coord[1];
	}

	std::vector<double> new_coords;
	new_coords.push_back(x);
	new_coords.push_back(y);

	return new_coords;
}

std::vector<double> GenerateCoordinates::RandomCoordinateAtR(std::vector<double> coord, double R){

	SetPhiGridStep((x_grid_step+y_grid_step)/2);
	GeneratePossiblePhiValues();

	int number_of_phi_values = possible_phi_values.size();

	double x = -1;
	double y = -1;
	while (x < 0 || x > x_max || y < 0 || y > y_max){
		int pos_phi = RandomInt(number_of_phi_values-1);
		double phi = possible_phi_values[pos_phi];

		x = R*std::cos(phi)+coord[0];
		y = R*std::sin(phi)+coord[1];
	}

	std::vector<double> new_coords;
	new_coords.push_back(x);
	new_coords.push_back(y);

	return new_coords;
}

std::vector<double> GenerateCoordinates::RandomCoordinateNotInr(std::vector<std::vector<double>> all_coords, double r){
	std::vector<double> coord;

	bool too_close = true;	
	while (too_close == true){
		too_close = false;
	                                                                                                                  
		coord = RandomCartesianCoordinates();
		int number_of_coords = all_coords.size();	
	                                                                                                                  
		for (int i=0; i<number_of_coords; i++){
			double dist = std::sqrt(std::pow(coord[0]-all_coords[i][0], 2) + std::pow(coord[1]-all_coords[i][1], 2));
			if (dist < r || too_close == true){
				too_close = true;
			}
		}
	}
	return coord;
}

std::vector<double> GenerateCoordinates::RandomContinuousCoordinateNotInr(std::vector<std::vector<double>> all_coords, double r){
	std::vector<double> coord;

	bool too_close = true;	
	while (too_close == true){
		too_close = false;
	                                                                                                                  
		coord = RandomContinuousCartesianCoordinates();
		int number_of_coords = all_coords.size();	
	                                                                                                                  
		for (int i=0; i<number_of_coords; i++){
			double dist = std::sqrt(std::pow(coord[0]-all_coords[i][0], 2) + std::pow(coord[1]-all_coords[i][1], 2));
			if (dist < r || too_close == true){
				too_close = true;
			}
		}
	}
	return coord;
}

std::vector<double> GenerateCoordinates::RandomContinuousCoordinateNotInr(std::vector<std::vector<double>> all_coords, double r, int thread){
	std::vector<double> coord;

	bool too_close = true;	
	while (too_close == true){
		too_close = false;
	                                                                                                                  
		coord = RandomContinuousCartesianCoordinates(thread);
		int number_of_coords = all_coords.size();	
	                                                                                                                  
		for (int i=0; i<number_of_coords; i++){
			double dist = std::sqrt(std::pow(coord[0]-all_coords[i][0], 2) + std::pow(coord[1]-all_coords[i][1], 2));
			if (dist < r || too_close == true){
				too_close = true;
			}
		}
	}
	return coord;
}

std::vector<double> GenerateCoordinates::RandomContinuousCoordinateNotInr3D(std::vector<std::vector<double>> all_coords, double r, int thread){
	std::vector<double> coord;

	bool too_close = true;	
	while (too_close == true){
		too_close = false;
	                                                                                                                  
		coord = RandomContinuousCartesianCoordinates3D(thread);
		int number_of_coords = all_coords.size();	
	                                                                                                                  
		for (int i=0; i<number_of_coords; i++){
			double dist = std::sqrt(std::pow(coord[0]-all_coords[i][0], 2) + std::pow(coord[1]-all_coords[i][1], 2) + std::pow(coord[2]-all_coords[i][2], 2));
			if (dist < r || too_close == true){
				too_close = true;
			}
		}
	}
	return coord;
}

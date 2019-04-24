#include "spin_dynamics_plusplus.h"

SpinDynamics::SpinDynamics(){}


void SpinDynamics::SetSizeOfMatrix(int size){ size_of_matrix = size; }
int SpinDynamics::GetSizeOfMatrix(){ return size_of_matrix; }

void SpinDynamics::SetSpinOperators(std::vector<std::vector<Eigen::MatrixXcd>> spin_ops){ spin_operators = spin_ops; }
std::vector<std::vector<Eigen::MatrixXcd>> SpinDynamics::GetSpinOperators(){ return spin_operators; }

//Generates the spin operator for a spin particle
std::vector<Eigen::MatrixXcd> SpinDynamics::deriveSpinOperator(float spin){

	std::complex<double> i(0.0,1.0);

	std::vector<float> ml = generateMlValues(spin);

	int number_of_spins = ml.size();
	Eigen::VectorXd basis [number_of_spins];
	

	for(int i=0; i<number_of_spins; i++){
		basis[i] = Eigen::VectorXd::Zero(number_of_spins);
		basis[i][i] = 1;
	}

	Eigen::MatrixXcd spin_z, spin_p, spin_m, spin_x, spin_y;
	spin_z = spin_p = spin_m = spin_x = spin_y =  Eigen::MatrixXcd::Zero(number_of_spins, number_of_spins);

	for(int row=0; row<spin_z.rows(); row++){
		for(int col=0; col<spin_z.cols(); col++){
			spin_z(row, col) = ml[col]*basis[row].dot(basis[col]);

			float Sp = sqrt(spin*(spin+1) - ml[col]*(ml[col]+1));
			if(col-1!=-1){
				spin_p(row, col) = Sp*basis[row].dot(basis[col-1]);
			} else {
				spin_p(row, col) = 0;
			}

			float Sm = sqrt(spin*(spin+1) - ml[col]*(ml[col]-1));

                        if(col+1<number_of_spins){
                        	spin_m(row, col) = Sm*basis[row].dot(basis[col+1]);
                        } else {
                        	spin_m(row, col) = 0;
                        }

		}
	}


	spin_x = 0.5*(spin_p+spin_m);
	spin_y = 0.5*i*(spin_p-spin_m);

	std::vector<Eigen::MatrixXcd> spin_op;
	spin_op.push_back(spin_x);
	spin_op.push_back(spin_y);
	spin_op.push_back(spin_z);

	return spin_op;
}

//Puts the spin operator into the correct Hilbert space by calculating the tensor product
//Arguments needed are the spin operator type ("x", "y", "z"), the position in the chain of tensor products
//and the total number of spins present.
std::vector<Eigen::MatrixXcd> SpinDynamics::spinOperatorCorrectSpace(int pos, int numberOfSpins, double spin){

	std::vector<Eigen::MatrixXcd> ops = SpinDynamics::deriveSpinOperator(spin);
	
	Eigen::MatrixXcd matLeft;
	Eigen::MatrixXcd matRight;
	int dimension;

	for (int i =0; i<3; i++){
		
		dimension = pow(2, pos);
		matLeft = Eigen::MatrixXcd::Identity(dimension, dimension);
                ops[i] = SpinDynamics::kroneckerProductComplex(matLeft, ops[i]);

		dimension = pow(2,(numberOfSpins-(pos+1)));
                matRight = Eigen::MatrixXcd::Identity(dimension, dimension);
                ops[i] = SpinDynamics::kroneckerProductComplex(ops[i], matRight);
	}

	return ops;
}


//Calculates the tensor product of 2 matrices, a and b
Eigen::MatrixXcd SpinDynamics::kroneckerProductComplex(Eigen::MatrixXcd a, Eigen::MatrixXcd b){

	Eigen::MatrixXcd c=Eigen::MatrixXcd::Zero(a.rows()*b.rows(), a.cols()*b.cols());
	for(int i=0; i<a.cols();i++){
		for(int j=0; j<a.rows();j++){
			c.block(i*b.rows(), j*b.cols(), b.rows(), b.cols()) = a(i,j)*b;
		}
	}

        return c;
}

//If given S will generate m_l using the rules -S, -S+1, ..., S-1, S
std::vector<float> SpinDynamics::generateMlValues(float spin){
	std::vector<float> ml;

        float ml_value=-1.0*spin;
        ml.push_back(ml_value);
        while (ml_value<spin){
                ml_value+=1.0;
		ml.insert(ml.begin(), ml_value);
        }
	return ml;

}

//Takes in an array of spins and returns a vector of spin operators
//Does this need to be redesined to calculate derive spin operators only once?
void SpinDynamics::spinsToSpinOperators(float spins []){

	int number_of_spins = sizeof(spins)/sizeof(spins[0]);
	std::vector<std::vector<Eigen::MatrixXcd> > spin_ops;
	for (int spin_counter = 0; spin_counter<number_of_spins; spin_counter++){
		spin_ops.push_back(spinOperatorCorrectSpace(spin_counter, number_of_spins, spins[spin_counter]));
	}
	spin_operators = spin_ops;
}

Eigen::MatrixXcd SpinDynamics::zeeman(std::vector<Eigen::MatrixXcd> electron_spin_ops, double magnetic_field[3]){

	Eigen::MatrixXcd h_zeeman = Eigen::MatrixXcd::Zero(size_of_matrix, size_of_matrix);
	for(int i=0; i<3; i++){
 
		Eigen::MatrixXcd zeeman = -1*0.17608597087294475*magnetic_field[i]*(electron_spin_ops[i]);
		h_zeeman += zeeman;
        }

        return h_zeeman;

}

Eigen::MatrixXcd SpinDynamics::hyperfine(std::vector<Eigen::MatrixXcd>  electron_spin_ops, std::vector<Eigen::MatrixXcd> nuclear_spin_ops, double hyperfine_constant_A){

	Eigen::Matrix3cd hyperfine_matrix = hyperfine_constant_A*0.002*M_PI*Eigen::Matrix3cd::Identity(3,3);

	Eigen::MatrixXcd h_hyperfine = calculateCombinationMatrix(hyperfine_matrix, electron_spin_ops, nuclear_spin_ops);

	return h_hyperfine;
}

Eigen::MatrixXcd SpinDynamics::calculateCombinationMatrix(Eigen::MatrixXcd constants_matrix, std::vector<Eigen::MatrixXcd> spin1, std::vector<Eigen::MatrixXcd> spin2){
	int size_of_matrix = spin1[0].cols();
	Eigen::MatrixXcd combination_matrix = Eigen::MatrixXcd::Zero(size_of_matrix,size_of_matrix);

	for(int i=0; i<3; i++){
                for(int j=0; j<3; j++){
                        Eigen::MatrixXcd spin_op_dot_prod = spin1[i].adjoint()*spin2[j];
                        combination_matrix+=constants_matrix(i,j)*spin_op_dot_prod;
                }       
        }

	return combination_matrix;
}


//Currently takes in distance, it would be better to give the coordinates and calculate the distances
Eigen::MatrixXcd SpinDynamics::dipolar(std::vector<Eigen::MatrixXcd> spin1, std::vector<Eigen::MatrixXcd> spin2, Eigen::Vector3d r){

	Eigen::MatrixXcd dipolar(3,3);

	//calculates the numerator for the dipolar interaction (lots of constants and unit conversions)
	double dr3 = -4*M_PI*1e-7 * pow((2.0023193043617 * 9.27400968e-24), 2)/(4*M_PI*1e-30)/6.62606957e-34/1e6;
	double r_norm = r.norm();
	double d = dr3/pow(r_norm,3);
	Eigen::Vector3d e = r/r_norm;


	Eigen::MatrixXcd A = d*(3*kroneckerProductComplexSlow(e, e.transpose())-Eigen::MatrixXcd::Identity(3,3)) *2e-3*M_PI;

	//Calculates the dipolar matrix, a 3x3 matrix with x*x, x*y, x*z on the first row etc...
	Eigen::MatrixXcd h_dipolar = calculateCombinationMatrix(A, spin1, spin2);
	return h_dipolar;
}

Eigen::MatrixXcd SpinDynamics::kroneckerProductComplexSlow(Eigen::MatrixXcd a, Eigen::MatrixXcd b){

        Eigen::MatrixXcd c(a.rows()*b.rows(), a.cols()*b.cols());
        for(int i=0; i<a.rows();i++){
                for(int j=0; j<a.cols();j++){
                        for(int k=0; k<b.rows();k++){
                                for(int l=0; l<b.cols();l++){
                                        c(b.rows()*i+k, b.cols()*j+l) = a(i,j)*b(k,l);
                                }
                        }
                }
        }
        
        return c; 
 
}

//Calculates the distances between points, the returned values are in the format (for a 3 coord system)
// [01, 02, 12]
// And for 4 coordinates:
// [01, 02, 03, 12, 13, 23]
// This is due to the nature of the recursion (meaning no distance is calculated twice)
std::vector<Eigen::Vector3d> SpinDynamics::calculateDistances(std::vector<Eigen::Vector3d> coordinates){
	int number_of_electrons_left = coordinates.size();
	std::vector<Eigen::Vector3d> distances;

	if (number_of_electrons_left == 2){
		distances.push_back(coordinates[0]-coordinates[1]);
		return distances;
	} else {
		//Slicing the vector
		std::vector<Eigen::Vector3d>::const_iterator begin = coordinates.begin();
		std::vector<Eigen::Vector3d>::const_iterator last = coordinates.begin() + coordinates.size();
		std::vector<Eigen::Vector3d> nu_coordinates(begin+1, last);

		distances = calculateDistances(nu_coordinates);

		for (int coord_index=0; coord_index<number_of_electrons_left-1; coord_index++){
			distances.push_back(coordinates[0] - coordinates[coord_index+1]);
		}
	}

	return distances;
}

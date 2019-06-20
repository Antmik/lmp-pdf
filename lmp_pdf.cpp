//Examples:
//lmp_stress filename dump.atom NVT n_steps 10 n_bins 100 kernel gaussian 0.1 density momentum viscosity
//lmp_stress filename dump.atom NVT n_steps 10 n_bins 100 kernel gaussian 0.1 all

#include <iostream>
#include <fstream>
#include <string>
#include <armadillo>
#include <math.h>
#include <time.h>

double cubic_bond_function(const double pos, const double R_cut ) {
    if (pos >= 0.0 && pos<R_cut) {
        return ( pos - pos*pos*pos/R_cut/R_cut + pos*pos*pos*pos/2/R_cut/R_cut/R_cut );
    }
    else if (pos >= 0.0 && pos> R_cut) {
        return R_cut/2;
    }
    else if (pos < 0.0 && pos> -R_cut) {
        return -( abs(pos) - abs(pos*pos*pos)/R_cut/R_cut + abs(pos*pos*pos*pos)/2/R_cut/R_cut/R_cut);
    }
    else if (pos < 0.0 && pos< -R_cut) {
        return -R_cut/2;
    }
    else{
        return 0.0 ;
    }
}


double kernel(arma::vec& dist, arma::vec& Prefactor, const double R_cut, const double flag_kernel) {
	double result;
	if (flag_kernel == 0) {
		result = arma::sum(0.5/(2.0* R_cut)*Prefactor % ( arma::sign(dist + R_cut)- arma::sign(dist - R_cut)));
			//exp(-(dist_sq) / (2 * R_cut*R_cut)) / sqrt(2 * arma::datum::pi *R_cut*R_cut));
	}
	else if (flag_kernel==1) {
//        arma::vec dist_sq = dist % dist;
		result = arma::sum(Prefactor % arma::exp(-(dist % dist) / (2 * R_cut*R_cut)) / sqrt(2 * arma::datum::pi *R_cut*R_cut));
	}
	else if (flag_kernel==2) {
        arma::vec dist_tmp = arma::abs(dist);
        
        arma::vec result_tmp= 1/R_cut- (3/(R_cut*R_cut*R_cut))*( dist_tmp % dist_tmp ) + (2/(R_cut*R_cut*R_cut*R_cut))*( (dist_tmp % dist_tmp) % dist_tmp);
        
        for (int i = 0; i < dist_tmp.n_elem; i++) {
            if (dist_tmp(i) > R_cut){
            result_tmp(i) = 0.0;
            }
        }
        result = arma::sum( Prefactor % result_tmp) ;
	}

	return result;
}

double kernel_derivative(arma::vec& dist, arma::vec& Prefactor, const double R_cut, const double flag_kernel) {
    double result;
    if (flag_kernel == 0) {
        std::cout << "No kernel derivative for constant type";
    }
    else if (flag_kernel==1) {
        arma::vec dist_sq = dist % dist;
        result = arma::sum(Prefactor % (- dist / (R_cut*R_cut) ) % arma::exp(-(dist_sq) / (2 * R_cut*R_cut)) / sqrt(2 * arma::datum::pi *R_cut*R_cut));
    }
    else if (flag_kernel==2) {
        arma::vec dist_tmp = dist;
        arma::vec result_tmp=dist;
        
        for (int i = 0; i < dist_tmp.n_elem; i++) {
            if (dist_tmp(i) > R_cut){
                result_tmp(i) = 0.0;
            }
            if (dist_tmp(i) >= 0.0 && dist_tmp(i)<R_cut) {
                result_tmp(i) = - 3*2*dist_tmp(i)/R_cut/R_cut/R_cut + 2*3*dist_tmp(i)*dist_tmp(i)/R_cut/R_cut/R_cut/R_cut;
            }
            else if (dist_tmp(i) >= 0.0 && dist_tmp(i)> R_cut) {
                result_tmp(i) = 0.0;
            }
            else if (dist_tmp(i) < 0.0 && dist_tmp(i)> -R_cut) {
               result_tmp(i) = - ( 3*2*dist_tmp(i)/R_cut/R_cut/R_cut + 2*3*dist_tmp(i)*dist_tmp(i)/R_cut/R_cut/R_cut/R_cut);
            }
            else if (dist_tmp(i) < 0.0 && dist_tmp(i)< -R_cut) {
                result_tmp(i) = 0.0;
            }
            else{
                result_tmp(i) = 0.0;
            }
        }
        
        result = arma::sum( Prefactor % result_tmp) ;
    }
    
    return result;
}

int compute_density(const int n_steps, const int n_bin, const double Lx, const double Lz, const double R_cut, arma::vec& Mass, arma::vec& density, arma::mat& y, arma::vec& bin_vec, const double flag_kernel) {
	for (int i_bin = 0; i_bin < n_bin; i_bin++) {
		for (int i_step = 0; i_step < n_steps; i_step++) {
			arma::vec dist = bin_vec(i_bin) - y.col(i_step);
			density(i_bin) += (1 / Lx / Lz * kernel ( dist , Mass, R_cut, flag_kernel)) / n_steps;
	}
	}
	return 0;
}

int compute_density_derivative(const int n_steps, const int n_bin, const double Lx, const double Lz, const double R_cut, arma::vec& Mass, arma::vec& density_derivative, arma::mat& y, arma::vec& bin_vec, const double flag_kernel) {
    for (int i_bin = 0; i_bin < n_bin; i_bin++) {
        for (int i_step = 0; i_step < n_steps; i_step++) {
            arma::vec dist = bin_vec(i_bin) - y.col(i_step);
            density_derivative(i_bin) += (1 / Lx / Lz * kernel_derivative ( dist , Mass, R_cut, flag_kernel)) / n_steps;
        }
    }
    return 0;
}

int compute_momentum(const int n_steps, const int n_bin, const double Lx, const double Lz, const double R_cut, arma::mat& Velx, arma::vec& momentum, arma::mat& y, arma::vec& bin_vec, const double flag_kernel) {
	for (int i_bin = 0; i_bin < n_bin; i_bin ++) {
		for (int i_step = 0; i_step < n_steps; i_step++) {
			arma::vec dist = bin_vec(i_bin) - y.col(i_step);
			arma::vec vel = Velx.col(i_step);
			momentum(i_bin) += (1 / Lx / Lz * kernel(dist, vel, R_cut, flag_kernel)) / n_steps;
		}
	}
	return 0;
}

int compute_momentum_derivative(const int n_steps, const int n_bin, const double Lx, const double Lz, const double R_cut, arma::mat& Velx, arma::vec& momentum_derivative, arma::mat& y, arma::vec& bin_vec, const double flag_kernel) {
    for (int i_bin = 0; i_bin < n_bin; i_bin ++) {
        for (int i_step = 0; i_step < n_steps; i_step++) {
            arma::vec dist = bin_vec(i_bin) - y.col(i_step);
            arma::vec vel = Velx.col(i_step);
            momentum_derivative(i_bin) += (1 / Lx / Lz * kernel_derivative(dist, vel, R_cut, flag_kernel)) / n_steps;
        }
    }
    return 0;
}


int compute_temperature(const int n_steps, const int n_bin, const double Lx, const double Lz, const double R_cut, arma::mat& velx, arma::mat& vely, arma::mat& velz, arma::mat& ave_velx, arma::vec& temp, arma::mat& y, arma::vec& bin_vec, const double flag_kernel) {
	for (int i_bin = 0; i_bin < n_bin; i_bin++) {
		for (int i_step = 0; i_step < n_steps; i_step++) {
			arma::vec dist = bin_vec(i_bin) - y.col(i_step);
			arma::vec vel = (velx.col(i_step) - ave_velx(i_bin)) % (velx.col(i_step) - ave_velx(i_bin)) + (vely.col(i_step) % vely.col(i_step)) + (velz.col(i_step) % velz.col(i_step));
			temp(i_bin) += (1 / 3.0 / Lx / Lz * kernel(dist, vel, R_cut, flag_kernel)) / n_steps;
		}
	}
	return 0;
}

int compute_stressK(const int n_steps, const int n_bin, const double Lx, const double Lz, const double R_cut, arma::mat& vel1, arma::mat& vel2, arma::mat& ave_vel1, arma::mat& ave_vel2, arma::vec& stressK, arma::mat& y, arma::vec& bin_vec, const double flag_kernel) {
	for (int i_bin = 0; i_bin < n_bin; i_bin++) {
		for (int i_step = 0; i_step < n_steps; i_step++) {
			arma::vec dist = bin_vec(i_bin) - y.col(i_step);
			arma::vec vel = (vel1.col(i_step) - ave_vel1(i_bin)) % (vel2.col(i_step) - ave_vel2(i_bin));
			stressK(i_bin) -= (1 / Lx / Lz * kernel(dist, vel, R_cut, flag_kernel)) / n_steps;
		}
	}
	return 0;
}

int read_file(const std::string& namefile, const int n_steps, const int n_particles, arma::mat& type , arma::mat& x, arma::mat& y, arma::mat& z, arma::mat& vx, arma::mat& vy, arma::mat& vz) {
	std::string line;
	//ifstream myfile("example.txt");
	std::ifstream myfile(namefile);
	if (myfile.is_open())
	{
		int n_line = 0;
		int n_time = 0;
		while (std::getline(myfile, line) && n_time<n_steps)
		{
			if (n_line < 9) {
				n_line += 1;
			}
			else if (n_line == n_particles + 8) {
				n_time += 1;
				n_line = 0;
			}
			else {

				//split line
				std::vector<std::string> result;
				std::istringstream iss(line);
				for (std::string line; iss >> line;) {
					result.push_back(line);
				}

				type(n_line - 8, n_time) = stod(result[1])-1;
				x(n_line - 8, n_time) =	stod(result[2]);
				y(n_line - 8, n_time) = stod(result[3]);
				z(n_line - 8, n_time) = stod(result[4]);
				vx(n_line - 8, n_time) = stod(result[5]);
				vy(n_line - 8, n_time) = stod(result[6]);
				vz(n_line - 8, n_time) = stod(result[7]);
				n_line += 1;
				//std::cout << line << '\n';
			
			}
		}
		myfile.close();
	}
	else std::cout << "Unable to open file";
	return 0;
}

int read_parameters(const std::string& namefile, int &n_types, int &n_particles, double &Lx, double &Ly, double &Lz) {
	std::string line;
	std::ifstream myfile(namefile);
	if (myfile.is_open())
	{
		int n_line = 0;
		n_particles = 0;
		n_types = 1;
		while (std::getline(myfile, line) && n_line< n_particles+8)
		{
			if (n_line == 3) {
				n_particles = stoi(line);
				n_line += 1;
			}
			else if (n_line == 5) {
				std::vector<std::string> result;
				std::istringstream iss(line);
				for (std::string line; iss >> line;) {
					result.push_back(line);
				}
				Lx = stod(result[1]) - stod(result[0]);
				n_line += 1;
			}
			else if (n_line == 6) {
				std::vector<std::string> result;
				std::istringstream iss(line);
				for (std::string line; iss >> line;) {
					result.push_back(line);
				}
				Ly = stod(result[1]) - stod(result[0]);
				n_line += 1;
			}
			else if (n_line == 7) {
				std::vector<std::string> result;
				std::istringstream iss(line);
				for (std::string line; iss >> line;) {
					result.push_back(line);
				}
				Lz = stod(result[1]) - stod(result[0]);
				n_line += 1;
			}
			else if (n_line > 9) {
				std::vector<std::string> result;
				std::istringstream iss(line);
				for (std::string line; iss >> line;) {
					result.push_back(line);
				}
				if (stod(result[1]) > n_types) {
					n_types = stoi(result[1]);
				}
				n_line += 1;
			}
			else {
				n_line += 1;
			}
		}
		myfile.close();
	}
	else std::cout << "Unable to open file";
	return 0;
	}


int read_parameter_file(const std::string& namefile,  arma::mat& epsilon , double &Force_cut) {
	std::string line;
	std::ifstream myfile(namefile);
	if (myfile.is_open())
	{
		while (std::getline(myfile, line) )
		{
				//split line
				std::vector<std::string> result;
				std::istringstream iss(line);
				for (std::string line; iss >> line;) {
					result.push_back(line);
				}
				
				if (&result[0] != NULL){
					std::string str1= result[0];
					//std::cout << "str1" << str1 << "\n";
						if (str1.compare("pair_coeff")==0){
						std::string stri= result[1];
						std::string strj= result[2];
						if (stri.compare("*")==0 || strj.compare("*")==0){
							Force_cut=	stod(result[5]);
							std::cout << "All epsilon set to 1" << '\n';
						}
						else
						{
							int i=stoi(result[1])-1;
							int j=stoi(result[2])-1;
							double e_value=stod(result[3]);
							epsilon(i,j)=e_value;
							epsilon(j,i)=e_value;
							Force_cut=	stod(result[5]);
						}
					}
				}

		}
		myfile.close();
		std::cout << "Epsilon " << epsilon << "\n";
		std::cout << "Force_cut " << Force_cut << "\n";
	}
	else{
		std::cout << "Unable to open file";
	}
	return 0;
}


void uvec_push(arma::uvec & v, const double value) {
    arma::uvec av(1);
    av.at(0) = value;
    
    v.insert_rows(v.n_rows, av.row(0));
}

int main(int argc, const char **argv) {
	clock_t tStart = clock();
	
	std::string filename;
	std::string parameter_filename;
	int	n_steps;
	int	n_bin;

	int n_types=0;
	int n_particles=0;
	double	Lx=0;
	double	Ly=0;
	double	Lz=0;

	// Read fundamental command line parameters 
		std::vector<std::string> args(argv, argv + argc);
		//check filename
		if (args[1] != "filename") {
			std::cout << "Arg 1 wrong";
			return 0;
		}
		else {
			filename = args[2];
			parameter_filename = args[3]; //lammps input filename
		}
		//check n_steps
		if (args[4] != "n_steps") {
			std::cout << "Arg 4 wrong";
			return 0;
		}
		else {
			n_steps = stoi(args[5]);
		}

		//check bins
		if (args[6] != "n_bins") {
			std::cout << "Arg 6 wrong";
			return 0;
		}
		else {
			n_bin = stoi(args[7]);
		}

	//Read fundamental parameters from dump file
		read_parameters(filename, n_types, n_particles, Lx, Ly, Lz);

	double	Ly_start = 0;//4;
	double	Ly_end = Ly_start+Ly;//10;

	//flags for the kernel
	int flag_kernel = 0;
	int flag_density = 0;
	int flag_momentum = 0; // 1=x, 2=y, 3=z


	arma::Mat<double> type = arma::zeros(n_particles, n_steps);
	arma::Mat<double> x = arma::zeros(n_particles, n_steps);
	arma::Mat<double> y = arma::zeros(n_particles, n_steps);
	arma::Mat<double> z = arma::zeros(n_particles, n_steps);
	arma::Mat<double> vx = arma::zeros(n_particles, n_steps);
	arma::Mat<double> vy = arma::zeros(n_particles, n_steps);
	arma::Mat<double> vz = arma::zeros(n_particles, n_steps);
    
	arma::vec bin_vec = arma::linspace<arma::vec>(Ly_start, Ly_end, n_bin);
	bin_vec.save("bin_vec.txt", arma::arma_ascii);
    
    int pdf_n_bin=30;
    double pdf_bin_vec_start=0;
    double pdf_bin_vec_end=5;
    
    arma::umat hist_vy = arma::zeros<arma::umat>(n_bin,pdf_n_bin);
    arma::Mat<double> pdf_vy = arma::zeros(n_bin,pdf_n_bin);
    
    arma::vec pdf_bin_vec = arma::linspace<arma::vec>(pdf_bin_vec_start, pdf_bin_vec_end, pdf_n_bin);
    double dy=(pdf_bin_vec_end-pdf_bin_vec_start)/(pdf_n_bin-1);
    
    pdf_bin_vec.save("pdf_bin_vec.txt", arma::arma_ascii);
    
//    arma::vec Mass = arma::ones(n_particles);
//    double  R_cut; //characteristic length of the kernel function
//
//    //Paramaters for lennard-jones interaction
//        double Force_cut = 3.5; //cutting radius LJ
//        arma::Mat<double>  epsilon = arma::ones(n_types, n_types);
//
//        read_parameter_file(parameter_filename, epsilon, Force_cut);
//
//        int type_wall = n_types - 1;
//        double Force_cut_sq = Force_cut*Force_cut;
	///////////////////////////////////////////////////////////////

//    for (size_t i_arg = 11; i_arg < args.size(); i_arg++) {
//    }

	//read file
	std::cout << "Reading data\n";
	read_file(filename, n_steps, n_particles, type, x, y, z, vx, vy, vz);
	//std::cout << "x:\n" << x << "\n";

	//Compute Density
//    std::cout << "Computing density \n";
//    arma::vec density = arma::zeros(n_bin);
    
//    compute_density(n_steps, n_bin, Lx, Lz, R_cut, Mass, density, y, bin_vec, flag_kernel);
//    density.save("density.txt", arma::arma_ascii);


	//Compute Velocity
//    arma::vec momentum = arma::zeros(n_bin);
//    arma::vec velocity_x;
//
//    std::cout << "Computing momentum and average velocity in x direction\n";
//    compute_momentum(n_steps, n_bin, Lx, Lz, R_cut, vx, momentum, y, bin_vec, flag_kernel);
//    velocity_x = momentum / density;
//    momentum.save("momentum.txt", arma::arma_ascii);
//    velocity_x.save("velocity_x.txt", arma::arma_ascii);


    std::cout << "Computing velocity distribution in y direction\n";

    for (int i_bin = 0; i_bin < n_bin; i_bin++) {
        arma::vec vy_tmp(1);
        for (int i_step = 0; i_step < n_steps; i_step++) {
            for (int i_particle = 0; i_particle < n_particles; i_particle++) {
                if ( ( y(i_particle,i_step)>bin_vec(i_bin)-dy/2) and (y(i_particle,i_step)<bin_vec(i_bin)+dy/2) ) {
//                   arma::vec a(1);
//                    a(0)= vy(i_particle,i_step);
//                    printf( "a: %.2c\n", typeid(a));
//                    a.print();
//                    uvec_push( vy_tmp, a );
                    int rows= vy_tmp.n_rows;
                    vy_tmp.resize( rows+1);
                    vy_tmp(rows-1)=std::abs(vy(i_particle,i_step));
//                    vy_tmp.insert_cols( rows, a);

//                    printf( "vy: %.2f\n", (double)a);
//                    printf( "NROWS: %.2d\n", (int)vy_tmp.n_rows);
                }
            }
        }
        
//        vy_tmp.print("vy_tmp:");
        int rows= vy_tmp.n_rows;
        vy_tmp.resize( rows-1);
//        vy_tmp.print("vy_tmp:");
        
        hist_vy.row(i_bin) = arma::hist(vy_tmp,pdf_bin_vec).t();
        
        pdf_vy.row(i_bin) = 1/( dy * rows)* arma::conv_to<arma::mat>::from(hist_vy.row(i_bin));
//        hist_vy = arma::sum(hist_vy, 0);
    }
    
    
//    arma::mat Z = arma::trapz(pdf_bin_vec , pdf_vy.row(3).t());
//    Z.print();
    
    pdf_vy.save("pdf_vy.txt", arma::arma_ascii);
    
	printf("Time taken: %.10fs\n", (double)(clock() - tStart)/CLOCKS_PER_SEC);
		//std::cin.get();
	return 0;
}

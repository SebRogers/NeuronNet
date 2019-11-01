#include "random.h"

RandomNumbers::RandomNumbers(unsigned long int s):
seed(s) {
	if (seed==0){ 
		std::random_device rd;
		seed = rd();	
	}
		rng = std::mt19937 (seed);
}

//recoit un vecteur d'une certaine taille et le rempli avec la distribution unif.
void RandomNumbers::uniform_double(std::vector<double>& vec, double lower, double upper){
	 for (size_t i(0) ; i<vec.size(); i++){										
		vec[i]=uniform_double(lower, upper);
	 }
} 
    
//retourne un nombre de distribution unif.
double RandomNumbers::uniform_double(double lower, double upper){

	std::uniform_real_distribution<double> uniform_d(lower, upper);		// dans la bibliothèque random
	return uniform_d(rng);
}
	
//recoit un vecteur d'une certaine taille et le rempli avec la distribution norm.
void RandomNumbers::normal(std::vector<double>& vec, double mean, double sd){
	   for ( size_t i(0) ; i<vec.size(); i++){
			vec[i]=normal(mean, sd);
	   }
}
  
//retourne un nombre de distribution unif.  
double RandomNumbers::normal(double mean, double sd){
	
	std::normal_distribution<double> normal_d(mean, sd); 			// dans la bibliothèque random
	return normal_d(rng);
}
 
//recoit un vecteur d'une certaine taille et le rempli avec la distribution poiss.   
void RandomNumbers::poisson(std::vector<int>& vec, double mean){
	for ( size_t i(0) ; i<vec.size(); i++){
			vec[i]=poisson(mean);
	   }
}

//retourne un nombre de distribution poisson.  
int RandomNumbers::poisson(double mean){
	std::poisson_distribution<> poiss(mean);
	return poiss(rng);
}


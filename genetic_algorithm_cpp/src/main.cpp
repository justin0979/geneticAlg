// Program: Genetic Algorithm
// Author: Justin Mangawang
// Date: February 2017 

#include "genetic_algorithm.cpp"
#include <time.h>

int main(int argc, char** argv)
{
	GA g;
	srand(time(NULL));
  std::cout << std::fixed;
  std::cout << std::setprecision(1);
	
	g.geneticAlgorithm();
	
	return 0;
}


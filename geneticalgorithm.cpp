// Program: Genetic Algorithm
// Author: Justin Mangawang
// Date: February 2017 

#include <iomanip>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <time.h>
#include <typeinfo>

class GA
{
public:
	GA();
	void geneticAlgorithm();
	// Function to randomly generate a population of 20 chromosomes.  Based on the Pco, a portion
	// of the population will be replicated into the nextGenChromosomes array with the remainder
	// of the population being crossed-over and the children stored in the same array.
	// The function will generate chromosome 1010101010, record the iterations it took, re-run
	// 19 more times for each Pco and finally calulate the average number of iterations it took
	// to find the target chromosome.  The function is run 100 times (20 runs per Pco).
	// Postcondition:  Average number of iteraions for each Pco is calculated.
	void printInitialPopulation(int arr[20][11]);  
	// Function to printout the specified population.
	// Postcondition:  Prints population to screen.
	int getTargetCounter() {return targetCounter;}; 
	// Function to return targetCounter, which ensures that 100 trials were run (20 runs per pco).
	// Postcondition:  Returns targetCounter.
	int getIteration(){return targetIteration;};
	// Function to return targetIteration, which records the number of generations it took to find 
	// the target chromosome.
	// Postcondition:  Returns targetIteration.
	
	// VARIABLES
	static const int popSize = 20;
	static const int chromSize = 11; // 0 will hold fitness value, 1 to 10 will hold chromosome
									 // data.
private:
	void gaWorkHorse(const int&, const int&, const float);
	// Function to conduct the main work of the genetic algorithm.
	void setInitialPop();
	// Function to randomly set the initial population for the genetic algorithm.
	void filler();
	// Function NOT USED, but if too few iterations occur, the remaining populations in fourGen are
	// filled with -1.
	void swap();
	// Function to swap the contents of the last two populations of fourGen, in the event that the
	// last generation is even.
	// Postcondition:  fourGen[2][20][11] = fourGen[3][20][11] && 
	// fourGen[3][20][11] = fourGen[2][2 0][11]
	void fillFourGen(const int& i);
	// Function to store the array containing the first two generations after initial population
	// and the last two generated generations.
	void printFirstTwoAndLastTwoGenerations();
	// Function to output the array containing the first two generations after the initial 
	// population.
	// and the last two generated generations.
	// Postcondition:  prints crossed-over and mutated populations.
	void printPopulationsOrderLabel(const int&);
	// Function to assist printFirstTwoAndLastTwoGenerations() by printing label for 2nd, 3rd, next to last and last 
	// populations.
	void printFirstTwoAndLastTwoGenerationsDecoration(const int&);
	// Funtion prints out ".)" + data for the 2nd, 3rd, next to last and last populations.
	void setPop(int arr[20][11]); 
	// Function to randomly generate the 20 chromosome population.
	// Postcondition:  pop array will have 20 randomly generated chromosomes.
	void selectAndStoreReplication(const float&); 
	// Function to randomly select 1 - Pco chromosomes and store those  in the nextGenChromosomes 
	// array.
	// Postcondition:  nextGenChromosomes will contain 1 - Pco replicated chromosomes.
	void setFitness(int arr[20][11], const int&, const int&);  
	// Function to compare each individual chromosome with the target chromosome and increment the 
	// fitness value of each gene match between the target and analyzed chromosomes.
	// Postcondition:  Index 0 of population arrays (pop, nextGenChromosomes) will be updated with 
	// chromosome fitness value.
	void fillCopyChromosomeList(int original[20][11], int copy[20][11]);
	// Function to copy a population array into another population array.
	void fillChromosomeRouletteVector();  
	// Function to fill chromosomeRoulette array based on the fitness value of each chromosome.  
	// i.e., Chromosome in index 1 has fitness value of 5, chromosomeRoulette's indices 0-4 will 
	// store the int value of 1. Chromosome in index 2 has fitness value of 2, 
	// chromosomeRoulette's indices 5-6 will store the int value 2; thus, the 
	// chromosomeRoulette contents from index 0 - 6 will be:  1 1 1 1 1 2 2.
	// Postcondition:  chromosomeRoulette vector is filled according to chromosome fitness values.
	void recombination(const float&); 
	// Function to randomly select chromosome pairs for cross-over.  
	void crossOver(const int&);
	// Function to conduct cross over of first half of first chromosome with back half of 
	// second chromosome and vice versa.  Children chromosomes are then stored in 
	// nextGenChromosomes array.
	// Precondition:  2 parent chromosomes from "original" generation.
	// Postcondition:  2 children chromosomes stored in nextGenChromosomes.
	void mutation();
	// Function to randomly select one chromosome from the newly generated nextGenChromosomes 
	// array and within that chromosome, a gene is randomly selected to change from a 0 to a 1 
	// or vice versa.
	// Precondition:  Newly crossed over nextGenChromosomes.
	// Postconditoin:  One chromosome with one "new" gene.
	void printFitness();
	// Function to output to screen the fitness values of every chromosome.
	void targetMatch(int arr[20][11]); 
	// Function to compare a chromosome population to the target chromosome.
	// Precondtion:  match = 0;
	// Postcondition:  if match found, return 1, else 0.
	void printPcoAvg(float*);
	// Function to calculate and print out to screen the average Pco after 20 runs.
  void printOnlyPcoSeventyOrZeroPercentLabel(float pco);
  void printChromosomeLineNumber(int line);
  void printChromosomeFitnessAndGenes(const int &i, const int &j);

	// VARIABLES
	int pop[popSize][chromSize];	// index 0 holds fitness value, indices 1 to 10 hold 0's or 
									// 1's.  Array holds
									// initial chromosome population throughout duration of program
	float popPercentage[popSize];	// CAN DELETE.
	std::vector<int> chromosomeRoulette; // holds number of chromo baseed on its fitness value.
	int nextGenChromosomes[popSize][chromSize]; // holds the recombinant and replicated chromosomes
	int copyChromosomeList[popSize][chromSize]; // holds copy of iterated chromosome pop; so, 
							// initial pop won't be erased and so every recombined pop will be 
							// stored
	int* compare; // array to ensure no replicated or co'd chrom's are duplicated.
	int end;
	bool match; // value to check if a chromosome matches the target 1010101010.
	int targetCounter; // used to count the number of 1010101010's occur.
	int fitness;
//	int totFitness; // didn't use.
	float chromGetsPicked;
	int iteration; // count number of iterations from 1 to target chromosome.
	int targetIteration; // records number of iterations it took to get target. 
	int storeIteration[5][20]; // Stores the Pco value and that Pco value's number of iterations 
							   // for 20 runs.
	int fourGen[4][20][11]; // stores crossed over and mutated pop's for 1st two and last 2 
							// iteratations.
};

/*************************MAIN FUNCTION******************************/

int main(int argc, char** argv)
{
	GA g;
	srand(time(NULL));
	
	g.geneticAlgorithm();
	
	return 0;
}

/***************************CONSTRUCTOR*******************************/
GA::GA()
{
	match = 0;
	end = 0;
	targetCounter = 0;
	iteration = 1;
	targetIteration = 0;
	fitness = 0;
	chromGetsPicked = 0.0;
	
	for(int i = 0; i < popSize; i++)
	{
		popPercentage[i] = 0.0;
		for(int j= 0; j < chromSize; j++)
		{
			pop[i][j] = 0;
			nextGenChromosomes[i][j] = 0;
			copyChromosomeList[i][j] = 0;
		}
	}
}

/***********************geneticAlgorithm()**********************
	Function to find 1010101010 with Pco.  Once match found, iterations will stop.
*/
void GA::geneticAlgorithm()
{
  const int numOfPcos = 5;
	float Pco[numOfPcos] = {0.7, 0.3, 0.5, 0.9, 0.0};
	setInitialPop();
	targetMatch(pop);
	for(int pcoIndex = 0; pcoIndex < numOfPcos; pcoIndex++)
	{
    printOnlyPcoSeventyOrZeroPercentLabel(Pco[pcoIndex]);

		for(int i = 0; i < 20; i++)
		{
			gaWorkHorse(pcoIndex, i, Pco[pcoIndex]);
			targetIteration = iteration;
			targetCounter++;
			iteration = 1;
			match = 0;
			storeIteration[pcoIndex][i] = getIteration() + 1;
			fillCopyChromosomeList(pop, copyChromosomeList);
		}

		if(Pco[pcoIndex]*10 == 7 || Pco[pcoIndex]*10 == 0)
		{
			if(end % 2 == 0 && end > 5)
			{
				swap();
			}
			printFirstTwoAndLastTwoGenerations();
		}
	}
	printPcoAvg(Pco);
}

void GA::printOnlyPcoSeventyOrZeroPercentLabel(float pco) {
  if(pco*10 == 7 || pco*10 == 0) {
			std::cout << "\nPco = " << pco << "\t" <<  "F(i) of 20 chromosomes"<< std::endl;
  }
  return;
}

/**************************setInitialPop()*************************
	Function to randomly select each "gene" in each "chromosome in the population.
*/
void GA::setInitialPop() {
	std::cout << "Initial Population.\n" << std::endl;
	std::cout << " i\tF(i)\tChromosome data" << "\n----\t----\t-------------------" << std::endl;
	setPop(pop); // stores original generated pop and saves it in pop array.
	printInitialPopulation(pop); // display chromosomes. 

}

/**************************gaWorkHorse()**************************
	Function to cycle through the template of the genetic algorithm.
*/
void GA::gaWorkHorse(const int& j, const int& i, const float Pco) {
	int counter = 0; // track which index stores the 1st two and last 2 generations.

	while(!match) // loop is the "work-horse" of the genetic algorithm.
	{	
		iteration++;
		fillChromosomeRouletteVector(); // creates percentage of any distinct chromosome to
										// be picked in replication/co.
		selectAndStoreReplication(Pco);

		recombination(Pco); // stores the recombinant chromosomes.
		targetMatch(nextGenChromosomes); // Checks if the target chromosome is generated.
		chromosomeRoulette.clear(); // empties vector for next recombinant iteration, so 
									// new data is NOT added to old.
		if(i == 0) // condition verifies stored data and printed fitness values come from 
				   // the 1st of 20 runs only.
		{
			if(counter < 2)
			{
				fillFourGen(counter);
			}
			if(counter > 1)
			{
				if(counter % 2 == 0)
				{
					fillFourGen(2);
				}
				else
				{
					fillFourGen(3);
				}
			}
			if(j==0 || j==4) // condition ensures that only one run will execute for 
							 // both Pco 0.7 and 0.0.
			{
				printFitness();
			}
		}	
		counter++;
	}
}

/***************************setPop()*******************************
	Function to randomly generate 20 "new chromosomes".
*/
void GA::setPop(int arr[popSize][chromSize])
{
	int randNum;

	for(int i = 0; i < popSize; i++)
	{
		for(int j = 1; j < chromSize; j++)
		{
			randNum = rand() % 2;
			arr[i][j] = randNum;
			setFitness(arr, i, j); // line 111
		}
		fitness = 0;
	}
	fillCopyChromosomeList(arr, copyChromosomeList);
}

/***************************printInitialPopulation()******************************
	Function to print out a population.
*/
void GA::printInitialPopulation(int arr[popSize][chromSize])
{
	for(int i = 0; i < popSize; i++)
	{
		std::cout << i + 1 << ".)\t";
		
		for(int j = 0; j < chromSize; j++)
		{
			std::cout << arr[i][j] << " ";
			if(j == 0)
			{
				std::cout << "\t";
			}
		}
		std::cout << std::endl;
	}
}

/***************************setFitness()******************************
	Function compares value in each index of pop array "chromosomes"
	to target "chromosome".  Where there is a match, fitness increments.
*/
void GA::setFitness(int arr[popSize][chromSize], const int& n, const int& index)
{
	if(index%2 == 0 && arr[n][index] == 0) // Checks if even index contains 0.
	{
		fitness++;
	}
	if(index%2 != 0 && arr[n][index] == 1) // Checks if odd index contains 1.
	{
		fitness++;
	}
	arr[n][0] = fitness; 
}

/**************************fillChromosomeRouletteVector()*************************
	Function to create a "roulette wheel" based on the fitness value of "chromosome" specified by
	pop[i][0].  The fitness value stored in pop[i][0] will determine the number of times that the 
	index holding the chromosome will be put into the vector.  
*/
void GA::fillChromosomeRouletteVector()

{
	for(int i = 0; i < popSize; i++)
	{
		for(int j = 0; j < copyChromosomeList[i][0]; j++)
		{
			chromosomeRoulette.push_back(i);
		}
		if(copyChromosomeList[i][0] == 0) // if fitness is 0, then chromosome will have only one 
										  // entry in Roulette.
		{
			chromosomeRoulette.push_back(i);
		}
	}
}

/*************************selectAndStoreReplication()*******************************
	Function to randomly choose replication chromosomes from population.  
	Function will check to ensure dulicate replication chromosomes are not allowed.
*/
void GA::selectAndStoreReplication(const float& Pco)
{
	int replicationNum = popSize * (1 - Pco);
	int random = 0;
	compare = new int[replicationNum]; // will store the generated random numbers to prevent 
									   // replicate duplications.
					                   // its sole urpose is for duplicate prevention.
	for(int i = 0; i < replicationNum; i++)
	{
		random = rand() % chromosomeRoulette.size();
		here:
		for(int j = 0; j < i; j++) // checks array for duplicate values.
		{
			if(random == compare[j])	
			{
				random = rand() % chromosomeRoulette.size();	// if duplicate value found, new 
																// number is generated.
				goto here;	// after new random generated, go back thru array to check again for 
							// duplicates.
			}
		}
		compare[i] = chromosomeRoulette[random];

		for(int k = 1; k < chromSize; k++) // STORES the replicated chromosomes
		{
			nextGenChromosomes[i][k] = copyChromosomeList[chromosomeRoulette[random]][k];
			setFitness(nextGenChromosomes, i, k);
		}
		fitness = 0;
	}
	delete compare;
}

/******************************fillCopyChromosomeList()************************************
	Function to make copy of a 2D array.
*/
void GA::fillCopyChromosomeList(int original[popSize][chromSize], int copy[popSize][chromSize])
{
	for(int i = 0; i < popSize; i++)
	{
		for(int j = 0; j < chromSize; j++)
		{
			copy[i][j] = original[i][j];
		}
	}
}

/*******************recombination****************************************8
	Function to prepare for crossover.  Chromosomes will be retrieved from copyChromosomeList 
	array.
*/
void GA::recombination(const float &Pco)
{
	int randomPosition = 0;
	int pairs = 0;
	int counter = 0;
	int coNum = Pco * popSize; // holds calculated number of total pop based on Pco.
	int temp = popSize - coNum; // holds index value in nextGenChromosomes, so new co chrom's will 
								// be after replicants. 
	compare = new int[coNum];
	
	for(int i = 0; i < coNum/2; i++) // coNum/2 b/c nested while loop iterates 2 times for 2 
									 // parents.
	{
		while(pairs++ < 2) // if coNum was not dived by 2, then there would be 2 * coNum parents 
						   // and children.
		{
			randomPosition = rand() % chromosomeRoulette.size();
			here:
			for(int x = 0; x < counter; x++) // for loop to prevent same chromosome crossing over 
											 // with itself.
			{
				if(compare[x] == chromosomeRoulette[randomPosition])
				{
					randomPosition = rand() % chromosomeRoulette.size();
					goto here;
				}
			}
			compare[counter++] = randomPosition;
		}
		pairs = 0;
	}
	crossOver(coNum);
	mutation();
	if(Pco == 0.7)
	{
		printFitness();
	}
	fillCopyChromosomeList(nextGenChromosomes, copyChromosomeList);
	delete compare;
}
/*********************crossOver()***********************************************************
	Function to conduct the actual crossover between two selected chromosomes.
*/
void GA::crossOver(const int& num)
{
	int index = 0;
	for(int i = (popSize - num); i < popSize; i+=2) // increment by 2 b/c two parents produce two 
													// children
	{
		for(int j = 1; j < 6; j++) // front half of each crossover pair
		{
			nextGenChromosomes[i][j] = copyChromosomeList[chromosomeRoulette[compare[index]]][j];
			nextGenChromosomes[i+1][j] = copyChromosomeList[chromosomeRoulette[compare[index+1]]][j];
		}
		for(int j = 6; j < chromSize; j++) // back half of each crossover pair
		{
			nextGenChromosomes[i][j] = copyChromosomeList[chromosomeRoulette[compare[index+1]]][j];
			nextGenChromosomes[i+1][j] = copyChromosomeList[chromosomeRoulette[compare[index]]][j];
		}
		index+=2;
	}
}

/*********************************mutation()**************************************************
	Function to randomly select one gene from one chromosome in the population to swap a 1 to 0 
	or a 0 to 1.  Afterwards, a new fitness level is calculated for the entire pop.
*/
void GA::mutation()
{
	int chromToMutate = rand() % popSize;
	int geneToMutate = (rand() % (chromSize - 1) + 1);
	if(nextGenChromosomes[chromToMutate][geneToMutate] == 1)
	{
		nextGenChromosomes[chromToMutate][geneToMutate] = 0;
	}
	else
	{
		nextGenChromosomes[chromToMutate][geneToMutate] = 1;
	}

	for(int i = 0; i < popSize; i++) // for loop to set the fitness of recombinant chromosomes.
	{
		for(int j = 1; j < chromSize; j++)
		{
			setFitness(nextGenChromosomes, i, j); // line 172
		}
		fitness = 0;
	}
}

/**********************************printFitness()*******************************************
	Function to print fitness levels of each chromosome to the screen.
*/
void GA::printFitness()
{
	end = iteration;
	std::cout << "iteration " << iteration << "\t";
	for(int i = 0; i < popSize; i++)
	{
		std::cout << nextGenChromosomes[i][0] << "  ";
	}
	std::cout << std::endl;
}

/**************************targetMatch()**********************************8
	Function to check if a fitness of 10 is found.
*/
void GA::targetMatch(int arr[popSize][chromSize])
{
	for(int i = 0; i < popSize; i++)
	{
		if(arr[i][0] == 10)
		{
			match = true;
		}
	}
}

/****************************printPcoAvg()*********************************
	Function to print out the average number of generations for specified Pco
*/

void GA::printPcoAvg(float* Pco)
{
	int avg = 0;
	int totalIterations = 0;
	for(int i = 0; i < 5; i++)
	{
		for(int j = 0; j < 20; j++)
		{
			totalIterations += storeIteration[i][j];
		}
		avg = totalIterations/20;

		std::cout << "Average number of generations for Pco " << Pco[i] 
				  << " is " << avg << "." << std::endl;
	}
}

/******************************************fillFourGen()*******************************************
	Function to store the first two generations after initial pop and last two generated 
	generations.
*/
void GA::fillFourGen(const int& i)
{
	for(int j = 0; j < popSize; j++)
	{
		for(int k = 0; k < chromSize; k++)
		{
			fourGen[i][j][k] = nextGenChromosomes[j][k];
		}
	}
}

/*********************************printFirstTwoAndLastTwoGenerations()*************************************************
	Function to output the first two generations after initial pop and the last two generations 
	generated.
*/
void GA::printFirstTwoAndLastTwoGenerations()
{
	std::cout << std::endl;
	if(end < 4)
	{
		for(int i = 0; i < (end - 1); i++)
		{
			printPopulationsOrderLabel(i);
			printFirstTwoAndLastTwoGenerationsDecoration(i);
		}
	}	
	else
	{
		for(int i = 0; i < 4; i++)
		{
			printPopulationsOrderLabel(i);
			printFirstTwoAndLastTwoGenerationsDecoration(i);
		}
	}
}

/**********************************************printPopulationsOrderLabel()*****************************
	Fuction printouts out the order of 2nd, 3rd and next to last and last pop listing.
*/
void GA::printPopulationsOrderLabel(const int& i) {
		
	if(i < 2)
	{
		if(i == 0)
		{
			std::cout << "The first generation after initial population:" << std::endl;
		}
		else
		{
			std::cout << "The second generation after the initial population:" << std::endl;
		}
	}
	else
	{
		if(i == 2)
		{
			std::cout << "The second to last generation:" << std::endl;
		}
		else
		{
			std::cout << "The last generation:" << std::endl;
		}
	}	
}

/**********************************************printFirstTwoAndLastTwoGenerationsDecoration()*****************************
	Function adds ".)" + data to printout.
*/
void GA::printFirstTwoAndLastTwoGenerationsDecoration(const int& i) {
	for(int j = 0; j < popSize; j++)
	{
    printChromosomeLineNumber(j + 1);
    printChromosomeFitnessAndGenes(i,j);

		std::cout << std::endl;
	}
	std::cout << std::endl;
}

void GA::printChromosomeLineNumber(int lineNumber) {
  std::cout << lineNumber << ".)\t";
}

void GA::printChromosomeFitnessAndGenes(const int &i, const int &j) {
  for(int k = 0; k < chromSize; k++)
  {
		std::cout << fourGen[i][j][k] << " ";
		if(k == 0)
		{
			std::cout << '\t';	
		}
  }
}

/********************************************************swap()*****************************
	Function to swap the conents of the last two indices of fourGen array.
*/
void GA::swap()
{
	int temp[popSize][chromSize];
	for(int i = 0; i < popSize; i++)
	{
		for(int j = 0; j < chromSize; j++)
		{
			temp[i][j] = fourGen[3][i][j];
			fourGen[3][i][j] = fourGen[2][i][j];
			fourGen[2][i][j] = temp[i][j];
		}
	}
}

#ifndef GENETIC_ALGORITHM
#define GENETIC_ALGORITHM
#include <iomanip>
#include <vector>

class GA
{
public:
	GA();
	void geneticAlgorithm();
	void printPopulation(int arr[20][11]);  
	// Function to return targetIteration, which records the number of generations it took to find 
	// the target chromosome.
	// Postcondition:  Returns targetIteration.
	int getIteration(){return targetIteration;};
	
	// VARIABLES
	static const int popSize = 20;
	static const int chromSize = 11; // 0 will hold fitness value, 1 to 10 will hold chromosome
									 // data.
private:
	void gaWorkHorse(const int&, const int&, const float);
	// Function to conduct the main work of the genetic algorithm.
	void fillInitialPop();
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
	void setFitness(int arr[20][11], const int&);  
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
	void chromosomeFitnessTargetMatch(int arr[20][11]); 
	// Function to compare a chromosome population to the target chromosome.
	// Precondtion:  match = 0;
	// Postcondition:  if match found, return 1, else 0.
	void printPcoAvg(float*);
	// Function to calculate and print out to screen the average Pco after 20 runs.
  void printOnlyPcoSeventyOrZeroPercentLabel(float pco);
  void printChromosomeLineNumber(int line);
  void printChromosomeFitnessAndGenes(int *arr);

	// VARIABLES
	int pop[popSize][chromSize];	// index 0 holds fitness value, indices 1 to 10 hold 0's or 
									// 1's.  Array holds initial chromosome population throughout duration of program
	float popPercentage[popSize];	// CAN DELETE.
	std::vector<int> chromosomeRoulette; // holds number of chromo baseed on its fitness value.
	int nextGenChromosomes[popSize][chromSize]; // holds the recombinant and replicated chromosomes
	int copyChromosomeList[popSize][chromSize]; // holds copy of iterated chromosome pop; so, 
							// initial pop won't be erased and so every recombined pop will be 
							// stored
	int* compare; // array to ensure no replicated or co'd chrom's are duplicated.
	int lastIteration;
	bool match; // value to check if a chromosome matches the target 1010101010.
	int targetCounter; // used to count the number of 1010101010's occur.
  const int FITNESS_VALUE_INDEX = 0;
	float chromGetsPicked;
	int iteration; // count number of iterations from 1 to target chromosome.
	int targetIteration; // records number of iterations it took to get target. 
	int storeIteration[5][20]; // Stores the Pco value and that Pco value's number of iterations 
							   // for 20 runs.
	int fourGen[4][20][11]; // stores crossed over and mutated pop's for 1st two and last 2 
							// iteratations.
};

#endif

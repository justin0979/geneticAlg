#include "genetic_algorithm.h"
#include <iostream>
#include <cstdlib>

GA::GA()
{
	match = 0;
	lastIteration = 0;
	targetCounter = 0;
	iteration = 1;
	targetIteration = 0;
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
  const int numOfPcos {5};
	float Pco[numOfPcos] = {0.7, 0.3, 0.5, 0.9, 0.0};

	fillInitialPop();
	chromosomeFitnessTargetMatch(pop);

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
			if(lastIteration % 2 == 0 && lastIteration > 5)
			{
				swap();
			}
			printFirstTwoAndLastTwoGenerations();
		}
	}
	printPcoAvg(Pco);
}

/**************************fillInitialPop()*************************
	Function to randomly select each "gene" in each "chromosome in the population.
*/
void GA::fillInitialPop() 
{
	std::cout << "Initial Population.\n" << std::endl;
	std::cout << " i\tF(i)\tChromosome data (genes)" << "\n----\t----\t-------------------" << std::endl;
	setPop(pop); 
	printPopulation(pop);

}

/***************************setPop()*******************************
	Function to randomly generate 20 "new chromosomes" with genes containing 0 or 1.
*/
void GA::setPop(int chromosomesArray[popSize][chromSize])
{
	int zeroOrOne {0};

	for(int chromosome = 0; chromosome < popSize; chromosome++)
	{
		for(int gene = 1; gene < chromSize; gene++)
		{
			zeroOrOne = rand() % 2;
			chromosomesArray[chromosome][gene] = zeroOrOne;
		}

		setFitness(chromosomesArray, chromosome);
	}
	fillCopyChromosomeList(chromosomesArray, copyChromosomeList);
}

/***************************setFitness()******************************
	Function compares value in each index of pop array "chromosomes"
	to target "chromosome".  Where there is a match => fitness increments.
*/
void GA::setFitness(int chromosomesArray[popSize][chromSize], const int& chromosome)
{
  int fitness {0};

  for(int gene = 0; gene < chromSize; gene++)
  {

	  if(gene%2 == 0 && chromosomesArray[chromosome][gene] == 0) // Checks if even gene contains 0.
	  {
	  	fitness++;
	  }
	  if(gene%2 != 0 && chromosomesArray[chromosome][gene] == 1) // Checks if odd gene contains 1.
	  {
	  	fitness++;
	  }
	  chromosomesArray[chromosome][FITNESS_VALUE_INDEX] = fitness; 
  }
}

/***************************printPopulation()******************************
	Function to print out a population.
*/
void GA::printPopulation(int arr[popSize][chromSize])
{
	for(int i = 0; i < popSize; i++)
	{
    printChromosomeLineNumber(i + 1);
		printChromosomeFitnessAndGenes(arr[i]);
		std::cout << std::endl;
	}
  std::cout << std::endl;
}

void GA::printChromosomeLineNumber(int lineNumber)
{
  std::cout << lineNumber << ".)\t";
}

void GA::printChromosomeFitnessAndGenes(int *chromosome) {
	for(int gene = 0; gene < chromSize; gene++)
	{
		std::cout << chromosome[gene] << " ";
		if(gene == 0)
		{
			std::cout << "\t";
		}
	}
}

void GA::chromosomeFitnessTargetMatch(int chromosomesArray[popSize][chromSize])
{
	for(int chromosome = 0; chromosome < popSize; chromosome++)
	{
		if(chromosomesArray[chromosome][FITNESS_VALUE_INDEX] == 10)
		{
			match = true;
		}
	}
}

void GA::printOnlyPcoSeventyOrZeroPercentLabel(float pco) 
{
  if(pco*10 == 7 || pco*10 == 0) {
			std::cout << "\nPco = " << pco << "\t" <<  "F(i) of 20 chromosomes"<< std::endl;
  }
  return;
}

/**************************gaWorkHorse()**************************
	Function to cycle through the template of the genetic algorithm.

  TODO: check validity of interationOf20Runs name
*/
void GA::gaWorkHorse(const int& pcoIndex, const int& iterationOf20Runs, const float pcoValue) 
{
	int generationCounter = 0; 

	while(!match) 
	{	
		iteration++;
		fillChromosomeRouletteVector(); // creates percentage of any distinct chromosome to
										// be picked in replication/co.
		selectAndStoreReplication(pcoValue);

		recombination(pcoValue); // stores the recombinant chromosomes.
		chromosomeFitnessTargetMatch(nextGenChromosomes); // Checks if the target chromosome is generated.
		chromosomeRoulette.clear(); // empties vector for next recombinant iteration, so 
									// new data is NOT added to old.
		if(iterationOf20Runs == 0) // condition verifies stored data and printed fitness
                               // values come from the 1st of 20 runs only.
		{
			if(generationCounter < 2)
			{
				fillFourGen(generationCounter);
			}
			if(generationCounter > 1)
			{
				if(generationCounter % 2 == 0)
				{
					fillFourGen(2);
				}
				else
				{
					fillFourGen(3);
				}
			}
			if(pcoIndex == 0 || pcoIndex == 4) // condition ensures that only one run will execute for 
							 // both Pco 0.7 and 0.0.
			{
				printFitness();
			}
		}	
		generationCounter++;
	}
}

/**************************fillChromosomeRouletteVector()*************************
	Function to create a "roulette wheel" based on the fitness value of "chromosome".  The fitness 
  value stored in pop[i][0] will determine the number of times that the chromosome will be put 
  into the vector. Chromosomes with fitness value of 0 will be put into the vector one time.
*/
void GA::fillChromosomeRouletteVector()
{
	for(int chromosome = 0; chromosome < popSize; chromosome++)
	{
		for(int vectorEntryTimes = 0; vectorEntryTimes < copyChromosomeList[chromosome][FITNESS_VALUE_INDEX]; vectorEntryTimes++)
		{
			chromosomeRoulette.push_back(chromosome);
		}
		if(copyChromosomeList[chromosome][FITNESS_VALUE_INDEX] == 0) 
		{
			chromosomeRoulette.push_back(chromosome); 
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
	int random {0};
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
		}
		setFitness(nextGenChromosomes, i);
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
	int randomPosition {0};
	int pairs {0};
	int counter {0};
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
	int index {0};
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
		setFitness(nextGenChromosomes, i); 
	}
}

/**********************************printFitness()*******************************************
	Function to print fitness levels of each chromosome to the screen.
*/
void GA::printFitness()
{
	lastIteration = iteration;
	std::cout << "iteration " << iteration << "\t";
	for(int i = 0; i < popSize; i++)
	{
		std::cout << nextGenChromosomes[i][0] << "  ";
	}
	std::cout << std::endl;
}

/****************************printPcoAvg()*********************************
	Function to print out the average number of generations for specified Pco
*/
void GA::printPcoAvg(float* Pco)
{
	int avg {0};
	int totalIterations {0};

	for(int pco = 0; pco < 5; pco++)
	{
		for(int j = 0; j < 20; j++)
		{
			totalIterations += storeIteration[pco][j];
		}
		avg = totalIterations/20;

		std::cout << "Average number of generations for Pco " << Pco[pco] 
				  << " is " << avg << "." << std::endl;
	}
}

/******************************************fillFourGen()*******************************************
	Function to store the first two generations after initial pop and last two generated 
	generations.
*/
void GA::fillFourGen(const int& generation)
{
	for(int chromosome = 0; chromosome < popSize; chromosome++)
	{
		for(int gene = 0; gene < chromSize; gene++)
		{
			fourGen[generation][chromosome][gene] = nextGenChromosomes[chromosome][gene];
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
	if(lastIteration < 4)
	{
		for(int generation = 0; generation < (lastIteration - 1); generation++)
		{
			printPopulationsOrderLabel(generation);
			printFirstTwoAndLastTwoGenerationsDecoration(generation);
		}
	}	
	else
	{
		for(int generation = 0; generation < 4; generation++)
		{
			printPopulationsOrderLabel(generation);
			printFirstTwoAndLastTwoGenerationsDecoration(generation);
		}
	}
}

/**********************************************printPopulationsOrderLabel()*****************************
	Fuction printouts out the order of 2nd, 3rd and next to last and last pop listing.
*/
void GA::printPopulationsOrderLabel(const int& generation) 
{
		
	if(generation < 2)
	{
		if(generation == 0)
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
		if(generation == 2)
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
void GA::printFirstTwoAndLastTwoGenerationsDecoration(const int& generation) 
{
	for(int chromosome = 0; chromosome < popSize; chromosome++)
	{
    printChromosomeLineNumber(chromosome + 1);
    printChromosomeFitnessAndGenes(fourGen[generation][chromosome]);

		std::cout << std::endl;
	}
	std::cout << std::endl;
}

/********************************************************swap()*****************************
	Function to swap the conents of the last two indices of fourGen array.
*/
void GA::swap()
{
  const int THIRD_GEN {3};
  const int SECOND_GEN {2};
	int temp[popSize][chromSize];

	for(int chromosome = 0; chromosome < popSize; chromosome++)
	{
		for(int gene = 0; gene < chromSize; gene++)
		{
			temp[chromosome][gene] = fourGen[THIRD_GEN][chromosome][gene];
			fourGen[THIRD_GEN][chromosome][gene] = fourGen[SECOND_GEN][chromosome][gene];
			fourGen[SECOND_GEN][chromosome][gene] = temp[chromosome][gene];
		}
	}
}


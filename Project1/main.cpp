#include <iostream>
#include <chrono>
#include "Simulator.h"

using namespace std;

int main()
{
	vector<double> sampleLandscape = { 0.0, 0.011, 1.396, 0.82, 0.487, 0.963, 1.487, 1.471, 1.006, 1.234, 1.429, 1.0, 1.011, 1.077, 1.51, 1.393 };
	vector<int64_t> samplePopulation = { 1000000000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	Simulator sim;

	vector<string> genotypes = sim.getGenotypeStrings();

	for (int i = 0; i < samplePopulation.size(); i++)
		cout << genotypes[i] << ' ' << samplePopulation[i] << endl;
	cout << "----------------" << endl;

	auto t1 = chrono::high_resolution_clock::now();
	vector<int64_t> finalPop = sim.simpleSimulation(sampleLandscape);
	auto t2 = chrono::high_resolution_clock::now();

	for (int i = 0; i < finalPop.size(); i++)
		cout << genotypes[i] << ' ' << finalPop[i] << endl;
	cout << "Simulation took " << chrono::duration_cast<chrono::milliseconds>(t2 - t1).count() << " ms" << endl;

	cout << sim.saveTrace("../data/cpp-out.csv") << endl;

	return 0;
}

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <chrono>
#include <random>
#include <cstdint>
#include <fstream>
using namespace std;

class Simulation
{
private:
	int m_loci;
	int m_genotypes;
	double m_probMut;
	int64_t m_carrCap;
	vector<vector<int64_t>> trace;

	default_random_engine* m_generator;
	uniform_real_distribution<double>* m_realDist;
	uniform_int_distribution<unsigned>* m_lociDist;

	void m_growth(const vector<double> &r, vector<int64_t> &N);
	void m_mutation(const vector<double> &r, vector<int64_t> &N);
	void m_death(const vector<double> &r, vector<int64_t> &N);
	void m_stepForward(const vector<double> &r, vector<int64_t> &N);
public:
	Simulation(int numLoci = 4, double probMut = 1.0e-8, int64_t carrCap = 1000000000);
	~Simulation();
	double drand();
	int lrand();
	int poisson(double mu);
	vector<string> getGenotypeStrings();
	void simpleSimulation(const vector<double> &r, vector<int64_t> &N, const int &t = 1200);
	bool saveTrace(const string &filename);
};

Simulation::Simulation(int numLoci, double probMut, int64_t carrCap)
	: m_loci(numLoci), m_probMut(probMut), m_carrCap(carrCap)
{
	m_genotypes = 1u << m_loci;
	auto seed = static_cast<unsigned>(chrono::system_clock::now().time_since_epoch().count());
	m_generator = new default_random_engine(seed);
	m_realDist = new uniform_real_distribution<double>(0.0, 1.0);
	m_lociDist = new uniform_int_distribution<unsigned>(1, m_loci);
}

Simulation::~Simulation()
{
	delete m_generator;
	delete m_realDist;
	delete m_lociDist;
}


double Simulation::drand() // random double between 0.0 and 1.0
{
	return m_realDist->operator()(*m_generator);
}

int Simulation::lrand() // random locus
{
	return m_lociDist->operator()(*m_generator);
}

int Simulation::poisson(double mu) // sample from a poisson distribution
{
	if (mu <= 0.0)
		return 0;
	poisson_distribution<int> dist(mu);
	return dist(*m_generator);
}

vector<string> Simulation::getGenotypeStrings()
{
	vector<string> genotypes;
	for (int g = 0; g < m_genotypes; g++)
	{
		string s;
		for (int i = m_loci - 1; i >= 0; i--)
			(g & (1u << i)) ? s.push_back('1') : s.push_back('0');
		genotypes.push_back(s);
	}
	return genotypes;
}

void Simulation::m_growth(const vector<double> &r, vector<int64_t> &N)
{
	for (size_t i = 0; i < N.size(); i++)
	{
		double newSize = floor(N[i] * exp(r[i]));
		if (newSize > 0.0)
		{
			if (drand() < (N[i] * exp(r[i])) - newSize)
				newSize++;
			N[i] = static_cast<int64_t>(newSize);
		}
		else
			N[i] = 0;
	}
}

void Simulation::m_mutation(const vector<double> &r, vector<int64_t> &N)
{
	vector<int64_t> oldN(N);
	for (size_t i = 0; i < N.size(); i++)
	{
		double avgNumMutants = oldN[i] * m_probMut;
		int numMutants = poisson(avgNumMutants);
		for (int m = 0; m < numMutants; m++)
		{
			unsigned locus = lrand();
			size_t mutant = static_cast<unsigned>(i) ^ 1u << (locus - 1);
			N[mutant]++;
		}
		N[i] -= numMutants;
	}
}

void Simulation::m_death(const vector<double> &r, vector<int64_t> &N)
{
	int64_t sumN = 0;
	for (size_t i = 0; i < N.size(); i++)
		sumN += N[i];
	if (sumN == 0) return;
	for (size_t i = 0; i < N.size(); i++)
	{
		double freq = static_cast<double>(N[i]) / sumN;
		double newSize = floor(freq * m_carrCap);
		if (newSize > 0.0 && drand() < (freq * m_carrCap) - newSize)
			newSize++;
		N[i] = static_cast<int64_t>(newSize);
	}
}

void Simulation::m_stepForward(const vector<double> &r, vector<int64_t> &N)
{
	m_growth(r, N);
	m_mutation(r, N);
	m_death(r, N);
}

void Simulation::simpleSimulation(const vector<double> &r, vector<int64_t> &N, const int &t)
{
	trace.push_back(N);
	for (int i = 0; i < t; i++)
	{
		m_stepForward(r, N);
		trace.push_back(N);
	}
}

bool Simulation::saveTrace(const string &filename)
{
	ofstream fs(filename);
	if (fs.fail())
		return false;
	vector<string> genotypes = getGenotypeStrings();
	auto iter = genotypes.begin();
	for (; iter != genotypes.end() - 1; iter++)
		fs << *iter << ',';
	fs << *iter << '\n';
	for (auto timestep : trace)
	{
		auto iter = timestep.begin();
		for (; iter != timestep.end() - 1; iter++)
			fs << *iter << ',';
		fs << *iter << '\n';
	}
	fs.close();
	return true;
}

int main()
{
	vector<double> sampleLandscape = { 0.0, 0.011, 1.396, 0.82, 0.487, 0.963, 1.487, 1.471, 1.006, 1.234, 1.429, 1.0, 1.011, 1.077, 1.51, 1.393 };
	vector<int64_t> samplePopulation = { 1000000000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
	
	Simulation sim;

	vector<string> genotypes = sim.getGenotypeStrings();

	for (int i = 0; i < samplePopulation.size(); i++)
		cout << genotypes[i] << ' ' << samplePopulation[i] << endl;
	cout << "----------------" << endl;

	auto t1 = chrono::high_resolution_clock::now();
	sim.simpleSimulation(sampleLandscape, samplePopulation);
	auto t2 = chrono::high_resolution_clock::now();

	for (int i = 0; i < samplePopulation.size(); i++)
		cout << genotypes[i] << ' ' << samplePopulation[i] << endl;
	cout << "Simulation took " << chrono::duration_cast<chrono::milliseconds>(t2 - t1).count() << " ms" << endl;
	
	cout << sim.saveTrace("../data/cpp-out.csv") << endl;
	
	return 0;
}
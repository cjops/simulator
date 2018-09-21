#include "Simulator.h"

#include <cmath>
#include <fstream>

using namespace std;

Simulator::Simulator(int numLoci, int seedGen, double probMut, int64_t carrCap)
	: m_loci(numLoci)
	, m_genotypes(1u << m_loci)
	, m_probMut(probMut)
	, m_carrCap(carrCap)
	, m_lociDist(1, m_loci)
{
	if (!setPopulation(seedGen))
		setPopulation(DEFAULT_SEED_GEN);
}

vector<int64_t> Simulator::simpleSimulation(const vector<double> &landscape, int t)
{
	if (!setLandscape(landscape))
		return {};
	return simpleSimulation(t);
}

bool Simulator::setPopulation(const string &seedGen)
{
	return setPopulation(genotypeStringToInt(seedGen));
}

bool Simulator::setPopulation(int seedGen)
{
	if (seedGen < 0 || seedGen >= m_genotypes)
		return false;
	vector<int64_t> population(m_genotypes, 0);
	population[seedGen] = m_carrCap;
	m_population = population;
	return true;
}

bool Simulator::setPopulation(const vector<int64_t> &population)
{
	if (population.size() != m_genotypes)
		return false;
	for (auto i : population)
		if (i < 0)
			return false;
	m_population = population;
	return true;
}

bool Simulator::setLandscape(const vector<double> &landscape)
{
	if (landscape.size() < m_genotypes)
		return false;
	m_landscape = landscape;
	return true;
}

void Simulator::m_growth()
{
	auto &N = m_population; auto &r = m_landscape;
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

void Simulator::m_mutation()
{
	auto &N = m_population; auto &r = m_landscape;
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

void Simulator::m_death()
{
	auto &N = m_population; auto &r = m_landscape;
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

void Simulator::m_stepForward()
{
	m_growth();
	m_mutation();
	m_death();
}

vector<int64_t> Simulator::simpleSimulation(int t)
{
	if (m_landscape.empty())
		return {};
	m_trace.push_back(m_population);
	for (int i = 0; i < t; i++)
	{
		m_stepForward();
		m_trace.push_back(m_population);
	}
	return m_population;
}

bool Simulator::saveTrace(const string &filename)
{
	ofstream fs(filename);
	if (fs.fail())
		return false;
	vector<string> genotypes = getGenotypeStrings();
	auto iter = genotypes.begin();
	for (; iter != genotypes.end() - 1; iter++)
		fs << *iter << ',';
	fs << *iter << '\n';
	for (auto timestep : m_trace)
	{
		auto iter = timestep.begin();
		for (; iter != timestep.end() - 1; iter++)
			fs << *iter << ',';
		fs << *iter << '\n';
	}
	fs.close();
	return true;
}

vector<string> Simulator::getGenotypeStrings() const
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

int Simulator::genotypeStringToInt(const string &s)
{
	unsigned g = 0;
	unsigned sz = static_cast<unsigned>(s.size());
	for (unsigned i = 0; i < sz; i++)
		if (s[i] == '1')
			g |= 1u << (sz - i - 1);
	return g;
}

const vector<vector<int64_t>>& Simulator::getTrace() const
{
	return m_trace;
}

void Simulator::resetTrace()
{
	m_trace.clear();
}

int Simulator::getNumLoci() const
{
	return m_loci;
}

int Simulator::getNumGenotypes() const
{
	return m_genotypes;
}

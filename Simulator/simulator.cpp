#include "simulator.h"

#include <cmath>
#include <fstream>
#include <algorithm>
#include <iostream>

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
		if (newSize >= 0.0) // a negative growth rate might result in a negative
							// size which we must avoid
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
	auto &N = m_population;
	vector<int64_t> oldN(N);
	for (size_t i = 0; i < N.size(); i++)
	{
		double avgNumMutants = oldN[i] * m_probMut;
		int numMutants = poisson(avgNumMutants);
		for (int m = 0; m < numMutants; m++)
		{
			unsigned locus = lrand() - 1;
			size_t mutant = static_cast<unsigned>(i) ^ 1u << locus;
			if (m_firstAppearances[mutant] == -1)
				m_firstAppearances[mutant] = i;
			N[mutant]++;
		}
		N[i] -= numMutants;
	}
}

void Simulator::m_death()
{
	auto &N = m_population;
	int64_t sumN = 0;
	for (size_t i = 0; i < N.size(); i++)
		sumN += N[i];
	if (sumN == 0) return;
	for (size_t i = 0; i < N.size(); i++)
	{
		double freq = static_cast<double>(N[i]) / sumN;
		double newSize = floor(freq * m_carrCap);
		if (drand() < (freq * m_carrCap) - newSize)
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

void Simulator::m_setOptimalGenotype()
{
	m_optimalGenotype = max_element(m_landscape.begin(), m_landscape.end()) - m_landscape.begin();
}

void Simulator::m_resetCritTimes()
{
	fill(m_critTimes, m_critTimes+3, -1);
}

void Simulator::m_checkCritTimes()
{
	if (m_optimalGenotype == -1)
		return;
	int timestep = m_trace.size();
	if (m_critTimes[0] == -1 && m_population[m_optimalGenotype] > 0)
		m_critTimes[0] = timestep; // first appearance
	if (m_critTimes[1] == -1 && m_population[m_optimalGenotype] > 0.5 * m_carrCap)
		m_critTimes[1] = timestep; // dominance
	if (m_critTimes[2] == -1 && m_population[m_optimalGenotype] > 0.99 * m_carrCap)
		m_critTimes[2] = timestep; // fixation
}

void Simulator::m_resetFirstAppearances()
{
	vector<int> v(m_genotypes, -1);
	m_firstAppearances = v;
    for (int i = 0; i < m_population.size(); i++)
    {
        if (m_population[i] > 0)
            m_firstAppearances[i] = -2;
    }
}

vector<int64_t> Simulator::simpleSimulation(int t)
{
	if (m_landscape.empty())
		return {};
	m_setOptimalGenotype();
	m_resetCritTimes();
	m_resetFirstAppearances();
	m_setGreedyPath();
	m_trace.push_back(m_population);
	m_checkCritTimes();
	for (int i = 0; i < t; i++)
	{
		m_stepForward();
		m_trace.push_back(m_population);
		m_checkCritTimes();
	}
	m_setActualPath();
	return m_population;
}

vector<int64_t> Simulator::switchingSimulation(const switchingScheme &scheme, int opt)
{
	if (scheme.empty() || !setLandscape(scheme[0].second))
		return {};
	if (opt == -1)
		m_setOptimalGenotype();
	else
		m_optimalGenotype = opt;
	m_resetCritTimes();
	m_resetFirstAppearances();
	m_setGreedyPath();
	m_trace.push_back(m_population);
	m_checkCritTimes();
	for (auto const &iter : scheme)
	{
		setLandscape(iter.second);
		for (int i = 0; i < iter.first; i++)
		{
			m_stepForward();
			m_trace.push_back(m_population);
			m_checkCritTimes();
		}
	}
	m_setActualPath();
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

string Simulator::intToGenotypeString(int g) const
{
	string s;
	for (int i = m_loci - 1; i >= 0; i--)
		(g & (1u << i)) ? s.push_back('1') : s.push_back('0');
	return s;
}

vector<string> Simulator::getGenotypeStrings() const
{
	vector<string> genotypes;
	for (int g = 0; g < m_genotypes; g++)
		genotypes.push_back(intToGenotypeString(g));
	return genotypes;
}

int Simulator::genotypeStringToInt(const string &s) const
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

bool Simulator::setProbMut(double prob)
{
	if (prob > 1 || prob < 0)
		return false;
	m_probMut = prob;
	return true;
}

bool Simulator::setCarrCap(int64_t cap)
{
	if (cap < 0)
		return false;
	m_carrCap = cap;
	return true;
}

const int* Simulator::getCritTimes() const
{
	return m_critTimes;
}

void Simulator::m_setGreedyPath()
{
	if (m_optimalGenotype == -1)
		return;
	
	vector<int> path;
	double maxRate = 0.0;
	int g = 0;
	
	// start from genotype of nonzero abundance w/ highest growth rate
	// (most likely the seed)
	for (size_t i = 0; i < m_population.size(); i++)
	{
		if (m_population[i] > 0 && m_landscape[i] > maxRate)
		{
			maxRate = m_landscape[i];
			g = i;
		}
	}
	path.push_back(g);

	while (g != m_optimalGenotype)
	{
		if (path.size() > m_loci * 2)
			return; // avoid infinite loops
		maxRate = 0.0;
		for (int locus = 0; locus < m_loci; locus++)
		{
			int neighbor = static_cast<unsigned>(path.back()) ^ 1u << (locus);
			if (m_landscape[neighbor] > maxRate)
			{
				maxRate = m_landscape[neighbor];
				g = neighbor;
			}
		}
		if (m_landscape[g] < m_landscape[path.back()])
			break; // we have reached a local optimum
		path.push_back(g);
	}
	
	m_greedyPath = path;
}

void Simulator::m_setActualPath()
{
	vector<int> path;
	// start from most abundant genotype (most likely global optimum)
	int g = max_element(m_population.begin(), m_population.end()) - m_population.begin();
	path.push_back(g);
    g = m_firstAppearances[g];
	while (g != -1 && g != -2)
	{
		if (path.size() > m_loci * 2)
			return; // avoid infinite loops
		if (find(path.begin(), path.end(), g) != path.end())
			break; // we may not reach the seed
		path.push_back(g);
        g = m_firstAppearances[g];
	}
	reverse(path.begin(), path.end());
	m_path = path;
}

vector<int> Simulator::getActualPath() const
{
	return m_path;
}

vector<int> Simulator::getGreedyPath() const
{
	return m_greedyPath;
}

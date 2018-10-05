#ifndef SIMULATOR_INCLUDED
#define SIMULATOR_INCLUDED

#include <cstdint>
#include <string>
#include <vector>
#include <random>
#include <chrono>

#define RANDOM_SEED_VAL static_cast<unsigned>(std::chrono::system_clock::now().time_since_epoch().count())

const int		DEFAULT_NUM_LOCI	= 4;
const int		DEFAULT_SEED_GEN	= 0;
const double	DEFAULT_PROB_MUT	= 1.0e-8;
const int64_t	DEFAULT_CARR_CAP	= 1000000000;
const int		DEFAULT_TIMESTEPS	= 1200;

typedef std::vector<std::pair<int, std::vector<double>>> switchingScheme;

class Simulator
{
public:
	Simulator(
		int		numLoci = DEFAULT_NUM_LOCI,
		int		seedGen = DEFAULT_SEED_GEN,
		double	probMut = DEFAULT_PROB_MUT,
		int64_t carrCap = DEFAULT_CARR_CAP
	);
	std::vector<int64_t> simpleSimulation(
		const std::vector<double> &landscape,
		int t = DEFAULT_TIMESTEPS
	);
	std::vector<int64_t> simpleSimulation(int t = DEFAULT_TIMESTEPS);
	std::vector<int64_t> switchingSimulation(const switchingScheme &scheme, int opt=-1);
	bool saveTrace(const std::string &filename);
	double drand();
	int lrand();
	int poisson(double mu);
	std::string intToGenotypeString(int i) const;
	std::vector<std::string> getGenotypeStrings() const;
	int genotypeStringToInt(const std::string &s) const;
	bool setPopulation(const std::vector<int64_t> &population);
	bool setPopulation(int seedGen);
	bool setPopulation(const std::string &seedGen);
	bool setLandscape(const std::vector<double> &landscape);
	const std::vector<std::vector<int64_t>>& getTrace() const;
	void resetTrace();
	int getNumLoci() const;
	int getNumGenotypes() const;
	bool setProbMut(double prob);
	bool setCarrCap(int64_t cap);
	const int* getCritTimes() const;
	std::vector<int> getActualPath() const;
	std::vector<int> getGreedyPath() const;
private:
	const int m_loci;
	const int m_genotypes;
	double m_probMut;
	int64_t m_carrCap;
	std::vector<double> m_landscape;
	std::vector<int64_t> m_population;
	std::vector<std::vector<int64_t>> m_trace;
	int m_critTimes[3];
	void m_resetCritTimes();
	void m_checkCritTimes();
	int m_optimalGenotype = -1;
	void m_setOptimalGenotype();
	std::vector<int> m_path;
	std::vector<int> m_greedyPath;
	void m_setGreedyPath();
	std::vector<int> m_firstAppearances;
	void m_resetFirstAppearances();
	void m_setActualPath();

	std::default_random_engine m_generator{ RANDOM_SEED_VAL };
	std::uniform_real_distribution<double> m_realDist{ 0.0, 1.0 };
	std::uniform_int_distribution<unsigned> m_lociDist;

	void m_growth();
	void m_mutation();
	void m_death();
	void m_stepForward();
};

inline double Simulator::drand() // random double between 0.0 and 1.0
{
	return m_realDist(m_generator);
}

inline int Simulator::lrand() // random locus
{
	return m_lociDist(m_generator);
}

inline int Simulator::poisson(double mu) // sample from a Poisson distribution
{
	if (mu <= 0.0)
		return 0;
	std::poisson_distribution<int> dist(mu);
	return dist(m_generator);
}

#endif

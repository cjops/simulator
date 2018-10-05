#ifndef _SIMULATOR_NO_PYTHON
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#include <iostream>
#include <chrono>
#include <cmath>
#include "simulator.h"

using namespace std;

vector<double> sampleLandscape = { 0.0, 0.011, 1.396, 0.82, 0.487, 0.963, 1.487, 1.471, 1.006, 1.234, 1.429, 1.0, 1.011, 1.077, 1.51, 1.393 };
vector<int64_t> samplePopulation = { 1000000000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

static PyObject* getPyTrace(const Simulator &sim)
{
	PyObject* pDict = PyDict_New();
	const vector<vector<int64_t>> &trace = sim.getTrace();
	size_t numGen = sim.getNumGenotypes();
	size_t traceLength = trace.size();
	vector<string> genStrings = sim.getGenotypeStrings();
	for (int i = 0; i < numGen; i++)
	{
		int64_t* genArray = new int64_t[traceLength];
		for (int j = 0; j < traceLength; j++)
			genArray[j] = trace[j][i];
		npy_intp dims[] = { static_cast<npy_intp>(traceLength) };
		PyObject* pArray = PyArray_SimpleNewFromData(1, dims, NPY_INT64, genArray);
		PyArray_ENABLEFLAGS((PyArrayObject *)(pArray), NPY_ARRAY_OWNDATA);
		PyDict_SetItem(pDict, Py_BuildValue("s", genStrings[i].c_str()), pArray);
	}
	return pDict;
}

static PyObject *
simulator_run_simple(PyObject *self, PyObject *args)
{
	int timesteps;

	if (!PyArg_ParseTuple(args, "i", &timesteps))
		return NULL;
	Simulator sim;
	sim.simpleSimulation(sampleLandscape, timesteps);
	return getPyTrace(sim);
}

static PyObject *
simulator_simulate(PyObject *self, PyObject *args, PyObject *keywds)
{
	PyObject* landscapeObj;
    int timesteps = DEFAULT_TIMESTEPS;
    PyObject* populationObj = nullptr;
    PyObject* seedObj = nullptr;
    double	probMut = DEFAULT_PROB_MUT;
	int64_t carrCap = DEFAULT_CARR_CAP;
	PyObject* durationsObj = nullptr;
	int frequency = -1;

    static char* kwlist[] = {"landscapes", "timesteps", "starting_population", "seed",
                             "prob_mutation", "carrying_cap", "durations", "frequency", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|iOOdLOi", kwlist,
                                     &landscapeObj, &timesteps, &populationObj, &seedObj,
                                     &probMut, &carrCap, &durationsObj, &frequency))
        return NULL;
	
	// interpret landscape(s)
	
	Py_ssize_t numGenotypes;
	int numLoci;
	Py_ssize_t numLandscapes;
	bool switching;
	vector<double> landscape;
	switchingScheme scheme;
	
	if (!PyList_Check(landscapeObj) || PyList_Size(landscapeObj) == 0)
		return NULL;
	PyObject* firstItem = PyList_GetItem(landscapeObj, 0);
	if (PyList_Check(firstItem) && PyList_Size(landscapeObj) > 1)
	{ 	// multiple landscapes within a list
		switching = true;
		numGenotypes = PyList_Size(firstItem);
		numLoci = static_cast<int>(ceil(log2(numGenotypes)));
		numLandscapes = PyList_Size(landscapeObj);
		if (durationsObj && PyList_Size(durationsObj) == numLandscapes)
		{
			for (int i = 0; i < numLandscapes; i++)
			{
				PyObject* ls = PyList_GetItem(landscapeObj, i);
				cout << "landscape" << i << endl;
				int duration = PyLong_AsLong(PyList_GetItem(durationsObj, i));
				vector<double> landscape(numGenotypes);
				for (Py_ssize_t j = 0; j < numGenotypes; j++)
				{
					cout << PyFloat_AsDouble(PyList_GetItem(ls, j)) << " ";
					landscape[j] = PyFloat_AsDouble(PyList_GetItem(ls, j));
				}
				cout << endl;
				for (auto val : landscape)
					cout << val << " ";
				cout << endl;
				scheme.push_back({duration, {landscape}});
			}
		}
		else if (frequency != -1)
		{
			// to be implemented later
			return NULL;
		}
		else
			return NULL;
	}
	else
	{	// a single landscape
		switching = false;
		if (PyList_Check(firstItem)) // single landscape within a list
			landscapeObj = firstItem;
			
		numGenotypes = PyList_Size(landscapeObj);
		if (numGenotypes == 0)
			return NULL;
		numLoci = static_cast<int>(ceil(log2(numGenotypes)));
		landscape.resize(numGenotypes);
		for (Py_ssize_t i = 0; i < numGenotypes; i++)
			landscape[i] = PyFloat_AsDouble(PyList_GetItem(landscapeObj, i));
	}
		
	// initialize simulator object
	Simulator sim(numLoci);
	
	// interpret starting population
	
	if (populationObj)
	{
		vector<int64_t> population(numGenotypes);

		for (Py_ssize_t i = 0; i < numGenotypes; i++)
			population[i] = PyLong_AsLongLong(PyList_GetItem(populationObj, i));

		sim.setPopulation(population);
	}
	else if (seedObj)
	{
		if (PyLong_Check(seedObj))
		{
			int seed = PyLong_AsLong(seedObj);
			sim.setPopulation(seed);
		}
		else if (PyUnicode_Check(seedObj))
		{
			const char* seed = PyUnicode_AsUTF8(seedObj);
			sim.setPopulation(seed);
		}
	}
	
	// set constants
	
	sim.setProbMut(probMut);
	sim.setCarrCap(carrCap);
	
	// run simulation
	if (switching == true)
		sim.switchingSimulation(scheme);
	else
		sim.simpleSimulation(landscape, timesteps);
	
	// build results object
	
	const int* critTimes = sim.getCritTimes();
	vector<int> greedyPath = sim.getGreedyPath();
	vector<int> actualPath = sim.getActualPath();
		
	PyObject* results = PyDict_New();
	PyDict_SetItemString(results, "trace", getPyTrace(sim));
	PyDict_SetItemString(results, "T_1", Py_BuildValue("i", critTimes[0]));
	PyDict_SetItemString(results, "T_d", Py_BuildValue("i", critTimes[1]));
	PyDict_SetItemString(results, "T_f", Py_BuildValue("i", critTimes[2]));
	
	PyObject* greedy = PyList_New(greedyPath.size());
	for (size_t i = 0; i < greedyPath.size(); i++)
	{
		string item = sim.intToGenotypeString(greedyPath[i]);
		PyList_SetItem(greedy, i, Py_BuildValue("s", item.c_str()));
	}
	PyDict_SetItemString(results, "greedy_path", greedy);
	
	PyObject* actual = PyList_New(actualPath.size());
	for (size_t i = 0; i < actualPath.size(); i++)
	{
		string item = sim.intToGenotypeString(actualPath[i]);
		PyList_SetItem(actual, i, Py_BuildValue("s", item.c_str()));
	}
	PyDict_SetItemString(results, "actual_path", actual);
	
	return results;
}

static PyMethodDef SimulatorMethods[] = {
	{"run_simple",  simulator_run_simple, METH_VARARGS,
	 "Run a simple simulation."},
	 {"simulate", (PyCFunction)simulator_simulate, METH_VARARGS | METH_KEYWORDS,
     "Print a simulation with many options."},
	{NULL, NULL, 0, NULL}        /* Sentinel */
};

static struct PyModuleDef simulatormodule = {
	PyModuleDef_HEAD_INIT,
	"simulator",	/* name of module */
	NULL,			/* module documentation, may be NULL */
	-1,				/* size of per-interpreter state of the module,
					   or -1 if the module keeps state in global variables. */
	SimulatorMethods
};

PyMODINIT_FUNC
PyInit_simulator(void)
{
	import_array();
	return PyModule_Create(&simulatormodule);
}

int
main(int argc, char *argv[])
{
	wchar_t *program = Py_DecodeLocale(argv[0], NULL);
	if (program == NULL) {
		fprintf(stderr, "Fatal error: cannot decode argv[0]\n");
		exit(1);
	}

	/* Add a built-in module, before Py_Initialize */
	PyImport_AppendInittab("simulator", PyInit_simulator);

	/* Pass argv[0] to the Python interpreter */
	Py_SetProgramName(program);

	/* Initialize the Python interpreter.  Required. */
	Py_Initialize();

	/* Optionally import the module; alternatively,
	   import can be deferred until the embedded script
	   imports it. */
	PyImport_ImportModule("simulator");

		PyMem_RawFree(program);
	return 0;
}
#endif

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
simulator_simulate(PyObject *self, PyObject *args, PyObject *keywds, bool switching=false)
{
	PyObject* landscapeObj;
    int timesteps = DEFAULT_TIMESTEPS;
    PyObject* populationObj = nullptr;
    PyObject* seedObj = nullptr;
    double	probMut = DEFAULT_PROB_MUT;
	int64_t carrCap = DEFAULT_CARR_CAP;

    static char* kwlist[] = {"landscape", "timesteps", "starting_population", "seed",
                             "prob_mutation", "carrying_cap", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, keywds, "O|iOOdL", kwlist,
                                     &landscapeObj, &timesteps, &populationObj, &seedObj,
                                     &probMut, &carrCap))
        return NULL;
	
	Py_ssize_t numGenotypes;
	int numLoci;
	
	if (PyList_Check(landscapeObj))
	{
		numGenotypes = PyList_Size(landscapeObj);
		if (numGenotypes == 0)
			return NULL;
		numLoci = static_cast<int>(ceil(log2(numGenotypes)));
	}
	else
		return NULL;
	
	vector<double> landscape(numGenotypes);
	
	for (Py_ssize_t i = 0; i < numGenotypes; i++)
		landscape[i] = PyFloat_AsDouble(PyList_GetItem(landscapeObj, i));
	
	Simulator sim(numLoci);
	
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
	
	sim.setProbMut(probMut);
	sim.setCarrCap(carrCap);
	
	sim.simpleSimulation(landscape, timesteps);

	return getPyTrace(sim);
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

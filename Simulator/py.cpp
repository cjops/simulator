#ifndef _SIMULATOR_NO_PYTHON
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#include <iostream>
#include <chrono>
#include "Simulator.h"

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

static PyMethodDef SimulatorMethods[] = {
	{"run_simple",  simulator_run_simple, METH_VARARGS,
	 "Run a simple simulation."},
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

/*
 * Dummy C extension to force platform-specific wheel
 * This ensures the package is treated as platlib instead of purelib
 */

#include <Python.h>

static PyObject* dummy_function(PyObject* self, PyObject* args) {
    Py_RETURN_NONE;
}

static PyMethodDef dummy_methods[] = {
    {"dummy", dummy_function, METH_NOARGS, "Dummy function"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef dummy_module = {
    PyModuleDef_HEAD_INIT,
    "orca._dummy",
    "Dummy module to force platform wheel",
    -1,
    dummy_methods
};

PyMODINIT_FUNC PyInit__dummy(void) {
    return PyModule_Create(&dummy_module);
}

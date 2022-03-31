#include "Python.h"
#include <numpy/arrayobject.h>
#include "masHalos.h"

static PyObject* masHalos_masHalos(PyObject* self, PyObject* args){
	int ng, nH;
	double boxsize;
	/* numpy objects containing coordinates
	and density on the grid */
	PyObject *pos_obj, *d_obj, *Nx_obj;

    if (!PyArg_ParseTuple(args, "iidOOO", &ng, &nH, &boxsize, &pos_obj, &d_obj, &Nx_obj))
       	return NULL;

    /* Interpret the input objects as numpy arrays. */
    PyObject *pos_array = PyArray_FROM_OTF(pos_obj, NPY_FLOAT, NPY_IN_ARRAY);
    PyObject *d_array = PyArray_FROM_OTF(d_obj, NPY_FLOAT, NPY_IN_ARRAY);
    PyObject *Nx_array = PyArray_FROM_OTF(Nx_obj, NPY_FLOAT, NPY_IN_ARRAY);

    /* If that didn't work, throw an exception. */
    if ( pos_array == NULL || d_array == NULL || Nx_array == NULL ) {
        Py_XDECREF(pos_array);
        Py_XDECREF(d_array);
        Py_XDECREF(Nx_array);
        return NULL;
    }

    /* Get pointers to the data as C-types. */
    float *pos = (float*)PyArray_DATA(pos_array);
    float *density = (float*)PyArray_DATA(d_array);
    float *Nx = (float*)PyArray_DATA(Nx_array);

    masHalos( ng, nH, boxsize, pos, density, Nx);

    /* Clean up */
    Py_DECREF(pos_array);
    Py_DECREF(d_array);
    Py_DECREF(Nx_array);

    // I return something for index counter
	return Py_BuildValue("i", 1);
}

static PyMethodDef myMethods[] = {
    {"masHalos",	masHalos_masHalos,	METH_VARARGS, "calculate mas of halo field with cic"},
    {NULL, NULL, 0, NULL}  // always
};

#if PY_MAJOR_VERSION < 3
/* --- Python 2 --- */
PyMODINIT_FUNC init_masHalos(void){
    (void) Py_InitModule("_masHalos", myMethods);

    /* Load `numpy` functionality. */
    import_array();
}

#else
/* --- Python 3 --- */
static struct PyModuleDef mas_module = {
  PyModuleDef_HEAD_INIT,
  "mas", 	/* m_name */
  NULL,         /* m_doc */
  -1,           /* m_size */
  myMethods,    /* m_methods */
};

PyMODINIT_FUNC
PyInit__masHalos(void)
{

  PyObject *m = NULL;

  /* Load `numpy` functionality. */
  import_array();

  /* Module initialization  */
  m = PyModule_Create(&mas_module);
  if (m == NULL) goto bad;

  return m;

 bad:
  return NULL;
}

#endif

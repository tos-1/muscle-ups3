#include "Python.h"
#include <numpy/arrayobject.h>
#include "halomodel.h"

// bind function to Python
static PyObject* halomodel_halomodel(PyObject* self, PyObject* args)
{
    int ng;
    double boxsize, redshift, Delta_v, mcrit, rho;
    // increase reference count
    PyObject *sift_obj, *psi3_obj, *cc_obj, *pos_obj, *path_obj, *hmf_obj;

    if ( !PyArg_ParseTuple(args, "idddddOOOOOO", &ng, &boxsize, &redshift, &Delta_v, &mcrit, &rho, 
			&sift_obj, &psi3_obj, &cc_obj, &pos_obj, &path_obj, &hmf_obj) )
      return NULL;

    /* Interpret the input objects as numpy arrays. */
    PyObject *sift_array = PyArray_FROM_OTF( sift_obj, NPY_INT32 , NPY_IN_ARRAY);
    PyObject *cc_array = PyArray_FROM_OTF( cc_obj, NPY_FLOAT32 , NPY_IN_ARRAY);
    PyObject *psi3_array = PyArray_FROM_OTF( psi3_obj, NPY_INT32 , NPY_IN_ARRAY);
    PyObject *pos_array = PyArray_FROM_OTF( pos_obj, NPY_FLOAT64 , NPY_IN_ARRAY);
    PyObject *path_string = PyObject_Str(path_obj);
    PyObject *hmf_string = PyObject_Str(hmf_obj);

    /* If that didn't work, throw an exception, then decrement the reference count.
    Py_XDECREF accounts for the possibility of NULL objects */
    if ( sift_array == NULL || cc_array == NULL || psi3_array == NULL || pos_array == NULL || path_string == NULL || hmf_string == NULL ){
      Py_XDECREF( sift_array );
      Py_XDECREF( cc_array );
      Py_XDECREF( psi3_array );
      Py_XDECREF( pos_array );
      Py_XDECREF( path_string );
      Py_XDECREF( hmf_string );
      return NULL;
    }

    /* Get pointers to the data as C-types. */
    int *sift = (int*)PyArray_DATA(sift_array);
    int *psi3 = (int*)PyArray_DATA(psi3_array);
    float *cc = (float*)PyArray_DATA(cc_array);
    double *pos = (double*)PyArray_DATA(pos_array);

    #if PY_MAJOR_VERSION >= 3
    const char *path = PyUnicode_AsUTF8(path_string);
    const char *hmf  = PyUnicode_AsUTF8( hmf_string);
    #else
    const char *path = PyString_AsString(path_string);
    const char *hmf  = PyString_AsString( hmf_string);
    #endif

    int nhalos = halomodel( ng, boxsize, redshift, Delta_v, mcrit, rho, sift, psi3, cc, pos, path, hmf );

    /* Reference count clean up */
    Py_DECREF( sift_array );
    Py_DECREF( psi3_array );
    Py_DECREF( cc_array );
    Py_DECREF( pos_array );
    Py_DECREF( path_string );
    Py_DECREF( hmf_string );

    PyObject *ret = Py_BuildValue("i", nhalos);
    return ret;
}

static PyMethodDef myMethods[] = {
    {"halomodel", halomodel_halomodel, METH_VARARGS, "build halo catalogue and implement halo model"},
    {NULL, NULL, 0, NULL}  // always
};


#if PY_MAJOR_VERSION < 3
/* --- Python 2 --- */
PyMODINIT_FUNC init_halomodel(void)
{
    (void) Py_InitModule("_halomodel", myMethods);

    /* Load `numpy` functionality. */
    import_array();
}

#else
/* --- Python 3 --- */
static struct PyModuleDef hm_module = {
  PyModuleDef_HEAD_INIT,
  "halomodel", 	/* m_name */
  NULL,         /* m_doc */
  -1,           /* m_size */
  myMethods,    /* m_methods */
};

PyMODINIT_FUNC
PyInit__halomodel(void)
{

  PyObject *m = NULL;

  /* Load `numpy` functionality. */
  import_array();

  /* Module initialization  */
  m = PyModule_Create(&hm_module);
  if (m == NULL) goto bad;

  return m;

 bad:
  return NULL;
}

#endif

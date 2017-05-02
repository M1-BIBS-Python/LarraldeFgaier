Backend libraries
=================

Fully supported
^^^^^^^^^^^^^^^

*Every component implements this backend*

NumPy
-----

`numpy` provides a powerful yet flexible implementation of n-dimension arrays,
and is used as a base by most of the other libraries. It is required by
`dockerasmus.pdb`, since `Protein` instances have properties and methods
that directly return `numpy.ndarray` objects.

NumPy can only perform computations on the CPU.


Theano
------

`theano` is a library which can defined, evaluate, compile and optimise
n-dimension array computations. It works seamlessly with NumPy array, but
is a lot more efficient than interpreted computations while adding only
a light overhead. Theano is very useful to compute the same mathematical
expression a lot of times, which is what happend with a scoring function.

Theano can perform computations on the GPU using a configuration file.
Read the `official documentation <http://deeplearning.net/software/theano/tutorial/using_gpu.html>`_
to find out how.


Partially supported
^^^^^^^^^^^^^^^^^^^

*Only some components implement these backend*

MXNet
-----

`mxnet` is a deep-learning library developed that supports both symbolic and
imperative programming in order to perform efficient computations. It is
maintained by the developers previously behind `minerva`, `cxxnet` and
`purine2`.


Tensorflow
----------




Considered
^^^^^^^^^^

*Could possibly be implemented as backends*

* PyTorch
* PyCUDA

.. hint::

  Know a library not listed here that works well for your other projects ?
  Submit a feature request to the `issue tracker <https://gitlab.com/althonos/dockerasmus/issues>`_
  or even better, fork the project, try implementing a backend, add a test,
  and submit a `merge request <https://gitlab.com/althonos/dockerasmus/merge_requests>`_ !

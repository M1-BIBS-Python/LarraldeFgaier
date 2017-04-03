Dockerasmus
=======================================
*docking itself through university*

|docs| |build| |coverage|

.. |docs| image:: http://readthedocs.org/projects/dockerasmus/badge/?version=latest
   :target: http://dockerasmus.readthedocs.io/en/latest/?badge=latest

.. |build| image:: https://gitlab.com/althonos/dockerasmus/badges/master/build.svg
   :target: https://gitlab.com/althonos/dockerasmus/pipelines?scope=branches

.. |coverage| image:: https://codecov.io/gl/althonos/dockerasmus/branch/master/graph/badge.svg?token=eNxJwF5lhn
   :target: https://codecov.io/gl/althonos/dockerasmus




Dockerasmus is a Python version-agnostic module that was created to
quickly solve docking problems, as part of a Python assignment from
the M1 BIBS of the Université Paris-Saclay.

Dockerasmus provides a generic implementation of a scoring function,
which can be used with several *components* to compute the score of
a docking conformation of two proteins. It is backend agnostic, and
any component can be rewritten with another library.










.. toctree::
   :maxdepth: 2

   api/index


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

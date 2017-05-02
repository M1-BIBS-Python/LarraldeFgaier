dockerasmus
=======================================
*docking itself through university*

|docs| |gl| |build| |coverage| |license| |grade|

.. |docs| image:: http://readthedocs.org/projects/dockerasmus/badge/?version=latest
   :target: http://dockerasmus.readthedocs.io/en/latest/?badge=latest

.. |build| image:: https://gitlab.com/althonos/dockerasmus/badges/master/build.svg
   :target: https://gitlab.com/althonos/dockerasmus/pipelines?scope=branches

.. |coverage| image:: https://img.shields.io/codecov/c/token/eNxJwF5lhn/gitlab/althonos/dockerasmus.svg
   :target: https://codecov.io/gl/althonos/dockerasmus

.. |gl| image:: https://img.shields.io/badge/repository-gitlab-orange.svg
   :target: https://gitlab.com/althonos/dockerasmus

.. |license| image:: https://img.shields.io/badge/license-GPLv3-blue.svg
   :target: https://choosealicense.com/licenses/gpl-3.0/

.. |grade| image:: https://api.codacy.com/project/badge/Grade/87e03271b04e4154a3b513bddb0d08bc
   :target: https://www.codacy.com/app/althonos/dockerasmus



Introduction
------------

``dockerasmus`` is a version-agnostic Python module that was created
to quickly solve docking problems, as part of a Python assignment from
the M1 BIBS of the Université Paris-Saclay.

``dockerasmus`` provides a generic implementation of a scoring function,
which can be used with several *components* to compute the score of
a docking conformation involving two proteins. It is backend agnostic, and
every scoring component can be rewritten with any library supporting
``numpy`` arrays.



Example
-------

Use ``dockerasmus`` to compute the score of the barnase-barstar complex
using the scoring function defined by `Cornell et al
<http://dx.doi.org/10.1021/ja00124a002>`_:

.. code:: python

    from dockerasmus.pdb import Protein
    from dockerasmus.score import ScoringFunction, components

    # Import the pdb files (supports gzipped files or plain .pdb)
    barnase = Protein.from_pdb_file("tests/data/barnase.native.pdb.gz")
    barstar = Protein.from_pdb_file("tests/data/barstar.native.pdb.gz")

    # Create a scoring function with two components
    scoring_function = ScoringFunction(components.LennardJones,
                                       components.Coulomb)

    # Call the scoring function on the barnase (receptor)
    # and the barstar (ligand)
    scoring_function(barnase, barstar)  # -84.94...


API
---

``dockerasmus`` provides several submodules:

* a parser & object model for the Protein Data Bank (``dockerasmus.pdb``)
* a scoring library (``dockerasmus.score``)
* a soft 3D engine for spatial transformations (``dockerasmus.spatial``)

See the `API reference <http://dockerasmus.readthedocs.io/en/latest/api/>`_
from the online documentation to get more details.


License
-------

``dockerasmus`` is fully open-source and is released under the GPLv3.

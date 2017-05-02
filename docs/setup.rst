Setup
=====

With PyPI
---------

``dockerasmus`` is distributed on PyPI, which means it is possible to
install through `pip <https://pip.pypa.io/>`_, the PyPA recommended tool
for managing Python packages.

To install the whole module, simply run one of the following commands:

.. code-block:: console

    # pip install dockerasmus         # install on the whole system
    $ pip install --user dockerasmus  # install for current user only


This will install the module with its requirements (`numpy` and `six`).

To use additional requirements (that make computations within
``dockerasmus``) more efficient, install `theano` and `scipy` with
either one of the following commands:

.. code-block:: console

    # pip install theano scipy
    $ pip install --user theano scipy



With GitLab
-----------

.. warning:

   The releases distributed on PyPI are stable snapshots of the project,
   following semantic versioning. The same cannot be said of the git
   repository, which hosts the development version of `dockerasmus`.
   Use with caution as they may be more unstable.

The easiest is to clone the repository, and then install it with PyPI
(this is preferable to installing the module with `setuptools` as the
`setuptools` install uses the deprecated *egg* format.)

.. code-block:: console

    $ git clone https://gitlab.com/althonos/dockerasmus
    # pip install .         # install on the whole system
    $ pip install --user .  # install for current user only

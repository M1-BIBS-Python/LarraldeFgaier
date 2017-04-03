# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import collections
import six
import copy
import gzip
import functools
import string
import numpy

from ..utils import iterators
from .chain import Chain
from .residual import Residual
from .atom import Atom


class Protein(collections.OrderedDict):
    __slots__ = ("id", "name", "_atom_charges")

    _CMAP_MODES = {
        'mass_center': lambda r1,r2: r1.distance_to(r2.mass_center),
        'nearest': lambda r1, r2: min(a1.distance_to(a2.pos)
                for a1 in r1.itervalues() for a2 in r2.itervalues()),
        'farthest': lambda r1, r2: max(a1.distance_to(a2.pos)
                for a1 in r1.itervalues() for a2 in r2.itervalues()),
    }

    @staticmethod
    def _parse_pdb_atom_line(line):
        """Return a raw `dict` with atom properties from a pdb atom line.

        Returns:
            dict: a dictionary which keys are: ``serial``, ``name``,
                ``chainID``, ``altLoc``, ``resName``, ``resSeq``,
                ``iCode``, ``x``, ``y`` and ``z``.
        """
        schema = {'serial': (6, 11), 'name': (12, 16), 'altLoc': (16, 17),
                  'resName': (17, 20), 'chainID': (21, 22), 'resSeq': (22, 26),
                  'iCode': (26, 27), 'x': (30, 38), 'y': (38, 46), 'z': (46, 54)}
        # Decode a binary string to a unicode/str object
        decode = lambda s: s.decode('utf-8')
        # callback to be called after the value field  is isolated from the line,
        # either to transtype or to decode a binary string
        callbacks = {'serial': int, 'name': decode, 'altLoc': decode,
                    'resName': decode, 'chainID': decode, 'resSeq': int,
                    'iCode': decode, 'x': float, 'y': float, 'z': float}
        return {key: callbacks.get(key)(line[i:j].strip())
                for key,(i,j) in schema.items()}

    @classmethod
    def from_pdb(cls, handle):
        """Create a new Protein object from a PDB file handle.

        Arguments:
            handle (file handle): a file-like object opened in
                binary read mode (must be line-by-line iterable).
        """

        protein = cls()
        for line in handle:
            if line.startswith(b"ATOM  "):

                atom = cls._parse_pdb_atom_line(line)

                if not atom['chainID'] in protein:
                    protein[atom['chainID']] = Chain(atom['chainID'])

                if not atom['resSeq'] in protein[atom['chainID']]:
                    protein[atom['chainID']][atom['resSeq']] = Residual(atom['resSeq'], atom['resName'])

                protein[atom['chainID']][atom['resSeq']][atom['name']] = Atom(
                    atom['x'], atom['y'], atom['z'], atom['serial'], atom['name'],
                    protein[atom['chainID']][atom['resSeq']],
                )
        return protein

    @classmethod
    def from_pdb_file(cls, path):
        """Create a new Protein object from a PDB file.

        Arguments:
            path (`str`): the path to a PDB protein file (supports gzipped
                or plain text PDB files).
        """
        if path.endswith('.gz'):
            open_function = gzip.open
        else:
            open_function = functools.partial(open, mode='rb')
        with open_function(path) as pdb_file:
            return cls.from_pdb(pdb_file)

    def __init__(self, id=None, name=None, chains=None):
        """Create a new Protein object.

        Arguments:
            id (`int`): the id of the protein.
            name (`str`): the name of the protein.
            chains (`dict` of `Chain`): a dictionary of the chains
                of the proteins referenced by their ``id``.
        """
        super(Protein, self).__init__(chains or {})
        self.id = id
        self.name = name

        # Memoize matrices and vectors
        self._atom_charges = None
        self._atom_epsilon = None
        self._atom_positions = None
        self._atom_radius = None

    def __contains__(self, item):
        if isinstance(item, Chain):
            return super(Protein, self).__contains__(item)
        elif isinstance(item, (Residual, Atom)):
            return any(item in chain for chain in self)
        elif isinstance(item, six.text_type):
            return any(item == chain.id for chain in self.itervalues())
        else:
            raise TypeError(
                "'in <Protein>' requires Chain, Residual, Atom or unicode"
                " as left operand, not {}".format(type(item).__name__)
            )

    def __getitem__(self, item):
        """Overloaded __getitem__ allowing slicing and individual atom access.

        Example:
            >>> complex = Protein.from_pdb_file("tests/data/1brs.pdb.gz")
            >>> barstar = complex[u'D':]
            >>> sorted(barstar.keys())
            [u'D', u'E', u'F']
            >>> complex[1]
            Atom 1(16.783, 48.812, 26.447)
            >>> barstar[1]
            Traceback (most recent call last):
               ...
            KeyError: u'Could not find Atom with id: 1'
        """
        if isinstance(item, slice):
            stop = item.stop or iterators.nth(iterators.wordrange(max(self.keys())), 1)
            start = item.start or min(self.keys())
            return Protein(chains={
                k:super(Protein, self).__getitem__(k)
                    for k in iterators.wordrange(start, stop)
            })
        elif isinstance(item, int):
            atom = next((atom for atom in self.iteratoms() if atom.id==item), None)
            if atom is None:
                raise KeyError("Could not find Atom with id: {}".format(item))
            return atom
        else:
            return super(Protein, self).__getitem__(item)

    @property
    def mass(self):
        """The mass of the protein.

        Warning:
            Computed as the sum of the masses of the residuals
            of the chain (it does not take the masses of the atoms
            in the peptidic bound into account).
        """
        return sum(chain.mass for chain in self.itervalues())

    @property
    def mass_center(self):
        r"""The position of mass center of the protein.

        .. math::

           mc &= \sum_{i}{\frac{w_i}{W}
             \begin{pmatrix} x_i \\ y_i \\ z_i \end{pmatrix}
           }

        where :math:`i` is the index of each atom, :math:`x_i`
        (resp. :math:`y_i`) (resp. :math:`z_i`) if the abcissa
        (resp. ordinate) (resp. height) of the atom :math:`i` in the
        worldspace, and :math:`W = \sum_i{w_i}` the approximated mass
        of the whole protein.

        Warning:
            Uses `Protein.mass`, so only the atoms on the residuals
            of each aminoacid are used for the computation.
        """
        mass = self.mass
        return sum(
            (atom.mass/mass)*atom.pos for chain in self.itervalues()
            for res in chain.itervalues() for atom in res.itervalues()
        )

    @property
    def radius(self):
        """The radius of the sphere the protein would fit in.

        Equals to the norm of the position of the atom of the protein farthest
        from its mass center.
        """
        origin = self.mass_center
        return max(
            atom.distance_to(origin)
                for chain in self.itervalues()
                    for residual in chain.itervalues()
                        for atom in residual.itervalues()
        )

    def atom_charges(self):
        """The vector of the charge of each atom of the protein.
        """
        if self._atom_charges is None:
            self._atom_charges = numpy.array([
                a.charge for a in self.iteratoms()
            ])
        return self._atom_charges

    def atom_epsilon(self):
        """The vector of the epsilon value of each atom of the protein.
        """
        if self._atom_epsilon is None:
            self._atom_epsilon = numpy.array([
                a.epsilon for a in self.iteratoms()
            ])
        return self._atom_epsilon

    def atom_positions(self):
        """The matrix of the positions of each atom of the protein.
        """
        if self._atom_positions is None:
            self._atom_positions = numpy.array([
                a.pos for a in self.iteratoms()
            ])
        return self._atom_positions

    def atom_radius(self):
        """The vector of the Van der Waals radius of each atom of the protein.
        """
        if self._atom_radius is None:
            self._atom_radius = numpy.array([
                a.radius for a in self.iteratoms()
            ])
        return self._atom_radius

    def contact_map(self, other, mode='nearest'):
        """Return a 2D contact map between residuals of ``self`` and ``other``.

        Arguments:
            other (`Protein`): the other protein with which to create
                a contact map (chains/residuals/atoms must have the same
                names in both proteins)

        Keyword Arguments:
            mode (`str`): how to compute the contact map. Available modes are:
              ``'nearest'`` (the distance between the two closest atoms of
              the two residuals), ``'farthest'`` (the distance between the two
              farthest atoms of the two residuals) or ``'mass_center'``
              (the distance between the mass center of the two residuals).
        """
        if not mode in self._CMAP_MODES:
            raise ValueError("Unknown mode: '{}'".format(mode))

        dim_x = max(r.id for chain in self.itervalues() for r in chain.itervalues())
        dim_y = max(r.id for chain in other.itervalues() for r in chain.itervalues())

        cmap = numpy.zeros((dim_x+1, dim_y+1))

        for c in self.itervalues():
            for r in c.itervalues():
                for other_c in other.itervalues():
                    for other_r in other_c.itervalues():
                        cmap[r.id, other_r.id] = self._CMAP_MODES[mode](r, other_r)

        return cmap

    def copy(self):
        """Return a deep copy of ``self``.
        """
        return Protein(self.id, self.name, {
            chain.id: Chain(chain.id, chain.name, {
                residual.id: copy.deepcopy(residual)
                    for residual in chain.itervalues()
            })
                for chain in self.itervalues()
        })

    def iteratoms(self):
        """Yield every atom in ``self``.

        Yields:
            `dockerasmus.pdb.Atom`: every atom of the protein,
            ordered by the id of their chain and the id
            of their residual.
        """
        for chain in self.itervalues():
            for residual in chain.itervalues():
                for atom in residual.itervalues():
                    yield atom

    def nearest_atom(self, pos):
        """Return the atom nearest to the position ``pos``.
        """
        return min(
            (atom for chain in self.values() for res in chain.values() for atom in res.values()),
            key = lambda a: a.distance_to(pos)
        )

    if six.PY3:
        def itervalues(self):
            return six.itervalues(self)

        def iteritems(self):
            return six.iteritems(self)

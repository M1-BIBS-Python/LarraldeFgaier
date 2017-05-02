# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals
from __future__ import division

import collections
import six
import copy
import gzip
import functools
import numpy

from ..utils import iterators
from .chain import Chain
from .residue import Residue
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
            `dict`: a dictionary which keys are: ``serial``, ``name``,
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

                if atom['chainID'] not in protein:
                    protein[atom['chainID']] = Chain(atom['chainID'])

                if atom['resSeq'] not in protein[atom['chainID']]:
                    protein[atom['chainID']][atom['resSeq']] = Residue(atom['resSeq'], atom['resName'])

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
        self._atom_pwd = None
        self._atom_positions = None
        self._atom_radius = None

    def __add__(self, other):
        """Return a new Protein complexed with ``other``.

        Example:
            >>> sorted(barnase.keys())
            [u'B']
            >>> complex = barnase + barstar
            >>> sorted(complex.keys())
            [u'B', u'D']

        Raises:
            ValueError: when there's an ID collision with
                the chain or the protein
            TypeError: when ``other`` is neither a `Chain`
                nor a `Protein`.
        """
        return self.copy().__iadd__(other)

    def __iadd__(self, other):
        """Complex the Protein with ``other``.

        Example:
            >>> prot = Protein()
            >>> prot += barnase
            >>> sorted(prot.keys())
            [u'B']
            >>> prot += barstar
            >>> sorted(prot.keys())
            [u'B', u'D']

        Raises:
            ValueError: when there's an ID collision with
                the chain or the protein
            TypeError: when ``other`` is neither a `Chain`
                nor a `Protein`.
        """
        if isinstance(other, Chain):
            if other.id in self:
                raise ValueError(
                    'Protein already contains a chain'
                    'with id {} !'.format(other.id)
                )
            self[other.id] = other
        elif isinstance(other, Protein):
            if not set(self.keys()).isdisjoint(other.keys()):
                common_keys = set(self.keys()).intersection(other.keys())
                raise ValueError(
                    'Protein already contains a chain with id'
                    ' {} !'.format(', '.join(common_keys))
                )
            self.update(other)
        else:
            raise TypeError(
                "unsupported operand type(s) for +: 'Protein'"
                " and '{}'".format(type(other).__name__)
            )
        return self


    def __contains__(self, item):
        if isinstance(item, six.text_type):
            return super(Protein, self).__contains__(item)
        elif isinstance(item, (Residue, Atom)):
            return any(item in chain for chain in self.itervalues())
        elif isinstance(item, Chain):
            return super(Protein, self).__contains__(item.id)
        else:
            raise TypeError(
                "'in <Protein>' requires Chain, Residue, Atom or unicode"
                " as left operand, not {}".format(type(item).__name__)
            )

    def __getitem__(self, item):
        """Overloaded __getitem__ allowing slicing and individual atom access.

        Example:
            >>> complex = Protein.from_pdb_file("tests/data/1brs.pdb.gz")
            >>> barstar = complex[u'D':]
            >>> sorted(barstar.keys())
            [u'D', u'E', u'F']
        """
        if isinstance(item, slice):
            stop = item.stop or iterators.nth(iterators.wordrange(max(self.keys())), 1)
            start = item.start or min(self.keys())
            return Protein(chains=collections.OrderedDict([
                (k, super(Protein, self).__getitem__(k))
                    for k in iterators.wordrange(start, stop)
            ]))
        else:
            return super(Protein, self).__getitem__(item)

    @property
    def mass(self):
        """The mass of the protein.

        Warning:
            Computed as the sum of the masses of the residues
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
            Uses `Protein.mass`, so only the atoms on the residues
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
                    for residue in chain.itervalues()
                        for atom in residue.itervalues()
        )

    def atom_charges(self):
        """The vector of the charge of each atom of the protein.
        """
        if self._atom_charges is None:
            self._atom_charges = numpy.array([
                a.charge for a in self.iteratoms()
            ])
        return self._atom_charges

    def atom_pwd(self):
        """The vector of the potential well depth of each atom of the protein.
        """
        if self._atom_pwd is None:
            self._atom_pwd = numpy.array([
                a.pwd for a in self.iteratoms()
            ])
        return self._atom_pwd

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
        """Return a 2D contact map between residues of ``self`` and ``other``.

        Arguments:
            other (`Protein`): the other protein with which to create
                a contact map (chains/residues/atoms must have the same
                names in both proteins)

        Keyword Arguments:
            mode (`str`): how to compute the contact map. Available modes are:
              ``'nearest'`` (the distance between the two closest atoms of
              the two residues), ``'farthest'`` (the distance between the two
              farthest atoms of the two residues) or ``'mass_center'``
              (the distance between the mass center of the two residues).
        """
        if mode not in self._CMAP_MODES:
            raise ValueError("Unknown mode: '{}'".format(mode))
        if not isinstance(other, Protein):
            raise TypeError("other must be a Protein,"
                            " not {}".format(type(other).__name__))

        dim_x = max(r.id for chain in self.itervalues() for r in chain.itervalues())
        dim_y = max(r.id for chain in other.itervalues() for r in chain.itervalues())

        cmap = numpy.zeros((dim_x+1, dim_y+1))

        for c in self.itervalues():
            for r in c.itervalues():
                for other_c in other.itervalues():
                    for other_r in other_c.itervalues():
                        cmap[r.id, other_r.id] = self._CMAP_MODES[mode](r, other_r)

        return cmap

    def atom(self, atom_id):
        """Get atom of ``self`` with id ``atom_id``.

        Raises:
            KeyError: when no Atom has the given id.
        """
        try:
            atom_id = int(atom_id)
        except ValueError:
            six.raise_from(
                TypeError("Invalid type: {}".format(type(atom_id).__name__)),
                None
            )

        atom_id = int(atom_id)
        atom = next((atom for atom in self.iteratoms() if atom.id==atom_id), None)
        if atom is None:
            raise KeyError("Could not find Atom with id: {}".format(atom_id))
        return atom

    def residue(self, res_id):
        """Get residue of ``self`` with id ``res_id``.

        Raises:
            KeyError: when no Residue has the given id.
        """
        try:
            res_id = int(res_id)
        except ValueError:
            six.raise_from(
                TypeError("Invalid type: {}".format(type(res_id).__name__)),
                None
            )

        res = next((res for chain in self.itervalues()
                    for res in chain.itervalues() if res.id==res_id), None,
        )
        if res is None:
            raise KeyError("Could not find Residue with id: {}".format(res_id))
        return res

    def copy(self):
        """Return a deep copy of ``self``.
        """
        return Protein(self.id, self.name, collections.OrderedDict([
            (chain.id, Chain(chain.id, chain.name, collections.OrderedDict([
                (residue.id, copy.deepcopy(residue))
                    for residue in chain.itervalues()
            ])))
                for chain in self.itervalues()
        ]))

    def iteratoms(self):
        """Yield every atom in ``self``.

        Yields:
            `Atom`: every atom of the protein, ordered by
            the id of their chain and the id of their residue.
        """
        for chain in self.itervalues():
            for residue in chain.itervalues():
                for atom in sorted(residue.itervalues(), key=lambda a: a.id):
                    yield atom

    def interface(self, other, distance=4.5):
        """Return couples of residues of ``self`` interfacing with ``other``.

        Yields:
            tuple - a (`Residue`, `Residue`) couple where the first element
                is a residue of `self` and the second a residue of `other`
                only if two or more of their respective aminoacids are closer
                than ``distance``.
        """
        if not isinstance(other, Protein):
            raise TypeError("Invalid type: {}".format(type(other).__name__))

        # Reuse the distance matrix computation of dockerasmus.score !
        from ..score.requirements import distance as dist
        mx_distance = dist(self, other)

        # Atoms of self and other in the same order as
        # self.atom_positions() and other.atom_position
        atoms_self = list(self.iteratoms())
        atoms_other = list(other.iteratoms())

        # Start the iterator
        done = set()
        for i in range(mx_distance.shape[0]):
            for j in range(mx_distance.shape[1]):
                if mx_distance[i,j] < distance:
                    res_self, res_other = atoms_self[i].residue, atoms_other[j].residue
                    if (res_self, res_other) not in done:
                        done.add((res_self, res_other))
                        yield res_self, res_other


    def nearest_atom(self, pos):
        """Return the atom nearest to the position ``pos``.
        """
        return min(
            (atom for chain in self.values() for res in chain.values() for atom in res.values()),
            key = lambda a: a.distance_to(pos)
        )

    def rmsd(self, other):
        """
        """
        #TODO: matricial computation
        rmsd, length = 0, 0
        if isinstance(other, Protein):
            # Compute pairwise RMSD
            positions_other = other.atom_positions()
            positions_self = self.atom_positions()
            length = len(positions_self)
            if length != len(positions_other):
                raise ValueError(
                "'other' does not have the same number of atoms !"
                )
            rmsd = numpy.sum((positions_self - positions_other)**2)
        elif isinstance(other, (list, numpy.ndarray)):
            # Compute reference-wise RMSD
            if len(other) != 3:
                raise ValueError(
                    "Invalid reference position dimension: {}".format(len(other))
                )
            vec_distance_squared = (self.atom_positions() - other)**2
            rmsd = numpy.sum(vec_distance_squared)
            length = len(vec_distance_squared)
        else:
            raise TypeError("other must be Protein, list or numpy.ndarray,"
                            " not {}".format(type(other).__name__))
        return (rmsd/length)**.5

    if six.PY3:
        def itervalues(self):
            return six.itervalues(self)

        def iteritems(self):
            return six.iteritems(self)

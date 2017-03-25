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

from .utils import parse_pdb_atom_line
from .chain import Chain
from .residual import Residual
from .atom import Atom


class Protein(collections.OrderedDict):

    _CMAP_MODES = {
        'mass_center': lambda r1,r2: r1.distance_to(r2.mass_center),
        'nearest': lambda r1, r2: min(a1.distance_to(a2.pos)
                for a1 in r1.itervalues() for a2 in r2.itervalues()),
        'farthest': lambda r1, r2: max(a1.distance_to(a2.pos)
                for a1 in r1.itervalues() for a2 in r2.itervalues()),
    }

    @classmethod
    def from_pdb(cls, handle):
        """Create a new Protein object from a PDB file handle.

        Arguments:
            handle (file handle): a file-like object opened in
                binary read mode.
        """

        protein = cls()
        for line in handle:
            if line.startswith(b"ATOM  "):

                atom = parse_pdb_atom_line(line)

                if not atom['chainID'] in protein:
                    protein[atom['chainID']] = Chain(atom['chainID'])

                if not atom['resSeq'] in protein[atom['chainID']]:
                    protein[atom['chainID']][atom['resSeq']] = Residual(atom['resSeq'], atom['resName'])

                protein[atom['chainID']][atom['resSeq']][atom['name']] = Atom(
                    atom['x'], atom['y'], atom['z'], atom['serial'], atom['name'], atom['resName'][-3:],
                )
        return protein

    @classmethod
    def from_pdb_file(cls, path):
        """Create a new Protein object from a PDB file.

        Arguments:
            path (str): the path to a PDB protein file (supports gzipped
                PDB files).
        """
        if path.endswith('.gz'):
            open_function = gzip.open
        else:
            open_function = functools.partial(open, mode='rb')
        with open_function(path) as pdb_file:
            return cls.from_pdb(pdb_file)

    def __init__(self, prot_id=None, prot_name=None, chains=None):
        super(Protein, self).__init__(chains or {})
        self.id = prot_id
        self.name = prot_name

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
        """Overloaded dict.__getitem__ allowing chain slicing

        Exemple:
            >>> complex = Protein.from_pdb_file("tests/data/1brs.pdb.gz")
            >>> barstar = complex['D':]
            >>> list(barstar.keys())
            ['D', 'E', 'F']
        """
        if isinstance(item, slice):
            stop = item.stop or utils.nth(utils.infinitewords(max(self.keys())), 1)
            start = item.start or min(self.keys())
            return Protein(chains={
                k:super(Protein, self).__getitem__(k)
                    for k in utils.infinitewords(start, stop)
            })
        else:
            return super(Protein, self).__getitem__(item)

    if six.PY3:
        def itervalues(self):
            return six.itervalues(self)

        def iteritems(self):
            return six.iteritems(self)

    def iteratoms(self):
        for chain in self.itervalues():
            for residual in chain.itervalues():
                for atom in residual.itervalues():
                    yield atom

    def copy(self):
        """Returns a deep copy of self
        """
        return Protein(self.id, self.name, {
            chain.id: Chain(chain.id, chain.name, {
                residual.id: copy.deepcopy(residual)
                    for residual in chain.itervalues()
            })
                for chain in self.itervalues()
        })

    @property
    def mass(self):
        return sum(chain.mass for chain in self.itervalues())

    @property
    def mass_center(self):
        mass = self.mass
        return sum(
            (atom.mass/mass)*atom.pos for chain in self.itervalues()
            for res in chain.itervalues() for atom in res.itervalues()
        )

    def contact_map(self, other, mode='nearest'):
        """Return a 2D contact map between residuals of `self` and `other`

        Arguments:
            other (Protein): the other protein with which to create
                a contact map (chains/residuals/atoms must have the same
                names in both proteins)

        Keyword Arguments:
            mode (str): how to compute the contact map. Available methods are:
                * 'nearest': the distance between the two closest atoms of
                    the two residuals. **[Default]**
                * 'farthest': the distance between the two farthest atoms
                    of the two residuals.
                * 'mass_center': the distance between the mass center of the
                    two residuals.
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

    @property
    def radius(self):
        """The radius of the sphere the protein would fit in

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

    def nearest_atom(self, pos):
        """Returns the atom nearest to the position `pos`
        """
        return min(
            (atom for chain in self.values() for res in chain.values() for atom in res.values()),
            key = lambda a: a.distance_to(pos)
        )

# coding: utf-8
from __future__ import absolute_import
from __future__ import unicode_literals

import collections
import six


def parse_pdb_atom_line(line):

    schema = {
        #'record name': (0, 6),
        'serial': (6, 11),
        'name': (12, 16),
        'altLoc': (16, 17),
        'resName': (17, 20),
        'chainID': (21, 22),
        'resSeq': (22, 26),
        'iCode': (26, 27),
        'x': (30, 38),
        'y': (38, 46),
        'z': (46, 54),
    }

    # Decode a binary string to a unicode/str object
    decode = lambda s: s.decode('utf-8')

    # callback to be called after the value field
    # is isolated from the line, either to transtype
    # or to decode a binary string
    callbacks = {
        'serial': int,
        'name': decode,
        'altLoc': decode,
        'resName': decode,
        'chainID': decode,
        'resSeq': int,
        'iCode': decode,
        'x': float,
        'y': float,
        'z': float,
    }

    return {
        key: callbacks.get(key)(line[i:j].strip())
            for key,(i,j) in schema.items()
    }



def parse(fh):

    # Use a defaultdict to avoid checking if the chain
    # dictionnary exists
    protein = collections.defaultdict(dict)

    for line in fh:
        if line.startswith("ATOM  "):

            atom = parse_pdb_atom_line(line)

            # Create a new residual dictionary if needed
            if not atom['resSeq'] in protein[atom['chainID']]:
                protein[atom['chainID']][atom['resSeq']] = {
                    'resname': atom['resName'], 'atomlist': [],
                }
                # Add the residual to the chain's reslist
                protein[atom['chainID']]['reslist'].append(atom['resSeq'])

            # Create the atom dictionnary with positions and id
            protein[atom['chainID']][atom['resSeq']][atom['name']] = {
                'x': atom['x'], 'y': atom['y'], 'z': atom['z'], 'id': atom['serial']
            }

            # Add the atom to the residual's atomlist
            protein[atom['chainID']][atom['resSeq']]['atomlist'].append(atom['name'])


    # Checks the chains
    protein['chains'] = list(protein)

    # Return a non-default dictionary
    return dict(protein)

#!/usr/bin/env python
# coding: utf-8
"""
Usage:
    score-conformations.py -i IN -r REF [-o OUT] [-q]
    score-conformations.py (-h | --help)

Positional Arguments:
    -i IN, --input IN           The path to the directory
                                containing ligand PDB files.
    -r REF, --reference REF     The path to the PDB file of
                                the receptor.
    -o OUT, --output OUT        The directory in which to
                                create the subfolders in which
                                to export the results.
                                [default: .]

Optional Arguments:
    -h, --help                  Print this message.
    -q, --quiet                 Do not display any output.
"""
# stdlib imports
import sys
import os
import glob
import operator
import shutil

# update sys.path to make dockerasmus importable
# locally although it is in the parent directory
SCRIPTDIR = os.path.dirname(os.path.abspath(__file__))
MAINDIR = os.path.dirname(SCRIPTDIR)
sys.path.insert(0, MAINDIR)

# Try importing non standard dependencies
try:
    import docopt
except ImportError:
    sys.exit('Could not import docopt - is it installed ?')
try:
    import progressbar
except ImportError:
    sys.exit('Could not import progressbar2 - is it installed ?')


# Import dockerasmus
from dockerasmus.pdb import Protein
from dockerasmus.score import ScoringFunction, components


if __name__ == "__main__":

    # Parse the CLI arguments
    args = docopt.docopt(__doc__)
    ref = Protein.from_pdb_file(args["--reference"])
    indir = args['--input']
    outdir = args['--output']

    # Create the output folder if needed
    for subdir in ('scoring_Cornell', 'scoring_maison'):
        if not os.path.isdir(os.path.join(outdir, subdir)):
            os.mkdir(os.path.join(outdir, subdir))

    # Create the scoring function
    sf_cornell = ScoringFunction(components.LennardJones, components.Coulomb)
    sf_custom = ScoringFunction(components.Fabiola, components.ScreenedCoulomb,
                                weights=[4, 1])

    # Compute the score for every protein in the input directory
    scores_cornell = {}
    scores_custom = {}

    # Get a list of files to score
    files = glob.glob(os.path.join(indir, '*.pdb'))

    # Create a progressbar if not quiet
    if not args["--quiet"]:
        pb = progressbar.ProgressBar()
        files = pb(files)

    # Score each file
    for filename in files:
        prot = Protein.from_pdb_file(filename)
        scores_cornell[os.path.basename(filename)] = sf_cornell(ref, prot)
        scores_custom[os.path.basename(filename)] = sf_custom(ref, prot)

    ###########
    # CORNELL #
    ###########

    # Write the ordered results in a .tsv file
    with open(os.path.join(outdir, 'scoring_Cornell', 'scores.tsv'), 'w') as f:
        f.write("file\tscore\n")
        for filename, score in sorted(scores_cornell.items(), key=operator.itemgetter(1)):
            f.write("{}\t{}\n".format(filename, score))

    # Copy the best result (lower is better !)
    shutil.copy(
        os.path.join(indir, min(scores_cornell.keys(), key=scores_cornell.get)),
        os.path.join(outdir, 'scoring_Cornell', "complexe_predit_score1.pdb")
    )

    ##########
    # CUSTOM #
    ##########

    # Write the ordered results in a .tsv file
    with open(os.path.join(outdir, 'scoring_maison', 'scores.tsv'), 'w') as f:
        f.write("file\tscore\n")
        for filename, score in sorted(scores_custom.items(), key=operator.itemgetter(1)):
            f.write("{}\t{}\n".format(filename, score))

    # Copy the best result (lower is better !)
    shutil.copy(
        os.path.join(indir, min(scores_custom, key=scores_custom.get)),
        os.path.join(outdir, 'scoring_maison', "complexe_predit_score2.pdb")
    )

#!/usr/bin/env python
# released under the GNU General Public License version 3.0 (GPLv3)

import os
import setuptools


def format_for_setup(requirement_file):
    """Build a list of requirements out of requirements.txt files.
    """
    requirements = []
    with open(requirement_file) as rq:
        for line in rq:
            line = line.strip()
            if line.startswith('-r'):
                other_requirement_file = line.split(' ', 1)[-1]
                requirements.extend(format_for_setup(os.path.join(
                    os.path.dirname(requirement_file),
                    other_requirement_file
                )))
            elif line:
                requirements.append(line)
    return requirements


## SETUPTOOLS VERSION
setuptools.setup(
    install_requires= format_for_setup('requirements.txt'),
    tests_require = format_for_setup(os.path.join('tests', 'requirements.txt')),
    extras_require = {
        'docs': format_for_setup(os.path.join('docs', 'requirements.txt')),
    },
    test_suite="tests",
)

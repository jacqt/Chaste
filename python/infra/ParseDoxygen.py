#!/usr/bin/env python


"""Copyright (c) 2005-2015, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""

"""
Parse the output of running Doxygen to find files that aren't properly documented.

The script takes arguments:
 <log_file>       Doxygen's normal output log.
 <error_log_file> Doxygen's error log.
 <output_dir>     The directory in which to generate summary files and
                  an index page.
"""

import os
import re
import sys

def munge_name(file_name, status, output_dir):
    """Return the file name to use for processed results.

    file_name is the name as it appears in the Doxygen logs.  This will be
    an absolute path.  We strip off the path to the Chaste root (i.e. make
    the name relative) and convert os.path.sep to '-'.  This assumes we are
    run from the Chaste root directory.

    The status code is then added, and our output_dir added on the front.
    """
    munged_file_name = file_name.replace(os.path.sep, '-')
    root_dir = os.getcwd().replace(os.path.sep, '-') + '-'
    if munged_file_name.startswith(root_dir):
        munged_file_name = munged_file_name[len(root_dir):]
    return os.path.join(output_dir, munged_file_name + '.' + status + '.0')

def parse_doxygen(log_file_name, error_log_file_name, output_dir):
    """Python function for script

    Parse the output of running Doxygen to find files that aren't properly documented.
    """

    # Remove any old output files/test results from output_dir
    for filename in os.listdir(output_dir):
        os.remove(os.path.join(output_dir, filename))

    # Get a list of source files
    source_files = set()
    for line in open(log_file_name, 'U'):
        if line.startswith('Parsing file '):
            source_files.add(line[13:-4]) # NB: line ends with a newline

    # Now work out which have problems
    error_re = re.compile(r'^(.+):(\d+):')
    file_name = None
    problem_files = {}
    counts = {}
    for line in open(error_log_file_name):
        m = error_re.match(line)
        if m:
            if m.group(1)[0] != '<':
                file_name = m.group(1)
            if not file_name in problem_files:
                problem_files[file_name] = []
                counts[file_name] = 1
            else:
                counts[file_name] += 1
        if file_name:
            problem_files[file_name].append(line)

    # Write output files for those with errors
    for problem_file, lines in problem_files.iteritems():
        content = ''.join(lines)
        if '<autogenerated>' in content:
            # Ignore autogenerated source files.  Note that not even an 'OK'
            # file will be created for them.
            continue
        status = '%d_%d' % (counts[problem_file], counts[problem_file])
        output_file = open(munge_name(problem_file, status, output_dir), 'w')
        output_file.write(content)
        output_file.close()

    # Now write 'OK' files for those without any errors
    ok_files = source_files - set(problem_files.iterkeys())
    for file_name in ok_files:
        output_file = open(munge_name(file_name, 'OK', output_dir), 'w')
        output_file.close()

    # And generate a summary page
    os.system('python python/DisplayTests.py '+output_dir+' DoxygenCoverage')

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print "Syntax error."
        print "Usage:",sys.argv[0],"<log file> <error log file> <test output dir>"
        sys.exit(1)

    parse_doxygen(sys.argv[1], sys.argv[2], sys.argv[3])

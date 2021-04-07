"""==========
This script extracts a single log-scaled trajectory
from a number of shorter log-scaled trajectories "fort.1xx"
and removes duplicated configurations from the lin-scaled trajectory.
"""
# **********************************************************************
# *** LAMELLAR DI-BLOCK COPOLYMER + NANODIMERS *************************
# *** Python script - data curation/minimalisation *********************
# *** Author: J. Javier Burgos Mármol -- 2021  *************************
# **********************************************************************
#BSD 3-Clause License
#
#Copyright (c) 2021, J. Javier Burgos-Mármol (Github username: jjavier-bm)
#All rights reserved.
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are met:
#
#1. Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
#2. Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
#3. Neither the name of the copyright holder nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
#CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
#OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
#OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

import argparse
import os
import math

def create_argument_parser():
    """Create a parser for the command line arguments used in crops-renumber"""

    parser = argparse.ArgumentParser(prog="LamAnalysis - uploadable trajectories",
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description="Script to extract single trajectories from high-stats files")
    parser.add_argument("lintraj_path",nargs=1, metavar="Sequence_filepath",
                        help="Input linear trajectory filepath.")
    parser.add_argument("Nparticles",nargs=1, type=int, metavar="N_particles",
                        help="Input number of nanodimers.")

    return parser

def main():
    Nconfigsperfile = 3
    Ntimesteps = 100000000
    period = 500000
    Nch = 12696
    Lch = 30
    fortfiles = [101,102,103,104,105,106,107,108,109,110,
                 120,130,140,150,160,170,180,190,200]

    parser = create_argument_parser()
    args = parser.parse_args()
    Npart = args.Nparticles[0]

    inlinpath = os.path.abspath(args.lintraj_path[0])
    if not os.path.isfile(inlinpath):
        raise OSError('input file not found')

    fileparts = []
    while True:
        if len(fileparts) == 0:
            parts = os.path.splitext(inlinpath)
            outlinpath = parts[0] + '.cleanversion' + parts[1] 
        else:
            parts = os.path.splitext(parts[0])
        if parts[1]=='.lin':
            outlogpath = parts[0] + '.log'
            for part in reversed(fileparts):
                outlogpath += part
            break
        else:
            fileparts.append(parts[1])
            if len(fileparts)==1:
                fileparts.append('.singlestat')
    # Linearly-(time-)spaced trajectory file. Remove duplicated configurations.
    with open(inlinpath,'r') as fin:
        with open(outlinpath,'w') as fout:
            l = 1
            c = 0
            while True:
                line=fin.readline().rstrip()
                if not line:
                    break
                if c > 0 and c % Nconfigsperfile == 0:
                    pass
                else:
                    fout.write(line+'\n')
                if l == 9+Nch*Lch+Npart*2:
                    l = 1
                    c += 1
                else:
                    l += 1

    # Logarithmically-(time-)spaced trajectory files. Single stat trajectory file.
    Nconfigsperfile = 9*(int(math.log10(Ntimesteps))-2)+1
    with open(outlogpath,'w') as fout:
        for n in fortfiles:
            inpath = os.path.join(os.path.dirname(inlinpath),'fort.'+str(n))
            with open(inpath,'r') as fin:
                l = 1
                c = 0 if n == 101 else 1
                while True:
                    line=fin.readline().rstrip()
                    if not line:
                        break
                    if n == 101 or c == Nconfigsperfile:
                        fout.write(line+'\n')
                    if l == 9+Nch*Lch+Npart*2:
                        l = 1
                        c += 1
                    else:
                        l += 1

    return

if __name__ == "__main__":
    import sys
    import traceback

    try:
        main()
        sys.exit(0)
    except Exception as e:
        if not isinstance(e, SystemExit):
            msg = "".join(traceback.format_exception(*sys.exc_info()))
            raise RuntimeError(msg)
        sys.exit(1)

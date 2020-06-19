from sys import argv as sys_argv
from argparse import ArgumentParser as argparse_ArgumentParser
from os import path as os_path
from os import makedirs as os_makedirs

# from sys import path as sys_path
# sys_path.insert(0, '/home/rpCache')

from Selenzy import analyse, preLoad

def add_arguments(parser):
    parser.add_argument('rxn',
                        help='Input reaction [default = rxn file]')
    parser.add_argument('-tar', type=float, default=20,
                        help='Number of targets to display in results [default = 20]')
    parser.add_argument('-d', type=float, default=0,
                        help='Use similiarity values for preferred reaction direction only [default=0 (OFF)]')
    parser.add_argument('datadir',
                        help='specify data directory for required databases files, please end with slash')
    parser.add_argument('outdir',
                        help='specify output directory for all output files, including final CSV file')
    parser.add_argument('-outfile',
                        help='specify non-default name for CSV file output')
    parser.add_argument('-NoMSA', action='store_true',
                        help='Do not compute MSA/conservation scores')
    parser.add_argument('-smarts', action='store_true',
                        help='Input is a reaction SMARTS string')
    parser.add_argument('-smartsfile', action='store_true',
                        help='Input is a reaction SMARTS file')
    parser.add_argument('-host', type=str, default='83333',
                        help='Host organism taxon id [default: E. coli]')
    # parser.add_argument('-sm', type=str, default='file',
    #                     help='Store mode for data [default: file]')
    return parser


def build_parser():
    parser = argparse_ArgumentParser(description='SeqFind script for Selenzy')
    parser = add_arguments(parser)
    return parser


def entrypoint(args=sys_argv[1:]):
    parser = build_parser()

    params = parser.parse_args(args)

    newpath = os_path.join(params.outdir)
    if not os_path.exists(newpath):
        os_makedirs(newpath)

    if params.smarts:
        rxnInput = ['-smarts', params.rxn]
    elif params.smartsfile:
        rxnInput = ['-smartsfile', params.rxn]
    else:
        rxnInput = ['-rxn', params.rxn]


    pc = preLoad(params.datadir)
    analyse(rxnInput,
            params.tar, params.datadir, params.outdir+'/', params.outfile,
            params.d, params.host,
            NoMSA=params.NoMSA,
            pc=pc)


##
#
#
if __name__ == "__main__":

    entrypoint(sys_argv[1:])

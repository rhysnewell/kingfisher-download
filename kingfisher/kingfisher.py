#!/usr/bin/env python3

__author__ = "Ben Woodcroft"
__copyright__ = "Copyright 2021"
__credits__ = ["Ben Woodcroft"]
__license__ = "GPL3+"
__maintainer__ = "Ben Woodcroft"
__email__ = "b.woodcroft near qut.edu.au"
__status__ = "Development"

import argparse
import os
import sys
import logging

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path
import kingfisher
from kingfisher.ena import DEFAULT_LINUX_ASPERA_SSH_KEY_LOCATION, DEFAULT_OSX_ASPERA_SSH_KEY_LOCATION


class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

class MyParser(argparse.ArgumentParser):
    FULL_HELP_FLAG = 'full-help'

    def format_help(self):
        prog = self.prog
        detail_str = 'See %s --%s for further options and further detail.\n' % (
            prog,
            self.FULL_HELP_FLAG)
        splits = prog.split()
        if len(splits) == 2:
            subcommand = splits[1]
        else:
            subcommand = None
        if subcommand == 'pipe':
            pass
            #TODO: Add short help with colour, see singlem for examples.

        else:
            return argparse.ArgumentParser.format_help(self)


def add_extraction_args(parser):
    parser.add_argument(
        '--output-format-possibilities', '--output_format_possibilities',
        nargs='+',
        help='Allowable output formats. If more than one is specified, downloaded data will processed as little as possible.',
        choices=['sra', 'fastq', 'fastq.gz','fasta','fasta.gz'],
        default=kingfisher.DEFAULT_OUTPUT_FORMAT_POSSIBILITIES)
    parser.add_argument(
        '--force',
        action='store_true',
        help='Re-download / extract files even if they already exist [default: Do not].')
    parser.add_argument(
        '--unsorted',
        action='store_true',
        help='Currently requires "--output-format-possibilities fasta" and "--stdout" and download from NCBI rather than ENA. \
            Output the sequences in a single file, where the reads are in arbitrary order. \
            Even pairs of reads may be split up, but it is possible to tell which pair \
            is which, and which is a forward and which is a reverse read from the name \
            [default: Do not].')
    parser.add_argument(
        '--stdout',
        action='store_true',
        help='Output sequences to STDOUT. Currently requires --unsorted [default: Do not].')
    return parser

def main():
    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('--debug', help='output debug information', action="store_true")
    parent_parser.add_argument('--version', help='output version information and quit',  action='version', version=kingfisher.__version__)
    parent_parser.add_argument('--quiet', help='only output errors', action="store_true")

    parser = argparse.ArgumentParser(parents=[parent_parser])
    subparsers = parser.add_subparsers(title="Sub-commands", dest='subparser_name')#, parser_class=MyParser)
    subparser_name_to_parser = {}

    def new_subparser(subparsers, parser_name, parser_description):
        subpar = subparsers.add_parser(parser_name,
                                       description=parser_description,
                                       help=parser_description,
                                       epilog=__author__,
                                       parents=[parent_parser])
        # subpar.add_argument(
        #     '--%s' % MyParser.FULL_HELP_FLAG,
        #     '--%s' % MyParser.FULL_HELP_FLAG.replace('-','_'),
        #     action='store_true', default=False, help='display all help options')
        subparser_name_to_parser[parser_name] = subpar
        return subpar

    get_description = 'Download and extract public sequence files from SRA or ENA'
    get_parser = new_subparser(subparsers, 'get', get_description)

    get_parser.add_argument(
        '--run-identifier','--run_identifier','-r',
        help='Run number to download/extract e.g. ERR1739691',
        required=True)
    get_parser.add_argument(
        '-m','--download_methods', '--download-methods',
        nargs='+',
        help='How to download .sra file. If multiple are specified, each is tried in turn until one works.',
        choices=['aws-http', 'prefetch', 'aws-cp', 'gcp-cp', 'ena-ascp','ena-ftp'], required=True)

    get_parser = add_extraction_args(get_parser)

    get_parser.add_argument(
        '--gcp-project', '--gcp_project',
        help='Downloading from Google Cloud buckets require a Google project to charge '
        '(they are requester-pays) e.g. \'my-project\'. This can alternately be set '
        'beforehand using \'gcloud config set project PROJECT_ID\' '
        '[default: value of `gcloud config get-value project` command]')
    get_parser.add_argument(
        '--gcp-user-key-file', '--gcp_user_key_file',
        help='Downloading from Google Cloud buckets requires a Google user to be setup. \
            Use this option to specify a JSON-formatted service account key, as per \
            https://cloud.google.com/iam/docs/creating-managing-service-account-keys \
            [default: not used]')
    get_parser.add_argument(
        '--aws-user-key-id', '--aws_user_key_id',
        help='Downloading from AWS requester pays buckets requires a key ID and secret key \
            [default: not used]')
    get_parser.add_argument(
        '--aws-user-key-secret', '--aws_user_key_secret',
        help='Downloading from AWS requester pays buckets requires a key ID and secret key \
            [default: not used]')
    get_parser.add_argument(
        '--allow-paid', '--allow_paid',
        help='Allow downloading from retriever-pays s3 and GCP buckets',
        action='store_true')
    get_parser.add_argument(
        '--allow-paid-from-gcp', '--allow_paid_from_gcp',
        help='Allow downloading from retriever-pays GCP buckets',
        action='store_true')
    get_parser.add_argument(
        '--allow-paid-from-aws', '--allow_paid_from_aws',
        help='Allow downloading from retriever-pays AWS buckets',
        action='store_true')
    get_parser.add_argument(
        '--ascp-ssh-key', '--ascp_ssh_key',
        help='\'linux\' or \'osx\' for default paths used in each OS ({} and {} respectively), \
            otherwise a path to the openssh key to used for aspera (i.e. the \
            -i flag of ascp) [default: \'{}\']'.format(
                DEFAULT_LINUX_ASPERA_SSH_KEY_LOCATION,
                DEFAULT_OSX_ASPERA_SSH_KEY_LOCATION,
                kingfisher.DEFAULT_ASPERA_SSH_KEY),
        default=kingfisher.DEFAULT_ASPERA_SSH_KEY)
    get_parser.add_argument(
        '--ascp-args', '--ascp_args',
        help='extra arguments to pass to ascp e.g. \'-k 2\' to resume with a \
        sparse file checksum [default: \'\']',
        default='')

    extract_description = 'extract .sra format files'
    extract_parser = new_subparser(subparsers, 'extract', extract_description)

    extract_parser.add_argument(
        '--sra',
        help='Extract this SRA file [required]',
        required=True)

    extract_parser = add_extraction_args(extract_parser)


    if (len(sys.argv) == 1 or sys.argv[1] == '-h' or sys.argv[1] == '--help'):
        print('')
        print('                ...::: Kingfisher v' + kingfisher.__version__ + ' :::...''')
        print('')
        print('   get      -> %s' % get_description)
        print('   extract  -> %s' % extract_description)


        print('\n  Use kingfisher <command> -h for command-specific help.\n'\
            '  Some commands also have an extended --full-help flag.\n')
        sys.exit(0)
    else:
        # Determine whether --full_help was specified before argument parsing.
        if ('--%s' % MyParser.FULL_HELP_FLAG in sys.argv or \
            '--%s' % MyParser.FULL_HELP_FLAG.replace('-','_') in sys.argv) and \
            sys.argv[1] in subparser_name_to_parser:
            subcommand = sys.argv[1]
            subparser = subparser_name_to_parser[subcommand]
            print(argparse.ArgumentParser.format_help(subparser))
            sys.exit(0)
        else:
            args = parser.parse_args()

    if args.debug:
        loglevel= logging.DEBUG
    elif args.quiet:
        loglevel= logging.ERROR
    else:
        loglevel= logging.INFO
    logging.basicConfig(
        level=loglevel, format='%(asctime)s %(levelname)s: %(message)s',
        datefmt='%m/%d/%Y %I:%M:%S %p')

    if args.subparser_name == 'get':
        kingfisher.download_and_extract(
            run_identifier = args.run_identifier,
            download_methods = args.download_methods,
            output_format_possibilities = args.output_format_possibilities,
            force = args.force,
            unsorted = args.unsorted,
            stdout = args.stdout,
            gcp_project = args.gcp_project,
            gcp_user_key_file = args.gcp_user_key_file,
            aws_user_key_id = args.aws_user_key_id,
            aws_user_key_secret = args.aws_user_key_secret,
            allow_paid = args.allow_paid,
            allow_paid_from_gcp = args.allow_paid_from_gcp,
            allow_paid_from_aws = args.allow_paid_from_aws,
            ascp_ssh_key = args.ascp_ssh_key,
            ascp_args = args.ascp_args
        )

    elif args.subparser_name == 'extract':
        output_files = kingfisher.extract(
            sra_file = args.sra,
            output_format_possibilities = args.output_format_possibilities,
            force = args.force,
            unsorted = args.unsorted,
            stdout = args.stdout,
        )
        logging.info("Output files: {}".format(', '.join(output_files)))
    else:
        raise Exception("Programming error")

    logging.info("Kingfisher done.")

if __name__ == '__main__':
    main()

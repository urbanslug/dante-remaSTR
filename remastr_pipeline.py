import argparse
import glob
import os
import subprocess


def load_arguments() -> argparse.Namespace:
    """
    Loads and parses arguments.
    :return: args - parsed arguments
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description='Python script to run the whole Dante REMASTR pipeline for several files')

    # add arguments
    parser.add_argument('bam_files', help='Glob path to bam file(s) that needs to be processed')
    parser.add_argument('nomenclature_file', help='Path to nomenclature file')
    parser.add_argument('output_dir', help='Path to directory where the output will be stored. Default=<current_dir>/results', nargs='?',
                        default=None)
    parser.add_argument('--dante-repo', help='Path to Dante repository. Default=use this location', default=None)
    parser.add_argument('--reference-fasta', help='Path to reference fasta. Default=/data/genome/human/grch38_decoy/grch38_decoy.fa',
                        default='/data/genome/human/grch38_decoy/grch38_decoy.fa')
    parser.add_argument('--param-file', type=str, help='Parameter file for inference of alleles. Default=MiSeq parameters from repo',
                        default=None)
    parser.add_argument('--skip-subpath', action='store_true', help='Do not append sub-path to output as <nomenclature_name>/<sample_name>')
    parser.add_argument('--dante-args', '-a', help='Send these additional arguments to the dante_remastr call. Default=None', default='')
    parser.add_argument('--use-old-defaults', '-d', action='store_true', help='Use old defaults for dante_remastr call')
    parser.add_argument('--visualize', '-v', action='store_true', help='Print all the outputs. Default is to print only the result table to stdout.')
    parser.add_argument('--skip-call', '-s', action='store_true', help='Skip the calling of scripts, just do a dry run.')
    parser.add_argument('--ignore-errors', '-e', action='store_true', help='Ignore errors and continue, otherwise end after first non 0 return code.')

    args = parser.parse_args()

    if args.dante_repo is None:
        args.dante_repo = os.path.dirname(os.path.realpath(__file__))
    if args.output_dir is None:
        args.output_dir = f'{os.getcwd()}/results'
    if args.param_file is None:
        args.param_file = f'{args.dante_repo}/src/inference/training_input/miseq_params.params'

    return args


def run_command(command: str, exit_on_error: bool, skip_call: bool):
    """
    Prints and runs the command
    :param command: str - command to run
    :param exit_on_error: bool - whether to exit on error
    :param skip_call: bool - dry run? do not call the command
    """
    print(command)
    if not skip_call:
        return_code = subprocess.call(command, shell=True)
        if exit_on_error and return_code != 0:
            exit(return_code)


if __name__ == '__main__':
    # load arguments
    args = load_arguments()

    # create the output dir if needed
    os.makedirs(args.output_dir, exist_ok=True)

    # select bam files
    bam_files = glob.glob(args.bam_files)

    print(f'Working on {len(bam_files)} samples...')
    print('conda activate remastr\n')

    nomenclature_name = os.path.basename(args.nomenclature_file).split('.')[0]

    # go through the bam files
    for bam_file in bam_files:
        # build remastr call
        sample_name = os.path.basename(bam_file).split('.')[0]
        output_dir = args.output_dir if args.skip_subpath else f'{args.output_dir}/{nomenclature_name}/{sample_name}'
        remastr_call = (f'{args.dante_repo}/remastr/target/release/remastr -f {args.reference_fasta} -m {args.nomenclature_file} '
                        f'-b {bam_file} -o {output_dir}/{sample_name}.tsv')

        # build dante_remastr call (report)
        additional_params = '--processes 24 --min-flank-len 5 --min-rep-len 6 --min-rep-cnt 2' if args.use_old_defaults else ''
        additional_params += args.dante_args + (' -v' if args.visualize else '')
        dante_call = (f'python {args.dante_repo}/dante_remastr.py --param-file {args.param_file} {additional_params} -o {output_dir} < '
                      f'{output_dir}/{sample_name}.tsv > {output_dir}/{sample_name}_res.tsv')

        # print and call commands
        if not args.skip_call:
            os.makedirs(output_dir, exist_ok=True)
        run_command(remastr_call, not args.ignore_errors, args.skip_call)
        run_command(dante_call, not args.ignore_errors, args.skip_call)

    # gather results
    output_dir = args.output_dir if args.skip_subpath else f'{args.output_dir}/{nomenclature_name}'
    call_gather = f'python {args.dante_repo}/gather_result_files.py {output_dir}'
    call_motif = f'python {args.dante_repo}/create_motif_report.py {output_dir}'
    call_visual = f'python {args.dante_repo}/create_visual_report.py {output_dir}'
    run_command(call_gather, not args.ignore_errors, args.skip_call)
    run_command(call_motif, not args.ignore_errors, args.skip_call)
    run_command(call_visual, not args.ignore_errors, args.skip_call)

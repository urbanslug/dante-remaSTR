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
                                     description='Python program to collect Dante result files')

    # add arguments
    parser.add_argument('input_dir', help='Path to directory with Dante reports')
    parser.add_argument('output_dir', help='Path to directory where the output will be stored. Default=<input_dir>/result_files', nargs='?',
                        default=None)

    args = parser.parse_args()
    if args.output_dir is None:
        args.output_dir = args.input_dir + '/result_files'

    return args


if __name__ == '__main__':
    # load arguments
    args = load_arguments()

    # create the output dir if needed
    os.makedirs(args.output_dir, exist_ok=True)

    # copy results
    for report in glob.glob(f'{args.input_dir}/*/report.html'):
        name = report.split('/')[-2]
        to_call = f'cp {report} {args.output_dir}/{name}.html'
        subprocess.call(to_call, shell=True)

    # copy alignments
    for alignment in glob.glob(f'{args.input_dir}/*/*/alignments.html'):
        name = alignment.split('/')[-3]
        motif = alignment.split('/')[-2].replace('(', '').replace(')', '').replace(' ', '_')
        alignment = alignment.replace('(', r'\(').replace(')', r'\)').replace(' ', r'\ ')

        os.makedirs(os.path.join(args.output_dir, name), exist_ok=True)

        to_call = f'cp {alignment} {args.output_dir}/{name}/{motif}.html'
        subprocess.call(to_call, shell=True)

        if not os.path.exists(f'{args.output_dir}/msa.min.gz.js'):
            to_call = f'cp {os.path.dirname(os.path.dirname(alignment))}/msa.min.gz.js {args.output_dir}/msa.min.gz.js'
            subprocess.call(to_call, shell=True)

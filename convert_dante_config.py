import os
import argparse
import textwrap

import yaml
import dotmap


def load_arguments() -> argparse.Namespace:
    """
    Loads and parses arguments.
    :return: args - parsed arguments
    """
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent("""Converts old-style Dante .yaml configuration files to new-style pipeline
                                     and configuration for a quick transition from old-style to new."""))

    # add arguments
    parser.add_argument('yaml', help='Path to the old-style .yaml configuration file')
    parser.add_argument('output_dir', help='Path to directory where the output will be stored. Default=./new_style_config/', nargs='?',
                        default='new_style_config')
    parser.add_argument('--dante-repo', help='Path to Dante repository. Default=use this location', default=None)
    parser.add_argument('--reference-fasta', help='Path to reference fasta. Default=/data/genome/human/grch38_decoy/grch38_decoy.fa',
                        default='/data/genome/human/grch38_decoy/grch38_decoy.fa')

    args = parser.parse_args()

    # if not specified use location of this file
    if args.dante_repo is None:
        args.dante_repo = os.path.dirname(os.path.realpath(__file__))

    return args


def convert_module(module: str) -> str:
    """
    Converts an old-style motif report module into a new-style one.
    :param module: module string (old-style)
    :return: module string (new-style)
    """
    parts = module.split('-')
    if len(parts) == 1:
        return f'{parts[0]}[1]'
    else:
        return f'{parts[1]}[{int(parts[0])}]'


if __name__ == '__main__':
    # load arguments
    args = load_arguments()

    # read the yaml file
    with open(args.yaml, 'r') as f:
        old_yaml = dotmap.DotMap(yaml.safe_load(f))

    # create the output dir if needed
    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    # generate parameters and convert to dante call
    bam_basename = os.path.basename(old_yaml.inputs[0].path).split('.')[0]
    nomenclature_file = f'{args.output_dir}/nomenclature_{bam_basename[:15]}.txt'
    remastr_call = (f'{args.dante_repo}/remastr/target/release/remastr -f {args.reference_fasta} -n {nomenclature_file} '
                    f'-b {old_yaml.inputs[0].path} -o {args.output_dir}/{bam_basename}.tsv')
    options = (f'--param-file {old_yaml.allcall.param_file} '
               f'-p {old_yaml.general.cpu} ')
    if 'bases' in old_yaml.motifs[0].postfilter[0]:
        options += (f'--min-flank-len {old_yaml.motifs[0].postfilter[0].bases.split(",")[0]} '
                    f'--min-rep-len {old_yaml.motifs[0].postfilter[0].bases.split(",")[1]} ')
    if 'repetitions' in old_yaml.motifs[0].postfilter[0]:
        options += f'--min-rep-cnt {old_yaml.motifs[0].postfilter[0].repetitions.split(",")[1]} '
    if 'max_errors' in old_yaml.motifs[0].postfilter[0]:
        options += f'--max-rel-error {old_yaml.motifs[0].postfilter[0].max_errors} '
    dante_call = f'{args.dante_repo}/dante_remastr.py {options} < {args.output_dir}/{bam_basename}.tsv > {args.output_dir}/{bam_basename}_res.tsv'

    # write the generated pipeline
    with open(f'{args.output_dir}/pipeline_{bam_basename[:15]}.sh', 'w') as pipeline_file:
        pipeline_file.write(f'{remastr_call}\n')
        pipeline_file.write(f'{dante_call}\n')

    # convert the yaml file to nomenclatures
    with open(nomenclature_file, 'w') as nomenclature_file:
        # generate modules/nomenclatures
        for motif in old_yaml.motifs:
            # generate string
            module_str = ''.join([convert_module(module.seq) for module in motif.modules[1:-1]])

            # check if chromosome_version is defined
            # if 'chromosome_version' in motif and False:
            #    nomenclature = f'{motif.chromosome_version}:g.{motif.ref_start}_{motif.ref_end}{module_str}'
            # else:
            nomenclature = f'{motif.chromosome}:g.{motif.ref_start}_{motif.ref_end}{module_str}'

            # write it to file
            nomenclature_file.write(f'{nomenclature}\n')

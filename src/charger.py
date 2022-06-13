"""
Wrapper for Fernanda charger scripts
"""
import argparse
import os
import re
import logging
import subprocess
import stat


logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.INFO)

parser = argparse.ArgumentParser()

parser.add_argument('vcf', type=str,
    help='Filepath of inputs vcf.gz')

parser.add_argument('--inheritance-gene-list', type=str,
    help='Inheritance gene list')

parser.add_argument('--pp2-gene-list', type=str,
    help='P2 gene list')

parser.add_argument('--pathogenic-variants', type=str,
    help='Pathogenic variants vcf')

parser.add_argument('--hotspot3d-clusters', type=str, default='logs',
    help='HotSpot3D clusters file')

parser.add_argument('--clinvar-alleles', type=str,
    help='Clinvar allele file')

parser.add_argument('--rare-threshold', type=float, default=.0005,
    help='rare threshold')

parser.add_argument('--out-dir', type=str,
    help='output dir')

parser.add_argument('--sample', type=str, default='sample',
    help='sample id')

parser.add_argument('--format-vcf-script', type=str,
    help='TinJasmine vcf formatting script')

parser.add_argument('--post-charger-script', type=str,
    help='Post charger fixup script')

parser.add_argument('--filter-charger-script', type=str,
    help='Filters variants from charger output')

args = parser.parse_args()


def prepare_tinjasmine_vcf(format_vcf_script, vcf, out_dir):
    pieces = [
        'python {f}'.format(f=format_vcf_script),
        '-i {i}'.format(i=vcf),
        '-O {o}'.format(o=out_dir),
    ]
    cmd = ' '.join(pieces)
    return cmd


def execute_charger(preprocessed_vcf, output_fp, inheritance_gene_list,
                pp2_gene_list, pathogenic_variants, hotspot3d_clusters,
                clinvar_alleles, rare_threshold=.0005):
    pieces = [
        'charger --include-vcf-details',
        '-f {f}'.format(f=preprocessed_vcf),
        '-o {o}'.format(o=output_fp),
        '-O -D',
        '--inheritanceGeneList {i}'.format(i=inheritance_gene_list),
        '--PP2GeneList {p}'.format(p=pp2_gene_list),
        '-z {z}'.format(z=pathogenic_variants),
        '-H {h}'.format(h=hotspot3d_clusters),
        '-l',
        '--mac-clinvar-tsv {m}'.format(m=clinvar_alleles),
        '--rare-threshold {t}'.format(t=rare_threshold)
    ]
    cmd = ' '.join(pieces)
    return cmd


def post_charger(post_charger_script, preprocessed_vcf, charger_out,
                 sample, out_dir):
    pieces = [
        'python {f}'.format(f=post_charger_script),
        '-i {i}'.format(i=preprocessed_vcf),
        '-c {c}'.format(c=charger_out),
        '-s {s}'.format(s=sample),
        '-O {o}'.format(o=out_dir),
    ]
    cmd = ' '.join(pieces)
    return cmd


def filter_charger(filter_charger_script, post_charger_out, sample,
                   out_dir, rare_threshold=.0005):
    pieces = [
        'python {f}'.format(f=filter_charger_script),
        '-c {c}'.format(c=post_charger_out),
        '-a {a}'.format(a=rare_threshold),
        '-o {s}'.format(s=sample),
        '-O {o}'.format(o=out_dir),
    ]
    cmd = ' '.join(pieces)
    return cmd


def run_charger(
        args, preprocess_out, charger_out, post_charger_out, filter_charger_out):
    root = args.vcf.split('/')[-1].split('.vcf')[0]

    logging.info('step 1: preprocess vcf')
    cmd = prepare_tinjasmine_vcf(args.format_vcf_script, args.vcf, preprocess_out)
    logging.info('executing command: {c}'.format(c=cmd))
    output = subprocess.check_output(cmd, shell=True)
    logging.info('step output: {o}'.format(o=output))

    logging.info('step 2: run charger')
    preprocessed_vcf = os.path.join(preprocess_out, '{r}.infoFixed.vcf'.format(r=root))
    charger_out_fp = os.path.join(charger_out, '{r}.charged.tsv'.format(r=root))
    cmd = execute_charger(
        preprocessed_vcf, charger_out_fp, args.inheritance_gene_list,
        args.pp2_gene_list, args.pathogenic_variants, args.hotspot3d_clusters,
        args.clinvar_alleles, rare_threshold=args.rare_threshold)
    logging.info('executing command: {c}'.format(c=cmd))
    output = subprocess.check_output(cmd, shell=True)
    logging.info('step output: {o}'.format(o=output))

    logging.info('step 3: postprocess charger')
    cmd = post_charger(
        args.post_charger_script, preprocessed_vcf, charger_out_fp,
        args.sample, post_charger_out)
    logging.info('executing command: {c}'.format(c=cmd))
    output = subprocess.check_output(cmd, shell=True)
    logging.info('step output: {o}'.format(o=output))

    logging.info('step 4: filter charger')
    post_charger_out_fp = os.path.join(
        post_charger_out, '{s}.charged2vcf.tsv'.format(s=args.sample))
    filter_charger(
        args.filter_charger_script, post_charger_out_fp, args.sample,
        filter_charger_out, rare_threshold=args.rare_threshold)
    logging.info('executing command: {c}'.format(c=cmd))
    output = subprocess.check_output(cmd, shell=True)
    logging.info('step output: {o}'.format(o=output))

    fp = os.path.join(
        filter_charger_out, '{s}.charged2vcf.filtered.tsv'.format(s=args.sample))
    logging.info('output written to {f}'.format(f=fp))


def main():
    preprocess_out = os.path.join(args.out_dir, '1.preprocess_vcf')
    charger_out = os.path.join(args.out_dir, '2.charger')
    post_charger_out = os.path.join(args.out_dir, '3.post_charger')
    filter_charger_out = os.path.join(args.out_dir, '4.filter_charger')

    for dir in [preprocess_out, charger_out, post_charger_out, filter_charger_out]:
        if not os.path.exists(dir):
            os.makedirs(dir)

    run_charger(
        args, preprocess_out, charger_out, post_charger_out, filter_charger_out)


if __name__ == '__main__':
    main()

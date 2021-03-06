#!/usr/bin/env python3

from pathlib import Path
import tempfile


#############
# FUNCTIONS #
#############

def busco_resolver(wildcards):
    return({
        'fasta': busco_inputs[wildcards.name]
        })


def resolve_path(x):
    return Path(x).resolve()


###########
# GLOBALS #
###########

asw_assembly = 'data/curated.fasta' # this is the purge_haplotigs output
raw_read_dir = 'data/rnaseq_reads'

# containers
funannotate = ('shub://TomHarrop/funannotate-singularity:funannotate_1.7.4'
               '@c5e7e94e1830f825ad4dabc1af29413c65c3fd13')
funannotate_conda = ('shub://TomHarrop/funannotate-singularity:'
                     'funannotate-conda_1.7.4')
busco = 'shub://TomHarrop/singularity-containers:busco_3.0.2'

########
# MAIN #
########

busco_inputs = {
    'predict_transcriptome':
    'output/020_funannotate/annotate_results/ASW.mrna-transcripts.fa',
    'original_transcriptome':
    'data/Trinity.fasta'
}

# get the trinity input
all_rnaseq_samples = glob_wildcards(
    f'{raw_read_dir}/{{sample}}_r1.fq.gz').sample


#########
# RULES #
#########

rule target:
    input:
        'output/020_funannotate/annotate_results/ASW.annotations.txt',
        expand('output/099_busco/run_{name}/full_table_{name}.tsv',
               name=list(busco_inputs.keys()))

# annotate
rule funannotate_annotate:
    input:
        'output/020_funannotate/predict_results/ASW.gff3',
        egg = 'output/030_eggnog/eggnog.emapper.annotations',
        ipr = 'output/040_interproscan/ASW.proteins.fa.xml',
        db = 'data/fundb_20200227',
    output:
        'output/020_funannotate/annotate_results/ASW.annotations.txt',
        'output/020_funannotate/annotate_results/ASW.mrna-transcripts.fa'
    params:
        predict_dir = resolve_path('output/020_funannotate/predict_results'),
        db = lambda wildcards, input: resolve_path(input.db),
        wd = resolve_path('output/020_funannotate'),
        egg = lambda wildcards, input: resolve_path(input.egg),
        ipr = lambda wildcards, input: resolve_path(input.ipr),
    log:
        'output/logs/funannotate_annotate.log'
    threads:
        workflow.cores
    singularity:
        funannotate_conda
    shell:
        'bash -c \''
        'funannotate annotate '
        '--input {params.predict_dir} '
        '--out {params.wd} '
        '-s ASW '
        '--eggnog {params.egg} '
        '--iprscan {params.ipr} '
        '--busco_db endopterygota '
        '-d {params.db} '
        '--cpus {threads} '
        '\' &> {log}'


# run iprscan manually
rule iprscan:
    input:
        'output/020_funannotate/predict_results/ASW.proteins.fa'
    output:
        'output/040_interproscan/ASW.proteins.fa.xml'
    log:
        'output/logs/iprscan.log'
    params:
        wd = 'output/040_interproscan',
        tmpdir = tempfile.mkdtemp()
    threads:
        workflow.cores
    # can't put on shub, see github.com/tomharrop/funannotate-singularity
    singularity:
        'interproscan_5.44-79.0.sif'
    shell:
        'interproscan.sh '
        '-dp '
        '-i {input} '
        '--tempdir {params.tmpdir} '
        '--output-dir {params.wd} '
        '--cpu {threads} '
        '&> {log}'

# run eggnog-mapper manually
rule eggnog_mapper:
    input:
        fa = 'output/020_funannotate/predict_results/ASW.proteins.fa',
        db = 'data/eggnog_proteins.dmnd'
    output:
        'output/030_eggnog/eggnog.emapper.annotations'
    params:
        fa = lambda wildcards, input: resolve_path(input.fa),
        wd = resolve_path('output/030_eggnog'),
        db = lambda wildcards, input: resolve_path(input.db),
        db_path = lambda wildcards, input: resolve_path(input.db).parent
    log:
        resolve_path('output/logs/eggnog_mapper.log')
    threads:
        workflow.cores
    singularity:
        funannotate_conda
    shell:
        'bash -c \''
        'cd {params.wd} || exit 1 ; '
        'emapper.py '
        '-m diamond '
        '-i {params.fa} '
        '-o eggnog '
        '--dmnd_db {params.db} '
        '--data_dir {params.db_path} '
        '--cpu {threads} '
        '\' &> {log}'


# update models - not working in 1.7.4
# Error: input file 'long-reads.mapped.fasta' not found 
# rule funannotate_update:
#     input:
#         'output/020_funannotate/predict_results/ASW.mrna-transcripts.fa',
#         fasta = ('output/010_prepare/repeatmasker/'
#                  'asw-cleaned_sorted.fasta.masked')
#     output:
#         'output/020_funannotate/update_results/idk'
#     params:
#         wd = resolve_path('output/020_funannotate')
#     log:
#         'output/logs/funannotate_update.log'
#     threads:
#         workflow.cores
#     singularity:
#         funannotate_conda
#     shell:
#         'bash -c \''
#         'funannotate update '
#         '-i {params.wd} '
#         '--cpus {threads} '
#         '\' &> {log}'

# try to predict
rule funannotate_predict:
    input:
        'output/020_funannotate/training/funannotate_train.transcripts.gff3',
        fasta = ('output/010_prepare/repeatmasker/'
                 'asw-cleaned_sorted.fasta.masked'),
        db = 'data/fundb_20200227'
    output:
        'output/020_funannotate/predict_results/ASW.gff3',
        'output/020_funannotate/predict_results/ASW.mrna-transcripts.fa',
        'output/020_funannotate/predict_results/ASW.proteins.fa'
    params:
        fasta = lambda wildcards, input: resolve_path(input.fasta),
        db = lambda wildcards, input: resolve_path(input.db),
        wd = resolve_path('output/020_funannotate')
    log:
        'output/logs/funannotate_predict.log'
    threads:
        workflow.cores
    singularity:
        funannotate_conda
    shell:
        'bash -c \''
        'cp /genemark/gm_key_64 ${{HOME}}/.gm_key ; '
        'funannotate predict '
        '-i {params.fasta} '
        '-s ASW '
        # '--transcript_evidence {input.trinity} '
        '-o {params.wd} '
        '-d {params.db} '
        '--cpus {threads} '
        '--augustus_species lbonariensis '
        '--optimize_augustus '
        '--busco_seed_species tribolium2012 '
        '--busco_db endopterygota '
        '--organism other '
        '--repeats2evm '
        '--max_intronlen 50000 '
        '\' &> {log}'

# run training algorithm
# DON'T RUN WITH CONTAINALL
rule funannotate_train:
    input:
        fasta = 'output/010_prepare/repeatmasker/asw-cleaned_sorted.fasta.masked',
        left = 'output/020_funannotate/rnaseq_r1.fq.gz',
        right = 'output/020_funannotate/rnaseq_r2.fq.gz',
    output:
        'output/020_funannotate/training/funannotate_train.transcripts.gff3',
    params:
        fasta = lambda wildcards, input: resolve_path(input.fasta),
        wd = resolve_path('output/020_funannotate'),
    log:
        'output/logs/funannotate_train.log'
    threads:
        workflow.cores
    singularity:
        funannotate_conda
    shell:
        'bash -c \''
        'cp /genemark/gm_key_64 ${{HOME}}/.gm_key ; '
        'funannotate train '
        '--input {params.fasta} '
        '--out {params.wd} '
        '--left {input.left} '
        '--right {input.right} '
        '--stranded RF '
        # '--trinity {params.wd}/trinity.fasta '
        # '--no_trimmomatic ' # disabling trimmomatic seems to disable normalising, bleuch
        '--max_intronlen 10000 '
        '--species ASW '
        '--cpus {threads} '
        '\' &> {log}'

rule combine_rnaseq_reads:
    input:
        left = expand(f'{raw_read_dir}/{{sample}}_r1.fq.gz',
                      sample=all_rnaseq_samples),
        right = expand(f'{raw_read_dir}/{{sample}}_r2.fq.gz',
                       sample=all_rnaseq_samples)
    output:
        left = temp('output/020_funannotate/rnaseq_r1.fq.gz'),
        right = temp('output/020_funannotate/rnaseq_r2.fq.gz')
    params:
        left = lambda wildcards, input: sorted(set(input.left)),
        right = lambda wildcards, input: sorted(set(input.right))
    singularity:
        funannotate
    shell:
        'cat {params.left} > {output.left} & '
        'cat {params.right} > {output.right} & '
        'wait'

# funannotate mask?
rule fa_mask:
    input:
        asw_assembly
    output:
        'output/010_prepare/repeatmasker/asw-cleaned_sorted.fasta.masked'
    log:
        'output/logs/fa_mask.log'
    threads:
        workflow.cores
    singularity:
        funannotate_conda
    shell:
        'bash -c \''
        'funannotate mask '
        '-i {input} '
        '-o {output} '
        '--cpus {threads} '
        '\' &> {log}'


# generic busco rule
rule busco:
    input:
        unpack(busco_resolver),
        lineage = 'data/endopterygota_odb9'
    output:
        'output/099_busco/run_{name}/full_table_{name}.tsv'
    log:
        resolve_path('output/logs/099_busco/{name}.log')
    params:
        wd = 'output/099_busco/',
        fasta = lambda wildcards, input: resolve_path(input.fasta),
        lineage = lambda wildcards, input: resolve_path(input.lineage),
        tmpdir = tempfile.mkdtemp()
    threads:
        workflow.cores
    priority:
        1
    singularity:
        busco
    shell:
        'cd {params.wd} || exit 1 ; '
        'run_BUSCO.py '
        '--force '
        '--tmp_path {params.tmpdir} '
        '--in {params.fasta} '
        '--out {wildcards.name} '
        '--lineage {params.lineage} '
        '--cpu {threads} '
        '--mode transcriptome '
        '&> {log}'


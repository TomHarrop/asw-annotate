#!/usr/bin/env python3

import multiprocessing
import pathlib2
import tempfile


#############
# FUNCTIONS #
#############

def busco_resolver(wildcards):
    return({
        'fasta': busco_inputs[wildcards.name]
        })


def resolve_path(x):
    return str(pathlib2.Path(x).resolve())


###########
# GLOBALS #
###########

asw_assembly = 'data/curated.fasta' # this is the purge_haplotigs output
raw_read_dir = 'data/rnaseq_reads'

# containers
# funannotate = ('shub://TomHarrop/funannotate-singularity:funannotate_1.6.0'
#                '@5d0496b71cc229fc31cf06953737f9c4038ee51a')
funannotate = 'shub://TomHarrop/funannotate-singularity:funannotate_8e2e0a1'
busco = 'shub://TomHarrop/singularity-containers:busco_3.0.2'

########
# MAIN #
########

busco_inputs = {
    'predict_transcriptome':
    'output/020_funannotate/predict_results/ASW.mrna-transcripts.fa',
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
        ('output/020_funannotate/predict_results/ASW.gff3'),
        expand('output/099_busco/run_{name}/full_table_{name}.tsv',
               name=list(busco_inputs.keys()))


# try to predict
rule funannotate_predict:
    input:
        'output/020_funannotate/training/funannotate_train.transcripts.gff3',
        fasta = ('output/010_prepare/repeatmasker/'
                 'asw-cleaned_sorted.fasta.masked'),
        db = 'data/fundb_20190729',
        trinity = 'data/Trinity.fasta'
    output:
        'output/020_funannotate/predict_results/ASW.gff3',
        'output/020_funannotate/predict_results/ASW.mrna-transcripts.fa'
    params:
        fasta = lambda wildcards, input: resolve_path(input.fasta),
        db = lambda wildcards, input: resolve_path(input.db),
        wd = resolve_path('output/020_funannotate')
    log:
        'output/logs/funannotate_predict.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        funannotate
    shell:
        'cp /genemark/gm_key_64 ${{HOME}}/.gm_key ; '
        'funannotate predict '
        '-i {params.fasta} '
        '-s ASW '
        '--transcript_evidence {input.trinity} '
        '-o {params.wd} '
        '-d {params.db} '
        '--cpus {threads} '
        '--augustus_species lbonariensis '
        '--busco_seed_species tribolium2012 '
        '--busco_db endopterygota '
        '--organism other '
        '--repeats2evm '
        '--max_intronlen 10000 '
        '&> {log}'

# run training algorithm
rule funannotate_train:
    input:
        fasta = ('output/010_prepare/repeatmasker/'
                 'asw-cleaned_sorted.fasta.masked'),
        left = 'output/020_funannotate/rnaseq_r1.fq.gz',
        right = 'output/020_funannotate/rnaseq_r2.fq.gz',
        trinity = 'data/Trinity.fasta'
    output:
        'output/020_funannotate/training/funannotate_train.transcripts.gff3',
    params:
        fasta = lambda wildcards, input: resolve_path(input.fasta),
        wd = resolve_path('output/020_funannotate'),
    log:
        'output/logs/funannotate_train.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        funannotate
    shell:
        'cp /genemark/gm_key_64 ${{HOME}}/.gm_key ; '
        'funannotate train '
        '--input {params.fasta} '
        '--out {params.wd} '
        '--left {input.left} '
        '--right {input.right} '
        '--stranded RF '
        '--trinity {input.trinity} '
        '--no_trimmomatic ' # disabling trimmomatic seems to disable normalising, bleuch
        '--max_intronlen 10000 '
        '--species ASW '
        '--cpus {threads} '
        '&> {log}'

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

# manually mask the assembly
rule funnanotate_mask:
    input:
        cons = 'output/010_prepare/repeatmasker/consensi.fa',
        fasta = 'output/010_prepare/repeatmasker/asw-cleaned_sorted.fasta'
    output:
        ('output/010_prepare/repeatmasker/'
         'asw-cleaned_sorted.fasta.masked')
    params:
        wd = resolve_path('output/010_prepare/repeatmasker'),
        lib = lambda wildcards, input: resolve_path(input.cons),
        fasta = lambda wildcards, input: resolve_path(input.fasta)
    log:
        resolve_path('output/logs/funnanotate_mask.log')
    threads:
        multiprocessing.cpu_count()
    singularity:
        funannotate
    shell:
        'cd {params.wd} || exit 1 ; '
        'RepeatMasker '
        '-engine ncbi '
        '-pa {threads} '
        '-lib {params.lib} '
        '-dir {params.wd} '
        '-gccalc -xsmall -gff -html '
        '{params.fasta} '
        '&> {log}'


# I don't think this can be done without the GIRINST library
# instead, annotate with dfam? https://www.dfam.org/help/tools
# rule funnanotate_mask_classify:
#     input:
#         'output/010_prepare/repeatmasker/families.stk',
#         'output/010_prepare/repeatmasker/consensi.fa'
#     output:
#         'output/010_prepare/repeatmasker/consensi.fa.classified'
#     params:
#         wd = resolve_path('output/010_prepare/repeatmasker'),
#     log:
#         resolve_path('output/logs/funnanotate_mask-classify.log')
#     singularity:
#         "funannotate_1.6.0-rminstall.sif"
#     shell:
#         'cd {params.wd} || exit 1 ; '
#         'RepeatClassifier '
#         '-engine ncbi '
#         '-consensi consensi.fa '
#         '-stockholm families.stk '
#         '&> {log}'


rule funnanotate_mask_model:
    input:
        'output/010_prepare/repeatmasker/ASW.translation'
    output:
        'output/010_prepare/repeatmasker/families.stk',
        'output/010_prepare/repeatmasker/consensi.fa'
    params:
        wd = resolve_path('output/010_prepare/repeatmasker'),
    log:
        resolve_path('output/logs/funnanotate_mask-model.log')
    threads:
        multiprocessing.cpu_count()
    singularity:
        funannotate
    shell:
        'cd {params.wd} || exit 1 ; '
        'RepeatModeler '
        '-database ASW '
        '-engine ncbi '
        '-pa {threads} '
        '-dir {params.wd} '
        # '-recoverDir {params.wd} '
        '&> {log}'

rule funnanotate_mask_build:
    input:
        fasta = 'output/010_prepare/asw-cleaned_sorted.fasta'
    output:
        'output/010_prepare/repeatmasker/ASW.translation',
        'output/010_prepare/repeatmasker/asw-cleaned_sorted.fasta'
    params:
        wd = resolve_path('output/010_prepare/repeatmasker'),
        fasta = lambda wildcards, input: resolve_path(input.fasta)
    log:
        resolve_path('output/logs/funnanotate_mask-build.log')
    threads:
        multiprocessing.cpu_count()
    singularity:
        funannotate
    shell:
        'cd {params.wd} || exit 1 ; '
        'cp {params.fasta} . ; '
        'BuildDatabase '
        '-name ASW '
        '-engine ncbi '
        '-dir {params.wd} '
        '&> {log} '

rule funnanotate_sort:
    input:
        'output/010_prepare/asw-cleaned.fasta'
    output:
        'output/010_prepare/asw-cleaned_sorted.fasta'
    log:
        'output/logs/funnanotate_sort.log'
    singularity:
        funannotate
    shell:
        'funannotate sort '
        '--input {input} '
        '--out {output} '
        '&> {log}'


rule funnanotate_clean:
    input:
        asw_assembly
    output:
        'output/010_prepare/asw-cleaned.fasta'
    log:
        'output/logs/funnanotate_clean.log'
    singularity:
        funannotate
    shell:
        'funannotate clean '
        '--input {input} '
        '--out {output} '
        '&> {log}'


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
        multiprocessing.cpu_count()
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


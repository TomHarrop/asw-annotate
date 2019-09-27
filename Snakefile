#!/usr/bin/env python3

import multiprocessing
import pathlib2


#############
# FUNCTIONS #
#############

def resolve_path(x):
    return str(pathlib2.Path(x).resolve())


###########
# GLOBALS #
###########

asw_assembly = 'data/flye_denovo_full.racon.fasta'

# containers
funannotate = 'shub://TomHarrop/funannotate-singularity:funannotate_1.6.0'

#########
# RULES #
#########

rule target:
    input:
        ('output/030_funannotate')

# try to predict
rule funannotate_predict:
    input:
        fasta = ('output/010_prepare/repeatmasker/'
                 'asw-cleaned_sorted.fasta.masked'),
        db = 'data/fundb_20190729',
        trinity = 'data/Trinity.fasta'
    output:
        directory('output/030_funannotate')
    log:
        'output/logs/funannotate_predict.log'
    threads:
        multiprocessing.cpu_count()
    singularity:
        'funannotate_1.6.0-ncbitools.sif'
    shell:
        'funannotate predict '
        '-i {input.fasta} '
        '-s ASW '
        '--transcript_evidence {input.trinity} '
        '-o {output} '
        '-d {input.db} '
        '--cpus {threads} '
        '&> {log}'


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
        '-recoverDir {params.wd} '
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

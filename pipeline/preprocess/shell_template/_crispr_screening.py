import logging
import multiprocessing
import math
from pipeline.common import path_check

cpu = math.floor(multiprocessing.cpu_count() * 0.8)


def mageck_count(samples, kwargs):
    """ Template of mageck count
    :param samples:
    :param suffix:
    :param library:
    :param sgrna_len: Length of sgRNA
    :param trim: Length for trim-5
    :param control_sgrna: Gene name of control sgRNA
    :param sample_labels: rename sample
    :return: command
    """
    # check directory
    path_check('mageck_count')

    # args parse
    try:
        suffix = kwargs['suffix']
        library = kwargs['library']
    except:
        lost_arg = list({'suffix', 'library'} - set(kwargs.keys()))
        logging.error('Template required parameter missing :' + ' '.join(lost_arg))
        return 'Template required parameter missing'
    sgrna_len = kwargs.get('sgrna_len')
    control_sgrna = kwargs.get('control_sgrna')
    control_gene = kwargs.get('control_gene')
    sample_labels = kwargs.get('sample_labels')

    # process use quality controled fastq
    filelist_1 = []
    filelist_2 = []
    sample_list = []
    for sample in samples:
        if type(sample) == str:
            filelist_1.append(sample + '_qc' + suffix[0])
            filelist_2.append(sample + '_qc' + suffix[1])
            sample_list.append(sample)
        elif type(sample) == tuple:
            sample_list.append(sample[0])
            rep_filelist_1 = []
            rep_filelist_2 = []
            for rep in sample:
                rep_filelist_1.append(rep + '_qc' + suffix[0])
                rep_filelist_2.append(rep + '_qc' + suffix[1])
            filelist_1.append(','.join(rep_filelist_1))
            filelist_2.append(','.join(rep_filelist_2))

    if not sample_labels:
        sample_labels = sample_list
    fastq_1 = '--fastq fastqs/' + ' fastqs/'.join(filelist_1)
    fastq_2 = '--fastq-2 fastqs/' + ' fastqs/'.join(filelist_2)
    sample_labels = '--sample-label ' + ','.join(sample_labels)
    library = '-l ' + library
    out = '-n mageck'
    if sgrna_len:
        sgrna_len = '--sgrna-len ' + str(sgrna_len)
    else:
        sgrna_len = ''
    if control_sgrna:
        control_sgrna = '--control-sgrna ' + control_sgrna
    else:
        control_sgrna = ''
    if control_gene:
        control_gene = '--control-gene ' + control_gene
    else:
        control_gene = ''
    command = 'mageck count %s %s %s %s %s %s %s %s --pdf-report --unmapped-to-file --norm-method control > mageck.log' % (
        fastq_1, fastq_2, sample_labels, library, sgrna_len, control_sgrna, control_gene, out)

    logging.debug(filelist_1)
    logging.debug(filelist_2)
    logging.debug(sample_list)
    command_list = [command]

    return command_list


def mageck_vispr(samples, **kwargs):
    """mageck-vispr config.yaml
    :param samples:
    :param suffix:
    :param library:
    :param sgrna_len: Length of sgRNA
    :param trim: Length for trim-5
    :param control_sgrna: Gene name of control sgRNA
    :param sample_labels: rename sample
    :return: command
    """
    # check directory
    path_check('mageck_vispr')

    # args parse
    try:
        suffix = kwargs['suffix']
        library = kwargs['library']
    except:
        lost_arg = list({'suffix', 'library'} - set(kwargs.keys()))
        logging.error('Template required parameter missing :' + ' '.join(lost_arg))
        return 'Template required parameter missing'

    sgrna_len = kwargs.get('sgrna_len')
    if not sgrna_len:
        sgrna_len = 'AUTO'
    else:
        sgrna_len = '"' + str(sgrna_len) + '"'
    trim = kwargs.get('trim')
    if not trim:
        trim = 'AUTO'
    else:
        trim = '"' + str(trim) + '"'
    control_sgrna = kwargs.get('control_sgrna')
    if control_sgrna:
        control_sgrna = """norm_method: control
control_sgrna: ../""" + control_sgrna
    else:
        control_sgrna = """# norm_method: control
        # control_sgrna: lib/hg19_library_1.aavs1.txt"""
    sample_labels = kwargs.get('sample_labels')
    if not sample_labels:
        sample_labels = samples
    filelist = []
    for i in range(len(samples)):
        sample = samples[i]
        sample_label = sample_labels[i]
        sample_file = """
    %s:
        - qc_fastqs/%s_qc%s
        - qc_fastqs/%s_qc%s""" % (sample_label, sample, suffix[0], sample, suffix[1])
        filelist.append(sample_file)

    config = """# General configuration:

# Path to library design file (csv or tab-separated txt format, columns: id, sequence, gene)
library: ../%s
# Species to use for linkouts in VISPR (e.g. mus_musculus, homo_sapiens, ...)
species: homo_sapiens
# Genome assembly to use for linkouts in VISPR (e.g. hg19, hg38, mm9, mm10, ...)
assembly: hg38

# Configuration of knockout target display in VISPR
targets:
    # if screening genes, set this to true for proper linkouts to GeneMANIA and Ensembl in VISPR
    genes: true
    # file with genes to hide per default in VISPR (optional, one gene per line)
    #controls: ribosomal_genes.txt

# Configuration of sgRNAs
sgrnas:
    # estimate sgRNA knockout efficiency during EM-procedure of MAGeCK-MLE
    update-efficiency: false
    # trim the 5 prime end to get rid of barcode sequences in the reads
    # if a number (instead of AUTO) is specified, use quotes; for example:
    # trim-5: "0"
    trim-5: %s
    # specify the length of the sgRNAs (without PAM sequence)
    len: %s
    # sequencing adapter that shall be removed from reads before processing with MAGeCK (optional)
    #adapter: ACGGCTAGCTGA
    #
    # Use pre-computed sgrnas to annotate the library? By default it's false. 
    # Only certain assemblies (hg19, hg38, mm9, mm10) and certain sgRNA length (19, 20) are supported
    annotate-sgrna: false
    # Use pre-computed sgrna efficiency as an initial value of knockout efficiency?
    # Need to set annotate-sgrna to true as well. 
    # If you need the sgRNA efficiency to be updated, set update-efficiency to true
    annotate-sgrna-efficiency: false
    # instead of downloading the sgrna annotation library from bitbucket,
    # provide either the name of the file used for annotation,
    # or the folder name where MAGeCK-VISPR will search the corresponding annotation library from that folder
    #annotation-sgrna-file: /dev/null
    #annotation-sgrna-folder: /src/exome_scan




# ATTENTION: You should and only should choose one type of input file between "samples" and "counts".
# Configuration of samples (Cannot set "counts" as valid at the same time!)
samples:
    # The following sample information was inferred from the given FASTQ files.
    # Adjust it according to your needs (e.g. providing descriptive sample names and grouping replicates together).
    %s


# Provide paired fastq files if pair-end sequencing data is available.
# paired:
#     # provide a label and a paths to the paired fastq files for each sample
#     A: path/to/A_R2.fastq
#     B: path/to/B_R2.fastq
#     C:
#        - path/to/C.1_R2.fastq
#         - path/to/C.2_R2.fastq

# Specify whether report valid alignments per pair (count when both end reads mapped), when paired fastq files are provided .
countpair: false

# Instead of providing fastq files, you can also provide your own normalized (or unnormalized) count matrix (Cannot set "samples" as valid at the same time!).
# If you do not want MAGeCK to normalize the counts, make sure to set up norm_methd to none
# Support QC from count matrix
# counts: rawcount/rawcount.txt

# Provide mageck count --day0-label (optional). Multiple labels should be seperated with comma.
# day0label: plasmid


# Provide normalization methods (none, control, median) and a list of negative control sgrnas (if norm_method=control). 
# These parameters will be applied to mageck_count, mageck_rra and mageck_mle modules.
# Note that if this option is not specified a default median normalization will be used.
# norm_method: control
%s

# Run mle with multi-thread. Default thread number is 1.
# When this parameter is set (e.g., threads: 4), make sure to specify the cores in running snakemake (snakemake --cores 4)
threads: %s

# Provide a batch matrix if the samples need to be batch corrected (optional).
# The format should be as follows (tab-separated):
# sample          batch   covariate 1 ...
# Sample 1        1       0           ...
# Sample 2        2       0           ...
# Sample 3        2       1           ...
#
# The first column must refer to the samples defined above.
# batchmatrix: path/to/batchmatrix.txt
#

# copy number variation correction
# Setting correct_cnv to true will allow CNV correction for both RRA and MLE
#
# to overwrite this option, you can use snakemake command line:
# snakemake --config correct_cnv=true cnv_norm=target/cnv/file
correct_cnv: false
cnv_norm: /dev/null

# the following cnv_cell_line is only required when performing CNV correction in RRA
# it's not required for MLE
cnv_cell_line: HL60


# Additional parameters to run MLE or RRA
# Note that additional_mle_rra_parameter is no longer used since MLE and RRA should be specified separately.
# additional_mle_parameter: --max-sgrnapergene-permutation 30
# additional_rra_parameter: --remove-zero both --additional-rra-parameters "--max-sgrnapergene-permutation 200"

# Configuration of experiments.
# An experiment defines the comparison that shall be analyzed with MAGeCK.
# You can define as many experiments as you want.
# You can define either MAGeCK-RRA or MAGeCK-MLE experiments, but cannot define both in a single configuration file.
experiments:
    # provide a descriptive name for your experiment (it will show up in VISPR)
    "mle":
        # This is a MAGeCK-MLE experiment.
        # Here, users can either specify a day0label, or a design matrix file (see http://mageck.sourceforge.net for details).
        # if day0label is specified, use an empty file (like /dev/null) as a design matrix file.
        # Sample names in the design matrix must refer to the samples defined above.
        designmatrix: /dev/null
    #"rra":
        # This is a MAGeCK-RRA experiment.
        # You must specify treatment and control samples.
        # The sample names must refer to the samples defined above in the
        # samples section.
        #treatment:
        #    - A
        #control:
        #    - B
        #    - C""" % (library, trim, sgrna_len, ''.join(filelist), control_sgrna, cpu)
    with open('mageck-vispr/config.yaml', 'w') as f:
        f.write(config)

    command = """cd mageck-vispr
    
              """
    return config

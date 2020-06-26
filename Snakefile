from snakemake.utils import min_version

##### set minimum snakemake version #####
min_version("5.18.0")

configfile: "configs/config.yaml"
configfile: "configs/callers.yaml"

container: "docker://continuumio/miniconda3:4.8.2"


def check_config(value, default=False, place=config):
    """ return true if config value exists and is true """
    return place[value] if (value in place and place[value]) else default

def read_samples():
"""
    Function to get names and dna fastq paths from a sample file
    specified in the configuration. Input file is expected to have 3
    columns: <unique_sample_id> <fastq1_path> <fastq2_path> or
    <unique_sample_id> <paired_bam_path> <bed_path>. Modify this function
    as needed to provide a dictionary of sample_id keys and either a tuple
    of strings: (fastq1, fastq2) OR a single string: paired_bam
"""
f = open(config['sample_file'], "r")
samp_dict = {}
for line in f:
    words = line.strip().split("\t")
    if len(words) == 2:
        samp_dict[words[0]] = (words[1], "")
    elif len(words) == 3:
        samp_dict[words[0]] = (words[1], words[2])
    else:
        raise ValueError('Your samples_file is not formatted correctly. Make sure that it has the correct number of tab-separated columns for every row.')
return samp_dict
SAMP = read_samples()

# the user can change config['SAMP_NAMES'] here (or define it in the config
# file) to contain whichever sample names they'd like to run the pipeline on
if 'SAMP_NAMES' not in config or not config['SAMP_NAMES']:
    config['SAMP_NAMES'] = list(SAMP.keys())
else:
    # double check that the user isn't asking for samples they haven't provided
    user_samps = set(config['SAMP_NAMES'])
    config['SAMP_NAMES'] = list(set(SAMP.keys()).intersection(user_samps))
    if len(config['SAMP_NAMES']) != len(user_samps):
        warnings.warn("Not all of the samples requested have provided input. Proceeding with as many samples as is possible...")


rule all:
	input:
		expand(
            config['out']+"/classify/{sample}_{type}/final.vcf.gz",
            sample=config['SAMP_NAMES'],
            type=[i for i in ["snp", "indel"] if check_config(i+"_callers")]
        )

# an internal variable we use to tell the other subworkflows not to import their configs
config['imported'] = True

include: "rules/prepare.smk"

config['predict'] = []
config['data'] = {}
for samp in config['SAMP_NAMES']:
	for i in ['snp', 'indel']:
		if check_config(i+"_callers"):
			sample_name = samp + "_" + i
			config['predict'].append(sample_name)
			config['data'][sample_name] = {
				'path': config['out']+"/merged_"+i+"/"+samp+"/final.tsv.gz",
				'merged': config['out']+"/merged_"+i+"/"+samp+"/merge.tsv.gz",
				'model': config[i+"_model"]
			}
config['out'] += "/classify"

include: "rules/classify.smk"

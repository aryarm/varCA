from snakemake.utils import min_version

##### set minimum snakemake version #####
min_version("5.5.0")

configfile: "config.yaml"


rule all:
    input:
        expand(
            config['out']+"/{sample}/results.tsv.gz",
            sample=[
                s for s in config['predict']
            ]
        )

rule norm_nums:
    """
        convert pseudo-numerical column values (like those in
        scientific notation or percent format) to simple numerics
    """
    input: lambda wildcards: config['data'][wildcards.sample]['path']
    output: config['out']+"/{sample}/norm.tsv.gz"
    conda: "env.yml"
    shell:
        "zcat {input} | scripts/norm_nums.awk -F $'\\t' -v 'OFS=\\t' | "
        "gzip >{output}"

rule apply_filters:
    """ apply filtering on the data according to the filtering expressions """
    input: rules.norm_nums.output if 'norm_numerics' in config and \
        config['norm_numerics'] else rules.norm_nums.input
    params: lambda wildcards: "\t".join(config['data'][wildcards.sample]['filter'])
    output: config['out']+"/{sample}/filter.tsv.gz"
    conda: "env.yml"
    shell:
        "zcat {input} | scripts/filter.bash {params:q} | gzip > {output}"

rule annotate:
    """ create a table of annotations at each site """
    input: lambda wildcards: rules.apply_filters.output \
        if 'filter' in config['data'][wildcards.sample] and \
        config['data'][wildcards.sample]['filter'] else \
        rules.apply_filters.input,
    output: temp(config['out']+"/{sample}/annot.tsv.gz")
    conda: "env.yml"
    shell:
        "paste <(zcat {input} | cut -f -2) <(scripts/classify.bash {input} {config[label]}) | gzip > {output}"

rule fillna:
    """
        prepare the caller for use by the classifier by
        1) extracting the columns desired by the user
        2) filling NA values with the defaults provided
    """
    input:
        tsv = lambda wildcards: rules.apply_filters.output \
        if 'filter' in config['data'][wildcards.sample] and \
        config['data'][wildcards.sample]['filter'] else \
        rules.apply_filters.input
    params:
        na_vals = lambda wildcards: [
            j for i in config['data'][wildcards.sample]['na'].items()
            for j in i
        ],
    output: config['out']+"/{sample}/fillna.tsv.gz"
    conda: "env.yml"
    shell:
        "scripts/fillna.bash {input.tsv} {params.na_vals:q} | gzip >{output}"

rule add_truth:
    """ add true labels as the last column in the training data """
    input:
        tsv = rules.fillna.output,
        annot = rules.annotate.output
    params:
        truth = lambda wildcards: '^'+config['data'][wildcards.sample]['truth']+"~" if 'truth' in config['data'][wildcards.sample] and config['data'][wildcards.sample]['truth'] else ""
    output: config['out']+"/{sample}/fillna.truth.tsv.gz"
    conda: "env.yml"
    shell:
        "paste <(zcat {input.tsv}) <(zcat {input.annot} | scripts/get_cols.bash {params.truth:q}) | gzip > {output}"

rule train:
    """ train the classifier """
    input: rules.add_truth.output
    params:
        balance = int(config['balance']) if 'balance' in config else 0,
        importance = config['out']+"/{sample}/variable_importance.tsv"
    output: config['out']+"/{sample}/model.rda"
    conda: "env.yml"
    shell:
        "Rscript scripts/train_RF.R {input} {output} {params.balance} {params.importance}"

rule predict:
    """ predict variants using the classifier """
    input:
        model = lambda wildcards: expand(rules.train.output, sample=config['train']),
        predict = lambda wildcards: expand(rules.fillna.output, sample=wildcards.sample)
    conda: "env.yml"
    output: temp(config['out']+"/{sample}/predictions.tsv")
    shell:
        "Rscript scripts/predict_RF.R {input.predict} {input.model} {output}"

rule join_results:
    """ join the predictions with the annotations """
    input:
        predict = rules.predict.output,
        annot = rules.annotate.output
    output: config['out']+"/{sample}/results.tsv.gz"
    conda: "env.yml"
    shell:
        "paste <(zcat {input.annot}) {input.predict} | gzip > {output}"

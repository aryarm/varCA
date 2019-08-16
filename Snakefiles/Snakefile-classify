from snakemake.utils import min_version

##### set minimum snakemake version #####
min_version("5.5.0")

configfile: "config.yaml"


rule all:
    input:
        expand(
            config['out']+"/{sample}/results.tsv.gz",
            sample=config['predict']
        ) + expand(
            config['out']+"/{sample}/prc/results.png",
            sample=[
                s for s in config['predict']
                if 'truth' in config['data'][s] and config['data'][s]['truth']
            ] if 'prcols' in config and isinstance(config['prcols'], dict) else []
        )

rule norm_nums:
    """
        convert pseudo-numerical column values (like those in
        scientific notation or percent format) to simple numerics
    """
    input: lambda wildcards: config['data'][wildcards.sample]['path']
    output: config['out']+"/{sample}/norm.tsv.gz"
    conda: "envs/env.yml"
    shell:
        "zcat {input} | scripts/norm_nums.awk -F $'\\t' -v 'OFS=\\t' | "
        "gzip >{output}"

rule apply_filters:
    """ apply filtering on the data according to the filtering expressions """
    input: rules.norm_nums.output if 'norm_numerics' in config and \
        config['norm_numerics'] else rules.norm_nums.input
    params: lambda wildcards: "\t".join(config['data'][wildcards.sample]['filter'])
    output: config['out']+"/{sample}/filter.tsv.gz"
    conda: "envs/env.yml"
    shell:
        "zcat {input} | scripts/filter.bash {params:q} | gzip > {output}"

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
    output: temp(config['out']+"/{sample}/fillna.tsv.gz")
    conda: "envs/env.yml"
    shell:
        "scripts/fillna.bash {input.tsv} {params.na_vals:q} | gzip >{output}"

rule annotate:
    """ create a table of annotations at each site """
    input: lambda wildcards: rules.apply_filters.output \
        if 'filter' in config['data'][wildcards.sample] and \
        config['data'][wildcards.sample]['filter'] else \
        rules.apply_filters.input
    params:
        label = config['label'] if 'label' in config else "."
    output: temp(config['out']+"/{sample}/annot.tsv.gz")
    conda: "envs/env.yml"
    shell:
        "paste <(zcat {input} | cut -f -2) <(scripts/classify.bash {input} '{params.label:q}') | gzip > {output}"

rule add_truth:
    """
        Add true labels as the last columns in the training data
        Also ensure that if a label is available for this dataset, it appears
        as the very last column
    """
    input:
        tsv = rules.fillna.output,
        annot = rules.annotate.output
    params:
        truth = lambda wildcards: '^'+config['data'][wildcards.sample]['truth']+"~" if 'truth' in config['data'][wildcards.sample] and config['data'][wildcards.sample]['truth'] else ""
    output: config['out']+"/{sample}/prepared.tsv.gz"
    conda: "envs/env.yml"
    shell:
        "paste "
        "<(zcat {input.annot} | cut -f 3- | scripts/cgrep.bash - -v '{params.truth:q}') "
        "<(zcat {input.annot} | cut -f 3- | scripts/cgrep.bash - '{params.truth:q}') | "
        "sed 's/^\t//' | paste <(zcat {input.tsv}) - | gzip > {output}"

rule train:
    """ train the classifier """
    input: rules.add_truth.output
    params:
        balance = int(config['balance']) if 'balance' in config else 0
    output:
        model = config['out']+"/{sample}/model.rda",
        importance = config['out']+"/{sample}/variable_importance.tsv"
    conda: "envs/env.yml"
    shell:
        "Rscript scripts/train_RF.R {input} {output.model} {params.balance} {output.importance}"

rule predict:
    """ predict variants using the classifier """
    input:
        model = lambda wildcards: expand(rules.train.output.model, sample=config['train']),
        predict = lambda wildcards: expand(rules.add_truth.output, sample=wildcards.sample)
    conda: "envs/env.yml"
    output: temp(config['out']+"/{sample}/predictions.tsv")
    shell:
        "Rscript scripts/predict_RF.R {input.predict} {input.model} {output}"

rule results:
    """
        join the predictions with the annotations
        also prefix the colnames of our method before merging
    """
    input:
        predict = rules.predict.output,
        annot = rules.annotate.output
    params:
        truth = lambda wildcards: config['data'][wildcards.sample]['truth'] if 'truth' in config['data'][wildcards.sample] and config['data'][wildcards.sample]['truth'] else "",
        label = config['label'] if 'label' in config else "."
    output: config['out']+"/{sample}/results.tsv.gz"
    conda: "envs/env.yml"
    shell:
        "cat {input.predict} | paste <(zcat {input.annot}) "
        "<(read -r head && echo \"$head\" | tr '\\t' '\\n' | "
        "sed 's/response/CLASS:{params.label}/' | sed 's/^/breakca~/' | "
        "paste -s && cat) | gzip > {output}"

rule prc_pts:
    """ generate single point precision recall metrics """
    input:
        annot = rules.annotate.output,
        predicts = rules.results.output
    params:
        truth = lambda wildcards: config['data'][wildcards.sample]['truth'] if 'truth' in config['data'][wildcards.sample] and config['data'][wildcards.sample]['truth'] else ""
    output: config['out']+"/{sample}/prc/pts/{caller}.txt"
    conda: "envs/prc.yml"
    shell:
        "paste "
        "<(zcat {input.annot} | scripts/cgrep.bash - -E '{params.truth}' | sed 's/\./0/') "
        "<(zcat {input.predicts} | scripts/cgrep.bash - '{wildcards.caller}~CLASS:') | "
        "tail -n+2 | scripts/metrics.py -o {output}"


def sort_col(caller):
    if caller in config['prcols']:
        return config['prcols'][caller], False
    elif "*"+caller in config['prcols']:
        return config['prcols']["*"+caller], True
    else:
        return "", False


rule prc_curves:
    """ generate the points for a precision recall curve """
    input:
        annot = rules.annotate.output,
        predicts = lambda wildcards: rules.results.output if wildcards.caller == 'breakca' else rules.add_truth.input.tsv
    params:
        truth = lambda wildcards: config['data'][wildcards.sample]['truth'] if 'truth' in config['data'][wildcards.sample] and config['data'][wildcards.sample]['truth'] else "",
        predict_col = lambda wildcards: 'prob.1' if wildcards.caller == 'breakca' else sort_col(wildcards.caller)[0],
        flip = lambda wildcards: ["", "-f"][sort_col(wildcards.caller)[1]]
    output: config['out']+"/{sample}/prc/curves/{caller}.txt"
    conda: "envs/prc.yml"
    shell:
        "paste "
        "<(zcat {input.annot} | scripts/cgrep.bash - -E '{params.truth}' | sed 's/\./0/') "
        "<(zcat {input.predicts} | scripts/cgrep.bash - -F '{wildcards.caller}~{params.predict_col}') | "
        "tail -n+2 | scripts/statistics.py -o {output} {params.flip}"


def sort_cols(strict=False):
    return [
        caller[caller.startswith("*") and len("*"):]
        for caller in config['prcols'].keys()
        if not strict or config['prcols'][caller]
    ]


rule prc:
    """ create plot containing precision recall curves """
    input:
        pts = lambda wildcards: expand(
            rules.prc_pts.output, sample=wildcards.sample,
            caller=['breakca']+sort_cols()
        ),
        curves = lambda wildcards: expand(
            rules.prc_curves.output, sample=wildcards.sample,
            caller=['breakca']+sort_cols(True)
        )
    params:
        pts = lambda _, input: [k for j in zip(['--'+i+"_pt" for i in ['breakca']+sort_cols()], input.pts) for k in j],
        curves = lambda _, input: [k for j in zip(['--'+i for i in ['breakca']+sort_cols(True)], input.curves) for k in j]
    output: config['out']+"/{sample}/prc/results.png"
    conda: "envs/prc.yml"
    shell:
        "scripts/prc.py {output} {params.pts} {params.curves}"

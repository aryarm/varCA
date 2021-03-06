#!/usr/bin/env python3

import sys
import argparse
import matplotlib
import pandas as pd
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpt


parser = argparse.ArgumentParser()
parser.add_argument(
    "-o", "--out", default=sys.stdout, help="the filename to save the data to"
)
parser.add_argument(
    "-c", "--callers", default="", help="a comma separated list of the caller names to include"
)
parser.add_argument(
    "--caller", action='store_false', help="whether to compare importance to each individual feature (default if not specified) or each caller"
)
parser.add_argument(
    "table", nargs="?", default=sys.stdin,
    help="a two column (variable/importance) table of importances w/ a header"
)
args = parser.parse_args()


# import the data
df = pd.read_csv(
    args.table, sep="\t", header=0
).sort_values(by="importance", ascending=False)
df.columns = ['variable', 'importance']


# # assign a caller to each variable
# callers = ['gatk-indel', 'varscan-indel', 'vardict-indel', 'breakca', 'delly', 'pindel', 'illumina-manta', 'illumina-strelka', 'pg-indel']
# curr_callers = args.callers.split(',')
# for caller in curr_callers:
#     indel_version = caller[:-len('snp')]+"indel"
#     if caller.endswith('-snp') and indel_version not in callers:
#     if caller not in callers:
#         callers.append(caller)
#     elif caller.endswith('-snp') and snp_version in curr_callers:
#         callers.append(snp_version)
callers = args.callers.split(',')

# create a caller column in the df
df['caller'] = [-1] * len(df)
for v in range(len(df['variable'])):
    for i in range(len(callers)):
        c = callers[i].replace("-", ".")
        # if we've found a match!
        if df['variable'][v].startswith(c):
            # record the id of the caller in callers
            df.at[v, 'caller'] = i
            # also strip the caller name out of the variable name
            # print(df.at['variable', v])
            df.at[v, 'variable'] = df['variable'][v][len(c):].strip('.')
            break

# should we create a feature plot or a caller plot?
if args.caller:
    # pick a color for each variable based on its caller
    colors = plt.cm.Set1(df['caller'])

    # create the plot
    plot = df[['variable', 'importance']].plot.bar(x='variable', color=[colors])

    # create a legend label for each color
    patches = [
        mpt.Patch(color=i[0], label=callers[i[1]])
        for i in set(zip([tuple(color) for color in colors], df['caller']))
    ]
    plot.legend(handles=patches)
    plt.xlabel('Feature')
    plt.ylabel('Importance')
    plt.gcf().set_size_inches(13, 13)
else:
    # aggregate importance by caller
    df = df[['importance', 'caller']].groupby('caller').sum().sort_values('caller')
    df['callern'] = callers
    # create the plot
    plot = df.plot.bar(x='callern', color=[plt.cm.Set1(df.index)])
    plot.legend().remove()
    plt.xlabel('Caller')
    plt.ylabel('Importance Sum')

# save the plot
plt.savefig(args.out, bbox_inches='tight', pad_inches=0.1, set_dpi=1000)

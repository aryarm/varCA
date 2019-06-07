#!/usr/bin/env python3

import sys

output = open(sys.argv[1], "w+")
files = [open(file) for file in sys.argv[2:]]


def lines():
    stored_line = [None]*len(files)
    curr_chrom = 1
    while True:
        line = [None]*len(files)
        curr_min = float('Inf')
        at_min = []
        count_stops = 0
        for i in range(len(files)):
            if stored_line[i] is None:
                # try to capture the line
                # you might not be able to if we're at the end of the file
                try:
                    line_i = next(files[i]).split('\t', 2)
                    line_i[1] = int(line_i[1])
                except StopIteration:
                    line_i = stored_line[i]
                    stored_line[i] = None
                    count_stops += 1
            else:
                line_i = stored_line[i]
                stored_line[i] = None
            # check if the POS for this file is the same as the current min POS
            if line_i[1] == curr_min[1] and line_i[0] == curr_min[0]:
                line[i] = line_i
            else:
                # houston, we might have a problem
                # do we have a new min POS?
                if line_i[1] < curr_min[1]:
                    if curr_min is float('Inf'):
                        # no big deal
                        line[i] = line_i
                    else:
                        # we've found a new min
                        # which is a problem
                        stored_line[i] = line_i
                        # move everything else to storage
                        
                    curr_min = line_i
                elif line_i[1] > curr_min[1]:
                    # houston, we have a problem
                    stored_line[i] = line_i
                # TODO: handle unequal chromosomes
        # create the merged line
        new_line = "" + 
        # check if we've reached the end of all files
        if count_stops == len(files):
            return
        yield line


header = next_line()


# output.close()
# for file in files:
#     file.close()

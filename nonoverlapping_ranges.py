#!/usr/bin/env python

# get the non-overlapping range of genomic regions. Designed especially for GIs predicted by IslandViewer
# print the results in individual files (1 per job id)
# !! WARNING !! DOES NOT WORK IF START OF REGION IS AFTER END OF REGION. NEEDS TO BE IMPLEMENTED
#
# Author: Claire Bertelli (claire.bertelli[@]gmail.com)
# Date: 08.2015
# ---------------------------------------------------------------------------

import argparse
import sys
import pandas
from itertools import chain
import numpy
import re


def join_ranges(data, distance):

    left, right = 1, -1

    sys.stdout.write("Checking if dataframe...")
    if isinstance(data, pandas.DataFrame):
        data = numpy.array(data)
        sys.stdout.write("DataFrame, changing to numpy array\n")

    sys.stdout.write("Sorting data...")
    data = sorted(chain.from_iterable(((start, left), (stop + distance, right)) for start, stop in data))

    sys.stdout.write("Data sorted. \nNow computing non overlapping ranges...\n")
    c = 0
    i = 0
    for value, label in data:
        if c == 0:
            x = value
        c += label
        if c == 0:
            yield "GI_%s" % i, x, value - distance
            i = i+1


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-f", '--file', type=str, nargs='+', help="File(s) containing the data frame of start and stop. By defaults, does not contain a header.")
    parser.add_argument("-d", '--distance', type=int, default=0, help="Authorized distance between two ranges to be considered as a single region. Default is 0, the regions must be overlapping." )
    parser.add_argument("-o", '--output', type=str, help="Name of the output file")
    args = parser.parse_args()
    for file in args.file:
        df = pandas.read_table(file, header=None, usecols=[1,2])

        nonovl_ranges = pandas.DataFrame(list(join_ranges(df, args.distance)))
        #nonovl_ranges.columns = ['start','end']
        
        print "Writing file in %s" % args.output
        nonovl_ranges.to_csv(args.output, sep='\t', index=False, header=None)

    # data to test function
    # print list(join_ranges([[138, 821], [900, 1150], [905, 935]]))
    # print list(join_ranges([[138, 821], [900, 910], [905, 915]], 80))

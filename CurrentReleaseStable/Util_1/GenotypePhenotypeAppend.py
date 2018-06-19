import csv
from itertools import izip

__author__ = "Michael J. Suggs // mjs3607@uncw.edu"


def transpose(marker_file):
    """Transposes an input marker file (in CSV format)

    :param marker_file: The marker file to be transposed
    :return:
    """
    a = izip(*csv.reader(open(marker_file, "rbU")))
    csv.writer(open("transposed_markers.csv".format(), "wb")).writerows(a)


def make_geno_dict(marker_file):
    """Returns all markers in a dictionary format. Note the input file must be transposed.

    :param marker_file: The transposed marker file
    :return: geno_dict, marker_names
    """
    geno_dict = {}

    with open(marker_file, "r") as mf:
        header = mf.readline()
        marker_names = header.split(",")[1:]
        for i in range(4):
            mf.next()

        for line in mf:
            line = line.split(",")
            geno_dict[line[0].lower()] = line[1:]

    return geno_dict, marker_names


def combine_files(geno_dict, marker_names, datafile):
    """Combines the marker and data files.

    :param geno_dict: Dictionary generated from make_geno_dict()
    :param marker_names: Marker names taken from the transposed file via the make_geno_dict() files
    :param datafile: The data file to be appended to
    :return:
    """
    with open(datafile, "rU") as df1, open("{}_appended.csv".format(df1.name), "w") as adf1:
        header = df1.readline().strip().split(',')
        line_column = header.index('Line')
        header += marker_names
        adf1.write(",".join(header))

        for line in df1:
            line = line.strip().split(',')
            if line[line_column].lower() in geno_dict.keys():
                new_line = line + geno_dict[line[line_column].lower().strip()]
            else:
                new_line = line + ['\n']
            adf1.write(",".join(new_line))



if __name__ == '__main__':
    transpose("IBM94markerset08seq.csv")
    geno_dict, marker_names = make_geno_dict("transposed_markers.csv")
    combine_files(geno_dict, marker_names, "Pyear2006fielddata.csv")
    combine_files(geno_dict, marker_names, "RandIBM642006PyearData.csv")

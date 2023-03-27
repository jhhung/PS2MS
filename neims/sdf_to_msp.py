import os
import sys

infile = open(sys.argv[1], 'r')
outfile = open(sys.argv[1][:-3] + "msp", 'w')
is_predict_spectrum = True if sys.argv[2].lower() == "true" else False
record = False
line = infile.readline()
name = ""
mw = ""
spectrum = []
while line:
    #print(line)
    if record:
        if len(line) != 1:
            line = line.strip().split(" ")
            spectrum.append((line[0], line[1]))
        else:
            record = False
    elif "$$$$" in line:
        outfile.write("Name: " + name)
        outfile.write("ID: " + name)
        outfile.write("MW: " + mw)
        outfile.write("Num Peaks: {0}\n".format(len(spectrum)))
        for s in spectrum:
            outfile.write("{0} {1}\n".format(s[0], s[1]))
        outfile.write("\n")
        spectrum = []
    else:
        if "<SMILES>" in line:
            line = infile.readline()
            name = line
        elif "<EXACT MASS>" in line:
            line = infile.readline()
            mw = line
        elif "<MASS SPECTRAL PEAKS>" in line and not is_predict_spectrum:
            record = True
        elif "<PREDICTED SPECTRUM>" in line and is_predict_spectrum:
            record = True
    line = infile.readline()

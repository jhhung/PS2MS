output_id_list = ["Name", "ID", "FP", "MW", "Num Peaks"]

def read_a_spectrum(infile):
    result = {}
    line = infile.readline().strip()
    while line:
        key, value = [item.strip() for item in line.split(':')]
        key = key.upper()
        result[key] = value
        if key == "NUM PEAKS":
            result["Spectrum"] = []
            for i in range(int(value)):
                result["Spectrum"].append([float(item) for item in infile.readline().strip().split(' ')])
        line = infile.readline().strip()
    return result

def write_a_spectrum(outfile, spectrum):
    for id in output_id_list:
        if id.upper() not in spectrum:
            continue
        outfile.write("{0}: {1}\n".format(id, spectrum[id.upper()]))
    for mass, intensity in spectrum["Spectrum"]:
        outfile.write("{0} {1}\n".format(mass, intensity))
    outfile.write("\n")

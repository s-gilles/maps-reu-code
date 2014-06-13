#!/usr/bin/python

# script configuration variables


# Load a CSV file organized by manifold and reorganize it by polynomial and volume.
# The result: dict poly ---> (dict vols ---> ((manifold name, tetrahedra), root), degree etc.)

def read_raw_csv(in_file):
    data = dict()
    for l in in_file.readlines():
            w = l.replace('"','').strip().split(',')
            data.setdefault(w[3],[dict(),w[4]])[0].setdefault(w[2],[list(),w[5]])[0].append(w[0:2])
            data[w[3]].extend(w[6:8])
    return data

# Writes a CSV file containing the mainfolds records as shown below.
def write_csv(out_file, data, append=False):
    if not append:
        out_file.write('InvariantTraceField,Volume,InvariantTraceFieldDegree,NumberOfComplexPlaces,Root,Disc,Name,Tetrahedra\n')
    for poly, polyrec in data.items():
        for vol, volrec in polyrec[0].items():
            for manrec in volrec[0]:
                out_file.write('"'+poly+'",')           # trace field
                out_file.write('"'+vol+'",')            # volume
                out_file.write('"'+polyrec[1]+'",')     # degree
                out_file.write('"'+polyrec[2]+'",')     # complex places
                out_file.write('"'+volrec[1]+'",')      # root
                out_file.write('"'+polyrec[3]+'",')     # discriminant
                out_file.write('"'+manrec[0]+'",')      # name
                out_file.write('"'+manrec[1]+'",')      # number of simplices
                out_file.write('\n')

# Test code
f = open('output.csv','r')
f.readline()
d = read_raw_csv(f)
g = open('newoutput.csv','w')
write_csv(g,d)

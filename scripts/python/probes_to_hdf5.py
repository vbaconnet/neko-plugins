import numpy as np
import h5py

import sys

try:
    fname = sys.argv[1]
except:
    print("Please provide the name of the csv file as an argument")
    exit(1)

# 
# Read the probes header: Nprobes,Nfields,field1,field2,...fieldN
#
with open(fname, "r") as f:
    l = f.readline()
    l = l.split(",") # format at this point is Nprobes,Nfields,f1,f2,f3...
    Nprobes = int(l.pop(0))
    Nfields = int(l.pop(0))
    fields = ",".join(l).split()[0] # remove \n from the names
    fields = fields.split(",")

print(fname, Nprobes, Nfields, fields)

#
# Read the x,y,z coordinates
#
print("Reading coords")
coords = np.genfromtxt(fname, skip_header=1, max_rows = Nprobes, delimiter = ",")
print(coords.shape)
print("Reading data")

names = ['t'] + fields
print(names)

#
# Read the rest of the data
#
data = np.genfromtxt(fname, skip_header = 1+Nprobes, delimiter = ",")
#data = pd.read_csv(fname, names = names, skiprows=Nprobes+1)

#for fname in ['3.csv', '4.csv', '5.csv']:
#    print("reading", fname)
#    temp = pd.read_csv(fname, names = names, skiprows=Nprobes+1)
#    data = pd.concat((data, temp))

print("Output to HDF5")

with h5py.File("probes.hdf5", "w") as f:
    for i in range(Nprobes):
        f.create_group(str(i))
        f["/"].attrs["Nprobes"] = Nprobes
        f["/"].attrs["Nfields"] = Nfields
        f[str(i)].attrs["x"] = coords[i,0]
        f[str(i)].attrs["y"] = coords[i,1] 
        f[str(i)].attrs["z"] = coords[i,2]
        for idx, key in enumerate(names):
            f[str(i)].create_dataset(key, data = data[i::Nprobes,idx])
print("Done")

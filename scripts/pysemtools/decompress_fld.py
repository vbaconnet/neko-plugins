# Import required modules
from mpi4py import MPI #equivalent to the use of MPI_init() in C
import numpy as np

from os.path import join
import sys

# Get mpi info
comm = MPI.COMM_WORLD

# Data types
from pysemtools.datatypes.msh import Mesh
from pysemtools.datatypes.coef import Coef
from pysemtools.datatypes.field import FieldRegistry

from pysemtools.datatypes.utils import create_hexadata_from_msh_fld
# Readers
from pysemtools.io.ppymech.neksuite import pynekread, pwritenek

# Adios2 wrappers
from pysemtools.io.adios2.compress import write_field, read_field

fname = join("compress",'compressed_field0.f00000')

# Instance the empty objects
msh = Mesh(comm, create_connectivity=False)

# Read the data
msh,fld = read_field(comm, fname)

out = create_hexadata_from_msh_fld(msh=msh, fld=fld)
pwritenek("yayaya.f00000", out, comm)

exit(0)

for i in range(int(sys.argv[1])):

    fname = "field0.f{:05d}".format(i)
    fname = join("snaps", fname)
    
    fld = FieldRegistry(comm)
    
    # Read the data
    pynekread(fname, comm, data_dtype=np.single, fld = fld)

    # Write the data in a subdomain and with a different order than what was read
    fout = './compressed_field0.f{:05d}'.format(i)
    fout = join("compress", fout)
    
    wrd_size = 4
    if i == 0:
        write_field(comm, msh=msh, fld=fld, fname=fout, wrd_size=wrd_size, write_mesh = True)
    else:
        write_field(comm, msh=msh, fld=fld, fname=fout, wrd_size=wrd_size, write_mesh = False)

    # Uncomment below to test if the read compressed file is correct
    #comm.Barrier()
    #msh2, fld2 = read_field(comm, fname = fout)
    #print(np.allclose(msh.x, msh2.x))
    #print(np.allclose(fld.fields['vel'][0], fld2.fields['vel'][0]))
    
if comm.Get_rank() == 0: print("DONE")

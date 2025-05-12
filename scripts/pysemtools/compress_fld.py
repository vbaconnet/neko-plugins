# Import required modules
from mpi4py import MPI #equivalent to the use of MPI_init() in C
import numpy as np

# Get mpi info
comm = MPI.COMM_WORLD

# Data types
from pysemtools.datatypes.msh import Mesh
from pysemtools.datatypes.coef import Coef
from pysemtools.datatypes.field import FieldRegistry

# Readers
from pysemtools.io.ppymech.neksuite import pynekread

fname = 'field0.f00000'

# Adios2 wrappers
from pysemtools.io.adios2.compress import write_field, read_field

# Instance the empty objects
msh = Mesh(comm, create_connectivity=False)
fld = FieldRegistry(comm)

# Read the data
pynekread(fname, comm, data_dtype=np.single, msh=msh, fld = fld)

# Write the data in a subdomain and with a different order than what was read
fout = 'compressed_field0.f00001'
wrd_size = 4

print("Before write_field")
write_field(comm, msh=msh, fld=fld, fname=fout, wrd_size=wrd_size, write_mesh=True)
print("After write_field")

comm.Barrier() # This is not needed, in general. Here just ensure we don't overlap the read and write operations

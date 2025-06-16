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
from pysemtools.datatypes.field import Field

from pysemtools.datatypes.utils import create_hexadata_from_msh_fld
# Readers
from pysemtools.io.ppymech.neksuite import pynekread, pwritenek, pynekwrite

# Adios2 wrappers
from pysemtools.io.adios2.compress import write_field, read_field


# ===== Parameters

compressed_snapshots_dir = "./"
decompressed_snapshots_dir = "./"

compressed_prefix = "compressed_field"
decompressed_prefix = "decompressed_field"

start_counter = 0
Nsnaps = 5

# =====

fld = Field(comm)
msh = Mesh(comm)

#
# Read compressed field
#

for i in range(start_counter, start_counter + Nsnaps):

    fname = join(compressed_snapshots_dir, f'{compressed_prefix}0.f'+'{:05d}'.format(i))

    # Read the data
    if i == start_counter:
        #_, fld = read_field(comm, fname)
        read_field(comm, fname, msh=msh, fld=fld, overwrite_fld = True)
    else:
        read_field(comm, fname, fld=fld, overwrite_fld = True)
        #msh, fld = read_field(comm, fname)

    for k in fld.fields.keys():
        print(k, len(fld.fields[k]))

    #
    # Write decompressed field
    #
    fout = join(decompressed_snapshots_dir, f"{decompressed_prefix}0.f" + "{:05d}".format(i))
    if comm.Get_rank() == 0: print("Writing to ", fout, "...")
    if i == start_counter:
        pynekwrite(fout, comm, msh=msh, fld=fld, wdsz=4, write_mesh=True)
    else:
        pynekwrite(fout, comm, msh=msh, fld=fld, wdsz=4, write_mesh=False)
    if comm.Get_rank() == 0: print("Writing to ", fout, "... done")

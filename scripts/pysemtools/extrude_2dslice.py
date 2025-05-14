# Import required modules
from mpi4py import MPI #equivalent to the use of MPI_init() in C
import matplotlib.pyplot as plt
import numpy as np

import sys

try:
    fname = sys.argv[1]
except:
    raise ValueError("Need to put the file name in argument! e.g. field0.f00100")


# Get mpi info
comm = MPI.COMM_WORLD

from pysemtools.io.ppymech.neksuite import pynekread, pynekwrite
from pysemtools.datatypes.msh import Mesh
from pysemtools.datatypes.coef import Coef
from pysemtools.datatypes.field import Field

msh = Mesh(comm, create_connectivity=True)
fld = Field(comm)
pynekread(fname, comm, data_dtype=np.single, msh=msh, fld=fld)
coef = Coef(msh, comm, get_area=False)

# Import helper function
from pysemtools.datatypes.utils import extrude_2d_sem_mesh

# Extrude the 2D mesh to 3D
msh3d, fld3d = extrude_2d_sem_mesh(comm, lz = msh.lx, msh = msh, fld = fld)

fld3d.fields["vel"].append( np.zeros_like( fld3d.fields["vel"][0] ) )

pynekwrite("2d_extruded_"+fname, comm, msh=msh3d, fld=fld3d)

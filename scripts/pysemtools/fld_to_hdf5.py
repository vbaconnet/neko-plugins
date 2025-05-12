from mpi4py import MPI
comm = MPI.COMM_WORLD

import h5py
from pynektools.io.wrappers import read_data
# Data types
from pynektools.datatypes.msh import Mesh
from pynektools.datatypes.coef import Coef
from pynektools.datatypes.field import FieldRegistry
from pynektools.io.ppymech.neksuite import pynekread

from os.path import join
import numpy as np

def read_mesh(fname, comm):
    """
    Reads a SEM mesh (all GLL points) from a field file.

    Needs a comm object from MPI.

    Returns arrays x,y,z and the mass matrix, which have shapes
    (nelv, lx, ly, lz)
    """
    
    # Initialize mesh object
    msh = Mesh(comm, create_connectivity = False)

    # Load the mesh from the field file
    pynekread(fname, comm, msh=msh)

    # Create the coef object to get the mass matrix
    coef = Coef(msh, comm)

    return msh.x, msh.y, msh.z, coef.B

def read_snapshot_velocities(fname, comm):
    """
    Reads (u,v,w) velocities from a snapshot.

    Returns an array containing flattened (u,v,w). 
    """

    data = read_data(comm, fname, ['vel_0', 'vel_1', 'vel_2'])

    # Flatten the arrays
    for key in data.keys():
        data[key] = data[key].flatten()

    return np.concatenate(( data["vel_0"], data["vel_1"], data["vel_2"] )).reshape((len(data["vel_0"]), 3))

in_dir = "../POD/interpolated_snapshots"
out_dir = "./"
fname = "resampled_field10240.f00000"
x,y,z,B = read_mesh(join(in_dir,fname), comm)

size = comm.Get_size()
rank = comm.Get_rank()
print(f"rank {rank} has size mesh {x.shape}")

write = True

if write:

    nel_lcl = x.shape[0]
    nel_glb = comm.allreduce(nel_lcl, op=MPI.SUM) 
    lx = x.shape[1]

    # Write the mesh
    print(f"{rank} Opening hdf5 file")
    if rank == 0: print(f"Global shape: {nel_glb},{lx}")

    fname = f"resampled_SEM_MESH.hdf5"
    f = h5py.File(join(out_dir, fname), "w", driver="mpio", comm=comm)    
    
    if rank == 0: print("Writing mass matrix...")
    dset = f.create_dataset("B", (nel_glb, lx, lx, lx))

    offset = rank * nel_lcl
    dset[ offset : offset + nel_lcl ] = B[:,:,:,:]

    f['B'].attrs['Details'] = "B is the mass matrix, computed using pyNekTools"
    if rank == 0:print("Writing mass matrix... done")
    
    if rank == 0: print("Writing mesh X...")
    dset = f.create_dataset("mesh/x", (nel_glb, lx, lx, lx))

    offset = rank * nel_lcl
    dset[ offset : offset + nel_lcl ] = x[:,:,:,:]
    
    if rank == 0: print("Writing mesh Y...")
    dset = f.create_dataset("mesh/y", (nel_glb, lx, lx, lx))

    offset = rank * nel_lcl
    dset[ offset : offset + nel_lcl ] = y[:,:,:,:]
    
    if rank == 0: print("Writing mesh Z...")
    dset = f.create_dataset("mesh/z", (nel_glb, lx, lx, lx))

    offset = rank * nel_lcl
    dset[ offset : offset + nel_lcl ] = z[:,:,:,:]
    
    f["mesh"].attrs["Details"] = "The group 'mesh' contains the (x,y,z) coordinates of all GLL points. Shape: nelv,lx,lx,lx"
    f["mesh"].attrs["nelv"] = nel_glb
    f["mesh"].attrs["lx"] = lx
    if rank == 0: print("Writing mesh... done")
    
    f.close()
    print(f"{rank} Done!")

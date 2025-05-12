"""

Use this program as python3 pod_pysem.py inputs.json

"""

# Adios2 wrappers
from pysemtools.io.adios2.compress import read_field


# Import general modules
import sys
import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"

from os.path import join

import json
import numpy as np

# -- Matplotlib config
import matplotlib.pyplot as plt

plt.rcParams['figure.figsize'] = [10, 6]
plt.rcParams['figure.dpi'] = 300
plt.rcParams['font.family'] = 'serif'
plt.rcParams['text.usetex'] = False
plt.rcParams['font.size'] = 20
default_cycle = plt.rcParams['axes.prop_cycle']
# --


from scipy.signal import welch

# Import MPI
from mpi4py import MPI #equivalent to the use of MPI_init() in C

# Split communicator for MPI - MPMD
worldcomm = MPI.COMM_WORLD
worldrank = worldcomm.Get_rank()
worldsize = worldcomm.Get_size()
col = 1
comm = worldcomm.Split(col,worldrank)
rank = comm.Get_rank()
size = comm.Get_size()

# Open input file to see path
f = open(sys.argv[1], "r")
params_file = json.loads(f.read())
f.close()

OUTPUT_DIR = params_file["IO"]["output_dir"]

# Read the POD inputs
pod_number_of_snapshots = int(params_file["number_of_snapshots"])
pod_batch_size = int(params_file["batch_size"])
pod_keep_modes = int(params_file["keep_modes"])
pod_write_modes = int(params_file["write_modes"])

pod_fields = params_file["fields"]
number_of_pod_fields = len(pod_fields)

dtype_string = params_file.get("dtype", "double")
backend = params_file.get("backend", "numpy")
if dtype_string == "single":
    dtype = np.float32
else:
    dtype = np.float64

# Import IO helper functions
from pysemtools.io.utils import get_fld_from_ndarray, IoPathData
# Import modules for reading and writing
from pysemtools.io.ppymech.neksuite import pynekread, pynekwrite, pwritenek
# Import the data types
from pysemtools.datatypes.msh import Mesh
from pysemtools.datatypes.coef import Coef
from pysemtools.datatypes.field import Field, FieldRegistry
from pysemtools.datatypes.utils import create_hexadata_from_msh_fld
# Import types asociated with POD
from pysemtools.rom.pod import POD
from pysemtools.rom.io_help import IoHelp

# Start time
start_time = MPI.Wtime()

# Read the data paths from the input file
mesh_data = IoPathData(params_file["IO"]["mesh_data"])
field_data = IoPathData(params_file["IO"]["field_data"])

# Instance the POD object
pod = POD(comm, number_of_modes_to_update = pod_keep_modes, global_updates = True, auto_expand = False, bckend = backend)

# Initialize the mesh file
path     = mesh_data.dataPath
casename = mesh_data.casename
index    = mesh_data.index
fname    = path+casename+'0.f'+str(index).zfill(5)
#msh = Mesh(comm)
#pynekread(fname, comm, data_dtype=dtype, msh = msh)
msh, _ = read_field(comm, fname)

# Initialize coef to get the mass matrix
coef = Coef(msh, comm)
bm = coef.B

# Instance io helper that will serve as buffer for the snapshots
ioh = IoHelp(comm, number_of_fields = number_of_pod_fields, batch_size = pod_batch_size, field_size = bm.size, field_data_type=dtype, mass_matrix_data_type=dtype)

# Put the mass matrix in the appropiate format (long 1d array)
mass_list = []
for i in range(0, number_of_pod_fields):
    mass_list.append(np.copy(np.sqrt(bm)))
ioh.copy_fieldlist_to_xi(mass_list)
ioh.bm1sqrt[:,:] = np.copy(ioh.xi[:,:])

#mean = np.zeros_like(ioh.xi)

# -- Define figures for plotting

plt.figure(1) # For singular values
plt.xlabel("Mode number")
plt.ylabel(r"$\sigma_i / \Sigma \sigma_i$")

plt.figure(2) # For reconstruction of percentages of sing. values
plt.xlabel("Number of snapshots")
plt.ylabel("Mode number")

plt.figure(3) # For reconstruction of percentages of sing. values
plt.xlabel("Frequency [-]")
plt.ylabel("Power density")

#

j = 0
rank_total_energy = 0.0
while j < pod_number_of_snapshots:

    # Recieve the data from fortran
    path     = field_data.dataPath
    casename = field_data.casename
    index    = field_data.index
    fname=path+casename+'0.f'+str(index + j).zfill(5)
    #fld = FieldRegistry(comm)
    #pynekread(fname, comm, data_dtype=dtype, fld = fld) 
    _, fld = read_field(comm, fname)

    # Get the required fields
    #u = fld.registry['u']
    #v = fld.registry['v']
    #w = fld.registry['w']
    u = fld.fields['vel'][0]
    v = fld.fields['vel'][1]
    w = fld.fields['vel'][2]

    # Compute snapshot energy
    snapshot_energy = np.sum(  u**2 * coef.B)
    snapshot_energy += np.sum( v**2 * coef.B)
    snapshot_energy += np.sum( w**2 * coef.B)
    rank_total_energy += snapshot_energy

    # Put the snapshot data into a column array
    ioh.copy_fieldlist_to_xi([u, v, w])

    # Load the column array into the buffer
    ioh.load_buffer(scale_snapshot = True)

    #mean += ioh.xi

    # Update POD modes
    if ioh.update_from_buffer:
        pod.update(comm, buff = ioh.buff[:,:(ioh.buffer_index)])       

        # Write the singular values and vectors for convergence testing
        if comm.Get_rank() == 0:
            suff_idx = j//pod_batch_size
            np.save(join(OUTPUT_DIR, f"singular_values_{suff_idx}"), pod.d_1t)
            print(f"Wrote {OUTPUT_DIR}/singular values_{suff_idx}")
            np.save(join(OUTPUT_DIR, f"right_singular_vectors_{suff_idx}"), pod.vt_1t)
            print(f"Wrote {OUTPUT_DIR}/right signular values_{suff_idx}")

            # Singular values alone
            plt.figure(1)
            sing_val_energy = np.cumsum(pod.d_1t[1:]) / np.sum(pod.d_1t[1:])
            plt.loglog( np.arange(1,len(pod.d_1t),1), sing_val_energy, "-o", 
                     label = f"{j+1} snapshots", markerfacecolor="white" )
            plt.legend()
            plt.savefig(join(OUTPUT_DIR, "singular_values.png"))

            # Find the # of modes necessary to reconstruct p percent of the energy
            plt.figure(2)
            ax = plt.gca()
            ax.set_prop_cycle(default_cycle)
            for c_idx, p in enumerate([0.3, 0.5, 0.75, 0.85, 0.95]):
                n = np.argmax( sing_val_energy > p ) # returns the index of the first element > p
                plt.plot(j+1, n, "o", label = f"{p*100:.1f}% reconstructed")
            if suff_idx == 0: plt.legend()
            plt.savefig(join(OUTPUT_DIR, "singular_values_reconstruct.png"))

            # Plot the power spectrum of all the modes to write
            for m in range(1,pod_write_modes):
                plt.figure(f"power_spectrum_{m}")
                f, Pxx_den = welch(pod.vt_1t[m,:], 1/0.021, nperseg=32)
                plt.semilogy(f, Pxx_den, label = f"{j+1} snapshots")
                plt.legend()
                plt.savefig(join(OUTPUT_DIR, f"power_spectrum_mode_{m}.png"))

        # Also write the total energy captured to check convergence
        total_energy = np.array(0.0, dtype=np.single)
        tmp_rank_total_energy = np.array(rank_total_energy, dtype=np.single)
        comm.Reduce(tmp_rank_total_energy, total_energy, op=MPI.SUM, root=0)
        
        if comm.Get_rank() == 0:
            print(f">>>>> BATCH {suff_idx}")
            print(f'Total energy in snapshots: {total_energy}')
            print(f'Energy captured in POD   : {np.sum(pod.d_1t**2)}')
            print(f'Energy of mode 0 (mean)  : {pod.d_1t[0]**2}')
            print(f'TKE in snapshots         : {total_energy - pod.d_1t[0]**2}')
            print(f'P_tke captured by POD    : {np.sum(pod.d_1t[1:]**2)}')
            print(f' => {np.sum(pod.d_1t[1:]**2) / (total_energy - pod.d_1t[0]**2) * 100} % ')

    j += 1

#mean /= j+1

# Check if there is information in the buffer that should be taken in case the loop exit without flushing
if ioh.buffer_index > ioh.buffer_max_index:
    ioh.log.write("info","All snapshots where properly included in the updates")
else: 
    ioh.log.write("warning","Last loaded snapshot to buffer was: "+repr(ioh.buffer_index-1))
    ioh.log.write("warning","The buffer updates when it is full to position: "+repr(ioh.buffer_max_index))
    ioh.log.write("warning","Data must be updated now to not lose anything,  Performing an update with data in buffer ")
    pod.update(comm, buff = ioh.buff[:,:(ioh.buffer_index)])

# Scale back the modes
pod.scale_modes(comm, bm1sqrt = ioh.bm1sqrt, op = "div")

# Rotate local modes back to global, This only enters in effect if global_update = false
pod.rotate_local_modes_to_global(comm)

# Finalize the total energy contained in our snapshots
total_energy = np.array(0.0, dtype=np.single)

rank_total_energy = np.array(rank_total_energy, dtype=np.single)
comm.Reduce(rank_total_energy, total_energy, op=MPI.SUM, root=0)
if comm.Get_rank() == 0: 
    print(f'Total energy in snapshots: {total_energy}')
    print(f'Energy captured in POD   : {np.sum(pod.d_1t**2)}')
    print(f'Energy of mode 0 (mean)  : {pod.d_1t[0]**2}')
    print(f'TKE in snapshots         : {total_energy - pod.d_1t[0]**2}')
    print(f'P_tke captured by POD    : {np.sum(pod.d_1t[1:]**2)}')
    print(f' => {np.sum(pod.d_1t[1:]**2) / (total_energy - pod.d_1t[0]**2) * 100} % ')

# Write the data out
for j in range(0, pod_write_modes):

    if (j+1) < pod.u_1t.shape[1]:

        ## Split the snapshots into the proper fields
        field_list1d = ioh.split_narray_to_1dfields(pod.u_1t[:,j])
        u_mode = get_fld_from_ndarray(field_list1d[0], msh.lx, msh.ly, msh.lz, msh.nelv)
        v_mode = get_fld_from_ndarray(field_list1d[1], msh.lx, msh.ly, msh.lz, msh.nelv)
        w_mode = get_fld_from_ndarray(field_list1d[2], msh.lx, msh.ly, msh.lz, msh.nelv)

        ## Create an empty field and update its metadata
        out_fld = Field(comm)
        out_fld.fields["scal"].append(u_mode)
        out_fld.fields["scal"].append(v_mode)
        out_fld.fields["scal"].append(w_mode)
        out_fld.update_vars()

        ## Create the hexadata to write out
        out_data = create_hexadata_from_msh_fld(msh = msh, fld = out_fld)

        ## Write out a file
        fname = "modes0.f"+str(j).zfill(5)
        pwritenek(join(OUTPUT_DIR,fname),out_data, comm)
        if comm.Get_rank() == 0: print("Wrote file: " + join(OUTPUT_DIR,fname))

# Write the singular values and vectors
if comm.Get_rank() == 0:
    np.save(join(OUTPUT_DIR,"singular_values"), pod.d_1t)
    print("Wrote singular values")
    np.save(join(OUTPUT_DIR,"right_singular_vectors"), pod.vt_1t)
    print("Wrote right signular values")
    with open("snapshot_energy.txt", "w") as f:
        print(f'{pod_number_of_snapshots} {total_energy}')
# End time
end_time = MPI.Wtime()
# Print the time
if comm.Get_rank() == 0:
    print("Time to complete: ", end_time - start_time)

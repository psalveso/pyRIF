# pyRIF
Using Rotamer Interaction Fields from RIFGen/Dock in python

### Installation
`pip install pyRIF`

### Requirements
[PyRosetta](http://www.pyrosetta.org/home) installed in python environment <br/>[RIFDock](https://github.com/rifdock/rifdock) for generating RIFs

### generating HDF5s of RIFs from RIFGen

include `-dump_rifgen_hdf5` in your `rifgen.flag` file when running RIFGen. this will produce `rif.h5` in your output directory. Then run `python tools/convert_rif_h5.py rif.h5`. This will produce `py_rif.h5`. This is the HDF5 file that will be needed to initialize `RotamerInteractionField()`. You will also need the path to the `rotamer_index_spec.txt` and `target.pdb.gz` created by RIFGen. These will be in the output directory you specified in `rifgen.flag`.

### Pose Example
Moves your input pose into the region of the RIF, by aligning the `target_selector` residues onto the model specified with `L_AA_RIF['target']`
```
import glob
import pyrosetta
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector

from pyRIF import RotamerInteractionField

pyrosetta.init()


# dictionary pointing to RIFGen outputs
L_AA_RIF = {
    'HDF5'   : '/path/to/py_rif.h5',
    'rots'   : '/path/to/rotamer_index_spec.txt',
    'target' : '/path/to/target.pdb.gz',
}


# residue selectors to select target and binder residues
target_selector = ChainSelector('A')
binder_selector = ChainSelector('B')


# create RIF object outside of loop
RIF = RotamerInteractionField(
    L_AA_RIF_kwargs=L_AA_RIF,
    residue_selector=binder_selector,
    target_selector=target_selector,
)

for PDB in glob.iglob('/path/to/pdbs/*.pdb'):
    pose = pyrosetta.io.pose_from_pdb(PDB)

    # apply the RIF
    STATUS, RIF_score, sequence_mapping = RIF.apply(pose)

    if STATUS:
        print(f'pass, {RIF_SCORE}\n{sequence_mapping}')
        # continue with remainder of protocol
    else:
        print('fail')

```

### Numpy Example
Assumes your XYZs are already within the region of the RIF
```
import glob
import numpy as np
import pyrosetta

from pyRIF import RotamerInteractionField
pyrosetta.init()


# dictionary pointing to RIFGen outputs
L_AA_RIF = {
    'HDF5'   : '/path/to/py_rif.h5',
    'rots'   : '/path/to/rotamer_index_spec.txt',
    'target' : '/path/to/target.pdb.gz',
}


# create RIF object outside of loop
RIF = RotamerInteractionField(
    L_AA_RIF_kwargs=L_AA_RIF,
)

binder_xyzs = np.random.rand(1000, 100, 3, 3)# [(1000 proteins), (100 residues), (N CA C), (X Y Z)]

for i in range(binder_xyzs.shape[0]):

    STATUS, RIF_score, sequence_mapping = RIF.search_xyzs(binder_xyzs[i, :, :, :])

    if STATUS:
        print(f'pass, {RIF_SCORE}\n{sequence_mapping}')
        # continue with remainder of protocol
    else:
        print('fail')

```
### Note about packing
Right now, packing is done by first pruning rotamers at a hit where the rotamer clashes in a voxel grid with the N CA C atoms of the binder (semi-equivalent to 1-body). Each pair of pruned rotamers at position i and j are then clashed checked against one another in the same voxel grid (semi-equivalent to 2-body). The lowest RIF_score set of rotamers at each hit position that survives these voxel clash checks is taken as the output sequence, and the sum of their RIF_scores is the output RIF_score for this pose/xyzs. This score is not the same as the RIF_score output by RIFDock.

### Note about homog
If you are getting an error on `import pyRIF`, this is due to a bug in the package `homog`. See this [pull-request](https://github.com/willsheffler/homog/commit/526c3f07c720f76333bc8be0cd64b436015ff509) for the fix. In short, on line 11 of `homog/util.py` change `nopython=1, fastmath=1` to `nopython=True, fastmath=True`.

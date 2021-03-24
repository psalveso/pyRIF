# pyRIF
Using Rotamer Interaction Fields from RIFGen/Dock in python


### Installation
`pip install pyRIF`

### Pose Example
Moves your input pose into the region of the RIF
```
import glob
import pyrosetta
from pyrosetta.rosetta.core.select.residue_selector import ChainSelector

from pyrif import RotamerInteractionField

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

from pyrif import RotamerInteractionField
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

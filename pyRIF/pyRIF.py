import numpy as np
import getpy as gp
import nerf
import h5py
import xbin

# pyrosetta
from pyrosetta.io import pose_from_pdb
from pyrosetta.rosetta.core.select.residue_selector import TrueResidueSelector
from pyrosetta.rosetta.core.select.residue_selector import PhiSelector
from pyrosetta.rosetta.core.conformation import ResidueFactory
from pyrosetta.rosetta.core.chemical import ChemicalManager

# for superimpose
from pyrosetta.rosetta.utility import vector1_numeric_xyzVector_double_t as rosetta_xyz_vec
from pyrosetta.rosetta.numeric import xyzMatrix_double_t
from pyrosetta.rosetta.numeric import xyzVector_double_t
from pyrosetta.rosetta.protocols.toolbox import superposition_transform
from pyrosetta.rosetta.protocols.toolbox import apply_superposition_transform





class RotamerInteractionField(object):
    '''
    methods:
        load_RIF_Table(py_rif.h5, rotamer_index_spec.txt) -> python dictionary of the RIF table
        RIF_stubs(Ns, CAs, Cs,)                           -> the (4,4) homogenous transforms of the corresponding cooridnates
        superimpose(target, moving)                       -> aligns moving onto target
        get_irots(L_offsets, D_offsets)                   -> given L and D offsets, returns the score_list and irot list of hit. Run prepacking
        fast_pack()                                       -> given lots of things, does voxel-based clash packing. returns a sequence map and total rif score
        search_xyzs(xyzs)                                 -> given numpy array of NCAC XYZs, returns True/False, RIF_score, sequence_mapping
        apply(pose)                                       -> given pose, gets XYZs, runs search_xyzs()
    '''
    def __init__(
        self,
        residue_selector=TrueResidueSelector(),
        target_selector=TrueResidueSelector(),
        min_RIF_hits=1,
        min_HOTSPOT_hits=0,
        min_RIF_score=0.0,
        HOTSPOTs_in_RIF=0,
        L_AA_RIF_kwargs={'HDF5' : None, 'rots' : None, 'target' : None},
        D_AA_RIF_kwargs={'HDF5' : None, 'rots' : None, 'target' : None},
        atom_selector=('N', 'CA', 'C'),
        fast_pack_kwargs={'pert_mag' : 0.5, 'shake_repeats' : 10, 'resolution' : 1.0, 'rotamer_clash_tolerance' : 0},
    ):
        self.D_three_to_four = {
            'DAL' : 'DALA',
            'DPH' : 'DPHE',
            'DIL' : 'DILE',
            'DLE' : 'DLEU',
            'DVA' : 'DVAL',
            'DCS' : 'DCYS',
            'DAS' : 'DASP',
            'DGU' : 'DGLU',
            'DHI' : 'DHIS',
            'DLY' : 'DLYS',
            'DME' : 'DMET',
            'DAN' : 'DASN',
            'DPR' : 'DPRO',
            'DGN' : 'DGLN',
            'DAR' : 'DARG',
            'DSE' : 'DSER',
            'DTH' : 'DTHR',
            'DTR' : 'DTRP',
            'DTY' : 'DTYR',
        }
        # load the RIFs
        self.L_AA_RIF = self.load_RIF_Table(**L_AA_RIF_kwargs)
        self.D_AA_RIF = self.load_RIF_Table(**D_AA_RIF_kwargs)
        # combine L and D iROTs in one location
        self.IROT_dict = {}
        self.IROT_dict['irot_DOFS'] = self.L_AA_RIF['irot_DOFS'].copy()
        self.IROT_dict['irot_DOFS'].extend(self.D_AA_RIF['irot_DOFS'].copy())
        self.IROT_dict['rotamer_lines'] = self.L_AA_RIF['rotamer_lines'].copy()
        self.IROT_dict['rotamer_lines'].extend(self.D_AA_RIF['rotamer_lines'].copy())
        assert len(self.IROT_dict['irot_DOFS']) == len(self.IROT_dict['rotamer_lines'])
        # set the parameters for the search
        self.min_RIF_hits = min_RIF_hits
        self.min_RIF_score = min_RIF_score
        self.min_HOTSPOT_hits = min_HOTSPOT_hits

        if HOTSPOTs_in_RIF > 0:
            self.hotspot_irots = np.arange(len(self.L_AA_RIF['rotamer_lines']))[-HOTSPOTs_in_RIF:]
        else:
            self.hotspot_irots = np.array([])


        self.fast_pack_params = fast_pack_kwargs

        # define selectors used to grab the binder vs the target
        self.atom_selector = atom_selector
        self.residue_selector = residue_selector
        self.target_selector = target_selector
        self.L_selector = PhiSelector()
        self.D_selector = PhiSelector()
        self.L_selector.set_select_positive_phi(False)
        self.D_selector.set_select_positive_phi(True)


    def load_RIF_Table(
        self,
        HDF5=None,
        rots=None,
        target=None,
    ):
        '''
        Given the path to the rif.h5, and path the rotamer lines file,
        loads the required data,
        returns a dictionary of this data structures
        '''
        RIF = {}
        rif_dict = gp.Dict(np.dtype('int64'), np.dtype('int64'), 0)

        if HDF5 == None or rots == None or target == None:
            RIF['rif_dict'] = rif_dict
            RIF['scores'] = [0]
            RIF['irots'] = [-1]
            RIF['binner'] = xbin.XformBinner()
            RIF['rotamer_lines'] = []
            RIF['sats'] = None
            RIF['irot_DOFS'] = []
            RIF['target'] = None

            return RIF

        with h5py.File(HDF5, 'r') as f:
            rif_dict[f['xbin_key'][:]] = f['offsets'][:]
            scores = np.zeros(f['scores'].shape, np.float16)
            irots = np.zeros(f['irots'].shape, np.int16)
            scores[:] = f['scores'][:]
            irots[:] = f['irots'][:]
            binner = xbin.XformBinner(*list(f['cart_ori_bound']))
            if 'sats' not in list(f.keys()):
                sats = None
            else:
                sats = np.zeros(f['sats'].shape, np.int16)
                sats[:] = f['sats'][:]

        rotamer_lines = []
        with open(rots) as f:
            for line in f:
                rotamer_lines.append(
                    line.strip().split('\t')
                )

        # initialize the things we need to build residues
        chemical_manager = ChemicalManager.get_instance()
        ResidueTypeSet = chemical_manager.residue_type_set( 'fa_standard' )

        backbone_atoms = [' C  ', ' N  ',  ' CA ', ' O  ', ' H  ', ' HA ']
        sc_dof_list = []
        for rotamer in rotamer_lines:
            # create the residue in pyrosetta
            residue = ResidueFactory.create_residue( ResidueTypeSet.name_map( self.D_three_to_four.get(rotamer[0], rotamer[0]) ) )
            # set the torions of the residue
            for i, torsion in enumerate(rotamer[1:]):
                chi_idx = i + 1
                if chi_idx > residue.nchi(): break
                torsion = float(torsion)
                residue.set_chi(chi_idx, torsion)

            # get the XYZs of the C N CA
            frame_xyzs = []
            for atom_name in backbone_atoms[:3]:
                frame_xyzs.append(
                    residue.xyz(atom_name),
                )

            # grab the coordinates of all side chain atoms
            side_chain_XYZs = []
            side_chain_atoms = []
            for i, atom in enumerate(residue.atoms()):
                atom_idx = i + 1
                if residue.atom_name(atom_idx) in backbone_atoms: continue
                side_chain_XYZs.append(residue.xyz(atom_idx))
                side_chain_atoms.append(residue.atom_name(atom_idx))

            # NeRF to get the DOFs
            frame_xyz = np.array(frame_xyzs)[np.newaxis, :, :]
            sc_xyz = np.array(side_chain_XYZs)[np.newaxis, :, :]
            res_xyz = np.concatenate((frame_xyz, sc_xyz), axis=1)
            sc_dof = nerf.iNeRF(res_xyz)
            sc_dof_list.append(sc_dof)

        assert len(rotamer_lines) == len(sc_dof_list)

        RIF = {}
        RIF['rif_dict'] = rif_dict
        RIF['scores'] = scores
        RIF['irots'] = irots
        RIF['binner'] = binner
        RIF['rotamer_lines'] = rotamer_lines
        RIF['sats'] = sats
        RIF['irot_DOFS'] = sc_dof_list
        RIF['target'] = pose_from_pdb(target)

        return RIF

    def superimpose(
        self,
        target,# pose object
        moving,# a second pose object
    ):
        """
        Given a pose of the target and a second pose object,

        CA alignment of moving onto target
        """
        # apply target_residue selector to moving
        target_res_idx = self.target_selector.apply(moving)

        #store the lengths of the two poses
        target_length = len(target.residues)
        moving_length = len(moving.residues)

        #### assuming that res 1, N in target are in moving
        assert target_length <= moving_length, 'RIFgen target model larger than pose, superimpose will fail'

        # make the rosetta object to hold xyz coordinates of the residue we want to align
        target_match_coords = rosetta_xyz_vec()
        moving_coords = rosetta_xyz_vec()

        # add the coordinates of the residues to be aligned to their respctive holders
        #iterate over all residues in target
        # we do this seperately incase user has weird chain ordering for target protein
        # in the moving pose
        for i in range(moving_length):
            position = i + 1
            if target_res_idx[position]:
                moving_coords.append(
                    moving.residues[position].xyz('CA')
                )
        for i in range(target_length):
            position = i + 1
            target_match_coords.append(
                target.residues[position].xyz('CA')
            )

        assert len(target_match_coords) == len(moving_coords), 'mismatch in number of CAs in between RIFgen target model and pose. Check definition of target_selector.'

        # make the rotation matrix and pose center rosetta objects
        rotation_matrix = xyzMatrix_double_t()
        target_center = xyzVector_double_t()
        moving_center = xyzVector_double_t()

        superposition_transform(
            moving_coords,
            target_match_coords,
            rotation_matrix,
            moving_center,
            target_center,
        )

        apply_superposition_transform(
            moving,
            rotation_matrix,
            moving_center,
            target_center,
        )

        return moving

    def RIF_stubs(self, Ns, CAs, Cs):
        """
        From bcov, how the stubs are calcualted in RIFGen/Dock
        """
        #by_res = npose.reshape(-1, R, 4)
        #Ns = by_res[:,N,:3]
        #CAs = by_res[:,CA,:3]
        #Cs = by_res[:,C,:3]
        assert Ns.shape == CAs.shape == Cs.shape
        e1 = (Cs+Ns)/2 - CAs
        e1 = e1 / np.linalg.norm(e1, axis=1)[...,None]
        e3 = np.cross( e1, Cs - CAs )
        e3 = e3 / np.linalg.norm(e3, axis=1)[...,None]
        e2 = np.cross( e3, e1 )
        e2 = e2 / np.linalg.norm(e2, axis=1)[...,None]
        frames = np.zeros((len(Cs), 4, 4))
        frames[:,:3,0] = e1
        frames[:,:3,1] = e2
        frames[:,:3,2] = e3
        frames[:,3,3] = 1.0
        t = frames[:,:3,:3] @ np.array([-1.952799123558066, -0.2200069625712990, 1.524857]) + CAs
        frames[:,:3,3] = t
        return frames

    def get_irots(
        self,
        L_offsets,
        D_offsets,

    ):
        '''
        give the L and D offsets,
        returns the score list, and irot index lest, matched to the bins
        '''
        # track hotspots
        hotspot_count = 0
        # make list of lists for holding rotamers and scores
        rotamer_list = [[] for i in range(len(L_offsets) + len(D_offsets))]
        score_list = [[] for i in range(len(L_offsets) + len(D_offsets))]

        # define what a failed hotspot match looks like
        failed_hotspot = np.array([-1, -1], dtype=np.int16)
        # create varaible to store nmber of hotspots that hit
        num_hotspots = 0

        # add the index to the dict, with an empty array
        num_L_irots = len(self.L_AA_RIF['irot_DOFS'])

        # this is to help us get pas situations where there is no L_aa hits...
        i = -1
        j = -1

        # start loading the irots and scores into rotamer_list and score_list
        # L AA, This could be in its own method. not doing that now, as its a few things I need to pass back and forth
        for i in range(len(L_offsets)):
            offset = L_offsets[i]
            while ( self.L_AA_RIF['irots'][offset] != -1):
                irot = self.L_AA_RIF['irots'][offset]
                score = self.L_AA_RIF['scores'][offset]
                rotamer_list[i].append(irot)
                score_list[i].append(score)
                # hot spot search
                if self.min_HOTSPOT_hits > 0:
                    if irot in self.hotspot_irots and self.L_AA_RIF['sats'] is not None:
                        if not np.array_equal(self.L_AA_RIF['sats'][offset, :], failed_hotspot):
                            #### this irot is a HOTSPOT
                            hotspot_count += 1
                            # clear all irots at this position except this hotspot
                            rotamer_list[i] = [irot]
                            score_list[i] = [score]
                            # break out of the while loop at this positon
                            break
                offset += 1

        # now look up and store hits in D_RIF
        # using j in the for loop, so that I can add the D hits to the rotamer-list and score_list
        for j in range(len(D_offsets)):
            offset = D_offsets[j]
            iter = 0
            while ( self.D_AA_RIF['irots'][offset] != -1):
                irot = self.D_AA_RIF['irots'][offset] + num_L_irots# this should push the irot index up to correct positon when i stack the two arrays
                score = self.D_AA_RIF['scores'][offset]
                rotamer_list[i+j+1].append(irot)
                score_list[i+j+1].append(score)
                offset += 1

        return rotamer_list, score_list, hotspot_count

    def fast_pack(
        self,
        list_of_rotamers=None,
        scores=None,
        hit_frames=None,
        macrocycle_xyz=None,
        rotamer_file=None,
        rotamer_dofs=None,
        macrocycle_bins=None,
        pert_mag=1.0,
        shake_repeats=10,
        resl=1.0,
        rotamer_clash_tolerance=0
    ):
        '''
        Does the packing of all RIF rotamers that hit...
        returns the total rif score, and the sequence map...

        This is VERY hacky...
        '''

        ###### define 1 body and 2 body packing
        def _build_RIF_hit(
            hit_frame,
            rotamer_dofs,
            macrocycle_bins,
            pert_mag=1.0,
            shake_repeats=10,
            resl=1.0,
            irot=0,
        ):
            """
            given macroycle XYZs,
            the frame of the hit,
            the rotamer torsions of the hits

            builds the CB -> Cn atoms using nerf, clash checks those against the macrocycle XYZs
            """
            # make everything the same axis
            hit_frame = hit_frame[np.newaxis, :, :]

            # redefine the hit frame
            C_N_CA_frame = hit_frame.copy()
            C_N_CA_frame[:,0,:] = hit_frame[:,2,:]
            C_N_CA_frame[:,1,:] = hit_frame[:,0,:]
            C_N_CA_frame[:,2,:] = hit_frame[:,1,:]

            SC_xyz = nerf.NeRF(rotamer_dofs, abcs=C_N_CA_frame)

            # clash check SC_xyzs against macrocycle_XYZ
            # we only send [:,3:,:] as first 3 atoms are backbone.
            SC_bins = bin_coordinates(
                SC_xyz[:,3:, :],
                pert_mag=pert_mag,
                shake_repeats=1000,
                resl=resl,
            )

            clash_count = check_bins(macrocycle_bins, SC_bins)
            return clash_count

        def _two_body_clash(
            hit_frame_i,
            hit_frame_j,
            rotamer_dofs_i,
            rotamer_dofs_j,
            pert_mag=1.0,
            shake_repeats=10,
            resl=1.0,
        ):
            '''
            build inversre rotamers i and j on their coresponding frame,
            returns voxel clash between them
            '''
            # make everything the same shape
            hit_frame_i = hit_frame_i[np.newaxis, :, :]
            hit_frame_j = hit_frame_j[np.newaxis, :, :]

            # redefine the hit frame
            C_N_CA_frame_i = hit_frame_i.copy()
            C_N_CA_frame_i[:,0,:] = hit_frame_i[:,2,:]
            C_N_CA_frame_i[:,1,:] = hit_frame_i[:,0,:]
            C_N_CA_frame_i[:,2,:] = hit_frame_i[:,1,:]

            C_N_CA_frame_j = hit_frame_j.copy()
            C_N_CA_frame_j[:,0,:] = hit_frame_j[:,2,:]
            C_N_CA_frame_j[:,1,:] = hit_frame_j[:,0,:]
            C_N_CA_frame_j[:,2,:] = hit_frame_j[:,1,:]

            i_xyz = nerf.NeRF(rotamer_dofs_i, abcs=C_N_CA_frame_i)
            j_xyz = nerf.NeRF(rotamer_dofs_j, abcs=C_N_CA_frame_j)

            # we only send [:,3:,:] as first 3 atoms are backbone.
            bins_i = bin_coordinates(
                i_xyz[:,3:,:],
                pert_mag=pert_mag,
                shake_repeats=1000,
                resl=resl,
            )
            bins_j = bin_coordinates(
                j_xyz[:,3:,:],
                pert_mag=pert_mag,
                shake_repeats=1000,
                resl=resl,
            )

            return check_bins(bins_i, bins_j)

        # first, make boolean np arrays that hold wether we can consider this rotamer at this positon
        packable = [np.repeat(1, len(x)).astype(np.bool) for x in list_of_rotamers]

        ## tests
        #for list in list_of_rotamers: print(f'LIST pre check: {list}')
        #for arr in packable: print(f'BOOL pre check: {arr}')


        # first mark rotamers false if they clash with backbone
        for i in range(len(list_of_rotamers)):

            for j, irot in enumerate(list_of_rotamers[i]):

                ### skip positons that are empty!!!!
                if list_of_rotamers[i] == []: continue

                rotamer_clash = _build_RIF_hit(
                    hit_frames[i, :, :],# the N CA C frame of the hit we are on
                    rotamer_dofs[irot],# the DOF for this hit
                    macrocycle_bins,# binned cooreds of backbone atoms in macrocycle
                    pert_mag=pert_mag,
                    shake_repeats=10,
                    resl=resl,
                    irot=irot,
                )

                if rotamer_clash > rotamer_clash_tolerance:
                    packable[i][j] = False

        ## test post 1-body, need to enure that there are still enough interactions possible!!!!!
        any_packable = 0
        for arr in packable: any_packable += np.sum(arr)

        if any_packable == 0: return 1000, None

        ## now test 2-body clashes
        for pos_i in range(len(list_of_rotamers)):
            for pos_j in range(len(list_of_rotamers)):

                # do not evaulate interaction between rotamers at same postion
                if pos_i == pos_j: continue
                # skip empty positions!!!!
                if list_of_rotamers[pos_i] == [] or list_of_rotamers[pos_j] == []: continue

                for i, rot_i in enumerate(list_of_rotamers[pos_i]):

                    # skip if pos_i 1-body clash
                    if not packable[pos_i][i]: continue

                    for j, rot_j in enumerate(list_of_rotamers[pos_j]):

                        # skip if pos_j 1-body clash
                        if not packable[pos_j][j]: continue
                        #print(f'{pos_i}: {i} {list_of_rotamers[pos_i]}')
                        #print(f'{pos_j}: {j} {list_of_rotamers[pos_j]}')
                        # evaluate 2-body between pos-i rot_i and pos_j rot_j
                        ij_clash_count = _two_body_clash(
                            hit_frames[pos_i, :, :],
                            hit_frames[pos_j, :, :],
                            rotamer_dofs[rot_i],
                            rotamer_dofs[rot_j],
                            pert_mag=pert_mag,
                            shake_repeats=10,
                            resl=resl,
                        )
                        #print(f'IJ CLASH: {ij_clash_count}')
                        if ij_clash_count > rotamer_clash_tolerance:
                            packable[pos_i][i] = False
                            packable[pos_j][j] = False

        ## test if there are enough 2 bodys to give us the RIF thrshold...
        any_packable = 0
        TOTAL_RIF_SCORE = 0
        SEQ_MAP = []
        for i, arr in enumerate(packable):
            if arr.size == 0: continue
            any_packable += np.sum(arr)
            #print(f'BOOL 2-body: {arr}')
            #print(f'argmax: {np.argmax(arr)}')

            lowst_E_packable_idx = np.argmax(arr)
            irot = list_of_rotamers[i][lowst_E_packable_idx]
            score = scores[i][lowst_E_packable_idx]
            AA = rotamer_file[irot][0]

            #print(f'list of rotamers @ i  : {list_of_rotamers[i]}')
            #print(f'scores of rotamers @ i: {scores[i]}')
            #print(f'LOWEST E packable IDX : {lowst_E_packable_idx}')
            #print(f'score                 : {score}')
            #print(f'AA                    : {AA}')

            TOTAL_RIF_SCORE += score
            SEQ_MAP.append(
                [i, AA]
            )

        #print(f'TOTAL RIF SCORE: {TOTAL_RIF_SCORE}')

        if any_packable == 0: return 1000, None

        # get the SCORE of the best interactions
        return TOTAL_RIF_SCORE, SEQ_MAP

    def search_xyzs(
        self,
        xyzs,
        L_idx=None,
        D_idx=None,
        binder_idx=None,
    ):
        '''
        Given numpy array of xyzss (Mres X 3 X 3) of [Mres x (N CA C) x (X Y Z)]
        returns if there are hits that meet score thresholds, and the sequence mapping
        '''
        if L_idx is None: L_idx = np.arange(xyzs.shape[0])
        if D_idx is None: D_idx = np.array([])# assuming most people wont search D rifs
        if binder_idx is None: binder_idx = np.arange(xyzs.shape[0])


        # 1. get the hashed stubs of the stubs, based on L or D rif
        L_stubs = self.RIF_stubs(xyzs[L_idx,0,:], xyzs[L_idx,1,:], xyzs[L_idx,2,:])
        D_stubs = self.RIF_stubs(xyzs[D_idx,0,:], xyzs[D_idx,1,:], xyzs[D_idx,2,:])
        L_bins = self.L_AA_RIF['binner'].get_bin_index(L_stubs)
        D_bins = self.D_AA_RIF['binner'].get_bin_index(D_stubs)

        # 2. look up if there are any hits for L and D RIFs
        L_matching_keys = self.L_AA_RIF['rif_dict'].contains(L_bins)
        D_matching_keys = self.D_AA_RIF['rif_dict'].contains(D_bins)

        # get the irot offsets
        L_keys = L_bins[L_matching_keys]
        L_offsets = self.L_AA_RIF['rif_dict'][L_bins]
        D_keys = D_bins[D_matching_keys]
        D_offsets = self.D_AA_RIF['rif_dict'][D_bins]

        # 3. test if we have reached the thresholds for number of hits at this point
        if np.count_nonzero(L_offsets) + np.count_nonzero(D_offsets) <= self.min_RIF_hits - 1: return False, None, None

        # 4. get the specific irots for each hit positon
        rotamer_list, score_list, hotspot_count = self.get_irots(L_offsets, D_offsets)
        # test on hotspot
        if self.min_HOTSPOT_hits > 0 and hotspot_count < self.min_HOTSPOT_hits: return False, None, None

        # 5. get the L and D frames to rebuild the irots
        LD_frames = np.concatenate(
            (xyzs[L_idx,:,:], xyzs[D_idx,:,:]),
            axis=0,
        )

        # 6. pack
        assert len(rotamer_list) == len(score_list) == LD_frames.shape[0]
        total_rif_score, sequence_mapping = self.fast_pack(
            list_of_rotamers=rotamer_list,
            scores=score_list,
            hit_frames=LD_frames,
            rotamer_dofs=self.IROT_dict['irot_DOFS'],
            rotamer_file=self.IROT_dict['rotamer_lines'],
            macrocycle_bins=bin_coordinates(xyzs[binder_idx,:,:].reshape(-1, 3)[np.newaxis, ...]),# reshape to (1 x Natoms x 3)
            pert_mag=self.fast_pack_params['pert_mag'],
            shake_repeats=self.fast_pack_params['shake_repeats'],
            resl=self.fast_pack_params['resolution'],
            rotamer_clash_tolerance=self.fast_pack_params['rotamer_clash_tolerance'],
        )

        # 7. filter on total rif score & number of hits post pack
        if total_rif_score > self.min_RIF_score or len(sequence_mapping) < self.min_RIF_hits: return False, None, None

        # 8. convert sequence mapping to residue index
        LD_idx = np.concatenate((L_idx, D_idx), axis=0)
        for i in range(len(sequence_mapping)):
            sequence_mapping[i][0] = LD_idx[i] + 1# add 1 to get us to Rosetta 1-indexed arrays

        # 9. return total_score, sequence_mapping
        return True, total_rif_score, sequence_mapping

    def apply(self, pose):
        '''
        given a rosetta pose object, gets numpy arrays, searches RIF

        if pass, returns True, rif_score, sequence_mapping
        if fail, returns False, None, None

        '''
        # 1. align pose to the RIF target pose
        if self.L_AA_RIF['target'] is not None:
            pose = self.superimpose(self.L_AA_RIF['target'] , pose)
        elif self.D_AA_RIF['target'] is not None:
            pose = self.superimpose(self.D_AA_RIF['target'] , pose)
        else:
            raise ValueError('RIF target pose not defined')

        # 2. get indecies for L and D residues in the pose, and binder
        L_idx = np.squeeze(
            np.argwhere(
                np.logical_and(
                    self.residue_selector.apply(pose),
                    self.L_selector.apply(pose)
                ),
            ),
            axis=1,
        )
        D_idx = np.squeeze(
            np.argwhere(
                np.logical_and(
                    self.residue_selector.apply(pose),
                    self.D_selector.apply(pose),
                ),
            ),
            axis=1,
        )
        binder_idx = np.squeeze(
            np.argwhere(
                self.residue_selector.apply(pose),
            ),
            axis=1,
        )

        # 3. get all of the xyzs of atoms in atom_selector in the pose
        xyzs = np.array([[residue.atom(atom).xyz() for atom in self.atom_selector] for residue in pose.residues])

        # 4. search xyzs in RIF
        return self.search_xyzs(
            xyzs,
            L_idx=L_idx,
            D_idx=D_idx,
            binder_idx=binder_idx,
        )


def check_bins(
    bin_set1,
    bin_set2,
):
    """
    Given:
        two getpy.Set of byte8array3s

    returns:
        the number (int) of items from set2 that are set1

    Assumes that both getpy.Sets were created at the same resolution,
    otherwise the check makes no sense.

    It does not matter which order the sets are given
    """
    return np.sum(bin_set1.contains(bin_set2.items()))

def bin_coordinates(
    xyz_array,
    pert_mag=0.5,
    shake_repeats=10000,
    resl = 0.5,
):
    """
    Given:
        a numpy array of XYZ coordinates for some set of atoms formated as (1xNx3)
        (optionally) a magnitde in Angstroms to shake the coords as a float
        (optionally) number of "trajecotries" to shake the input coords as an int
        (optionally) resolution (in angstroms) to bin the cooridantes as a float

    returns:
        a getpy Set containing the bins of the coordinates as byte8array3
    """
    # 1. Add 0 vector into XYZ -> XYZ0
    expanded_xyzs = np.empty((xyz_array.shape[0], xyz_array.shape[1], 4))
    expanded_xyzs[:,:,:3] = xyz_array[...]
    expanded_xyzs[:,:,-1] = 0.0

    # 2. make some number of copies of the coords
    shaken_coords = np.repeat(expanded_xyzs, shake_repeats, axis=0)
    # Second. make some number of random perturbations as NxMx3 matrix
    random_perturbation = np.random.uniform(
        low=-pert_mag,
        high=pert_mag,
        size=(shake_repeats*shaken_coords.shape[1]*3)
    ).reshape(shake_repeats, shaken_coords.shape[1], 3)

    # 3. apply the random perturbations to each "trajectory"
    shaken_coords[:,:,:3] += random_perturbation

    # 4. round the coordinates as view as int, flatten on the first axis
    rounded_coords = np.around(shaken_coords / resl, decimals=0).astype(np.int32).flatten().view(np.dtype('S16'))
    # 5. set the coordinates in a getpy set
    coord_set = gp.Set(np.dtype('S16'))
    coord_set.add(rounded_coords)

    return coord_set

#!/usr/bin/env python
# Written by Brian Coventry


import sys
import os
import getpy as gp# i have no idea why this works. Deals with environment specific SegFault on import xbin.
import xbin
import h5py

import numpy as np

import argparse
from numba import njit
import numba




parser = argparse.ArgumentParser()
parser.add_argument("rif_h5", type=str)
parser.add_argument("-drop_sats", action="store_true", default=False)


args = parser.parse_args(sys.argv[1:])


@njit
def is_in(array, value):
    for i in range(len(array)):
        if ( array[i] == value ):
            return True

    return False



@njit
def merge_duplicate_keys(keys, scores, irots, sats):
    new_keys = np.zeros(keys.shape, np.int64 )
    new_scores = np.zeros(scores.shape, np.float32)
    new_irots = np.zeros(irots.shape, np.int16)
    new_sats = np.zeros(sats.shape, sats.dtype)
    new_irots.fill(-1)
    new_sats.fill(-1)

    nrots = scores.shape[1]

    total_keys = 0

    known_keys = numba.typed.Dict.empty(
        key_type=numba.types.int64,
        value_type=numba.types.int64,
    )

    for ikey in range(len(keys)):
        key = keys[ikey]
        score = scores[ikey]
        irot = irots[ikey]
        sat = sats[ikey]
        if ( key not in known_keys ):
            new_keys[total_keys] = key
            new_scores[total_keys] = score
            new_irots[total_keys] = irot
            new_sats[total_keys] = sat

            known_keys[key] = total_keys
            total_keys += 1
            continue

        new_ikey = known_keys[key]


        # we have to merge sort the two lists
        old_score = new_scores[new_ikey].copy()
        old_irot = new_irots[new_ikey].copy()
        old_sat = new_sats[new_ikey].copy()

        r = 0
        r_old = 0
        r_this = 0

        # Start filling until we run out of old or this
        old_out = False
        this_out = False

        while ( r < nrots ):

            score_old = old_score[r_old]
            score_this = score[r_this]

            if ( score_old < score_this ):

                if ( not is_in( new_irots[new_ikey, :r], old_irot[r_old] ) ):

                    new_scores[new_ikey, r] = score_old
                    new_irots[new_ikey, r] = old_irot[r_old]
                    new_sats[new_ikey, r] = old_sat[r_old]

                    r += 1
                r_old += 1
                if ( r_old == nrots or old_irot[r_old] == -1 ):
                    old_out = True
                    break

            else:

                if ( not is_in( new_irots[new_ikey, :r], irot[r_this] ) ):
                    new_scores[new_ikey, r] = score_this
                    new_irots[new_ikey, r] = irot[r_this]
                    new_sats[new_ikey, r] = sat[r_this]

                    r += 1
                r_this += 1
                if ( r_this == nrots or irot[r_this] == -1 ):
                    this_out = True
                    break

        # fill remainer with old if this_out
        while ( r < nrots and this_out and r_old < nrots and old_irot[r_old] != -1 ):

            if ( not is_in( new_irots[new_ikey, :r], old_irot[r_old] ) ):
                new_scores[new_ikey, r] = old_score[r_old]
                new_irots[new_ikey, r] = old_irot[r_old]
                new_sats[new_ikey, r] = old_sat[r_old]

                r += 1
            r_old += 1

        # fill remainer with this if old_out
        while ( r < nrots and old_out and r_this < nrots and irot[r_this] != -1):

            if ( not is_in( new_irots[new_ikey, :r], irot[r_this] ) ):
                new_scores[new_ikey, r] = score[r_this]
                new_irots[new_ikey, r] = irot[r_this]
                new_sats[new_ikey, r] = sat[r_this]

                r += 1
            r_this += 1


    return new_keys[:total_keys], new_scores[:total_keys], new_irots[:total_keys], new_sats[:total_keys]



@njit
def to_linear_format(scores, irots, sats):
    max_len = scores.shape[0] * (scores.shape[1] + 1)

    offsets = np.zeros(len(scores), np.int64 )
    new_scores = np.zeros(max_len, np.float32)
    new_irots = np.zeros(max_len, np.int16)
    new_sats = np.zeros((max_len, sats.shape[-1]), sats.dtype)
    new_irots.fill(-1)
    new_sats.fill(-1)

    # leave the first one blank
    cur_off = 1

    for ikey in range(len(scores)):
        offsets[ikey] = cur_off

        for j in range(scores.shape[1]):
            if ( irots[ikey, j] == -1 ):
                break
            new_scores[cur_off] = scores[ikey, j]
            new_irots[cur_off] = irots[ikey, j]
            new_sats[cur_off] = sats[ikey, j]

            cur_off += 1

        # add a -1
        cur_off += 1

    return new_scores[:cur_off], new_irots[:cur_off], new_sats[:cur_off], offsets


def size(item):
    return np.prod(item.shape) * item.dtype.itemsize


with h5py.File(args.rif_h5) as f:

    cart, ori, bound = list(f['cart_ori_bound'])

    print("Cart resl: %5.2f"%(cart))
    print("Angl resl: %5.2f"%(ori))
    print("Bound:     %5.1f"%(bound))

    # convert from rif to xbin
    cart = cart*2/np.sqrt(3)/2*1.5

    hasher = xbin.XformBinner(cart, ori, bound)

    flats = f['bin_center']
    xforms = np.zeros((len(flats), 4, 4))
    xforms[:,:3,:3] = flats[:,:9].reshape(-1, 3, 3)
    xforms[:,:3,3] = flats[:,9:]
    xforms[:,3,3] = 1

    print("Converting %i xforms to keys"%len(xforms))
    keys = hasher.get_bin_index(xforms)

    uniq, counts = np.unique(keys, return_counts=True)
    non_unique = (counts - 1).sum()

    print(" %i / %i keys are not unique (%.1f%%)"%(non_unique, len(keys), non_unique/len(keys)*100))

    print("Pulling remaining data")

    scores = f['scores'] * np.array(1.0, dtype=np.float32)
    irots = f['irots'] * np.array(1, dtype=np.int16)

    ## SATS
    num_sats = len([x for x in list(f) if x.startswith("sat")])
    sat_dtype = np.int8
    if ( args.drop_sats ):
        num_sats = 0
    if ( num_sats > 0 ):
        sat_dtype = f['sat0'].dtype
    sats = np.zeros(list(scores.shape) + [num_sats], dtype=sat_dtype)
    for isat in range(num_sats):
        sats[:,:,isat] = f['sat%i'%isat]

    print("Found %i sat slots with dtype: %s"%(num_sats, str(sat_dtype)))


    print("Merging duplicate keys")
    new_keys, new_scores, new_irots, new_sats = merge_duplicate_keys( keys, scores, irots, sats )
    old_size = size(new_keys) + size(new_scores)/2 + size(new_irots) + size(new_sats)

    print("Converting to linear format")
    linear_scores, linear_irots, linear_sats, offsets = to_linear_format(new_scores, new_irots, new_sats)
    assert(len(offsets) == len(new_keys))
    new_size = size(new_keys) + size(offsets) + size(linear_scores)/2 + size(linear_irots) + size(linear_sats)

    print("Final size: %.2fGB -- %i%%   (old size: %.2fGB)"%(new_size/1e9, new_size / old_size * 100, old_size/1e9))

    new_name = os.path.join(os.path.dirname(args.rif_h5), "py_" + os.path.basename(args.rif_h5))
    print("Writing to: " + new_name)

    with h5py.File(new_name, "w") as newf:

        newf.create_dataset('cart_ori_bound', data=np.array([cart, ori, bound]))
        newf.create_dataset('xbin_key', data=new_keys)
        newf.create_dataset('offsets', data=offsets)
        newf.create_dataset('scores', data=linear_scores.astype(np.float16))
        newf.create_dataset('irots', data=linear_irots)

        if ( num_sats > 0 ):
            newf.create_dataset('sats', data=linear_sats)


    print("")
    print("======================= Bin error analysis ===============================")


    new_centers = hasher.get_bin_center(keys)

    cart_error = np.linalg.norm(new_centers[:,:3,3] - xforms[:,:3,3], axis=-1)

    traces = np.trace( new_centers[:,:3,:3] @ np.linalg.inv(xforms[:,:3,:3]), axis1=-1, axis2=-2 )
    cos_theta = ( traces - 1 ) / 2
    angl_error = np.degrees( np.arccos( cos_theta.clip(-1, 1) ) )



    def histogram(data, label):
        num_bins = 60
        height = 10
        space = 10
        lb = 0
        ub = np.percentile(data, 99)

        bins = np.linspace(lb, ub, num_bins)

        assignments = np.searchsorted(bins, data).clip(0, num_bins-1)


        counts = np.zeros(num_bins)
        np.add.at(counts, assignments, 1)

        scaled_counts = np.rint(counts * height / counts.max()).astype(int)

        for i in range(height):
            j = height - i
            mask = scaled_counts >= j
            string = "".join("*" if x else " " for x in mask)

            print(" "*space + string)

        print(" "*space + "-"*num_bins)

        to_print = " "*space
        for i in range(num_bins//10 + 1):
            col = i * 10
            value = bins[i]

            to_add = "%-4g"%value
            to_add += " "*(10 - len(to_add))
            to_print += to_add

        print(to_print)

        print("")
        label_size = len(label) // 2
        print(" "*(space + num_bins//2 - label_size//2) + label)




    histogram(cart_error, "Cart Error (A)")
    histogram(angl_error, "Angle Error (deg)")








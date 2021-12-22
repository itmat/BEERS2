import numpy as np
cimport numpy as np
np.import_array()

def get_frac_skipped(float rate, int max_skips, int molecule_count, int read_len):
    ''' Return an array of length (read_len, max_skips+1) that indicate how many
    of the `molecule_count` molecules have incurred exactly `i` skips by the `j`th
    base in (i,j) element
    '''
    if rate <= 0:
        res = np.zeros((read_len, max_skips+1))
        res[:,0] = 1
        return res

    cdef long skip_start, skip_end, molecule, skip_num, i

    # skip_locs[i,j] is the base number at which the ith molecule makes its jth skip
    # if past the end of read_len, then 'never' skips
    # Geometric distribution gives the number of throws until the first 'success' in
    # a bias coin-toss, so this determines where the skips happen.
    # Note: faster than performing a coin toss for each base if the skip rate is small, as usual
    cdef np.ndarray[long, ndim=2] skip_locs = np.random.geometric(rate, size=(molecule_count, max_skips))
    skip_locs = np.cumsum(skip_locs, axis=1) - 1
    # num_skipped[k,j] is the number of molecules that have had exactly j skips by the kth base
    cdef np.ndarray[long, ndim=2] num_skipped = np.zeros((read_len, max_skips), dtype=int)
    for molecule in range(molecule_count):
        for skip_num in range(max_skips):
            # We need to increment the counts of molecules that
            # have received exactly 'skip_num' skips over the
            # range of bases that this molecule has had
            # exactly that many skips.
            # The start and end of this range is determined by skip_locs
            skip_start = skip_locs[molecule, skip_num]
            if skip_start > read_len - 1:
                # Early stop, skip occurs after the end of the read
                break

            if skip_num == max_skips - 1:
                skip_end = read_len
            else:
                skip_end = skip_locs[molecule, skip_num+1]
                if skip_end > read_len:
                    skip_end = read_len
                
            # Increment the appropriate range
            for i in range(skip_start, skip_end):
                num_skipped[i, skip_num] += 1

    frac_skipped = num_skipped / float(molecule_count)
    # Add on the no-skipped ones
    frac_skipped = np.hstack([
        (np.ones(frac_skipped.shape[0]) * (1 - frac_skipped.sum(axis=1)))[:, None],
        frac_skipped
    ])
    return frac_skipped

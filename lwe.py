def get_mask(sample_range):
    """
    Get the mask for the sample_noise algorithm given the sample range
    """
    ctr=0
    while sample_range>0:
        ctr+=1
        sample_range = sample_range >> 1
    return (1 << (ctr+1)) - 1
    
def sample_noise(length, max_hamming = None, seed=None, sample_elem = None):
    """
    Python port of create_spter_vec from Round2 submission

    INPUT:

    - ``length`` - the length of the vector to return
    - ``max_hamming`` - the maximum hamming weight allowed for the response vector; for a given value of ``max_hamming``, exactly (``length``-``max_hamming``) elements of the vector will be set to 0. The other ``max_hamming`` elements will be chosen via the method given in ``sample_elem`` so potentially the actual hamming weight of the output vector could be strictly less than ``max_hamming`` if ``sample_elem`` ever outputs 0
    - ``sample_elem`` - should be in the form ``sample_elem`` = (``sample_func``, ``sample_min_val``, ``sample_max_val``) and samples individual elements
    
    """
    if max_hamming is None or max_hamming < 0 or max_hamming > length:
        max_hamming = length
    if seed is not None:
        from sage.misc.randstate import set_random_seed
        set_random_seed(seed)
    if sample_elem is None or len(sample_elem) != 3:
        # set default as +-1 randomly
        def temp_samp():
            return 2*randrange(0,2)-1
        sample_func=temp_samp
        sample_min_val=-1
        sample_max_val=1
    else:
        sample_func=sample_elem[0]
        sample_min_val = sample_elem[1]
        sample_max_val = sample_elem[2]

    sample_range = sample_max_val - sample_min_val
    # get mask 2^{k}-1
    sample_mask = get_mask(sample_range)
    
    h_arr=[]
    h_arr.extend([sample_func() - sample_min_val for i in range(max_hamming)])
    h_arr.extend([-sample_min_val for i in range(length-max_hamming)])

    rnd_arr = [(randrange(0,1 << 32) & ~sample_mask) ^ (h_arr[i] & sample_mask)    for i in range(length)]
    rnd_arr.sort()

    ret_vec=[(rnd_arr[i] & sample_mask) + sample_min_val for i in range(length)]
    return ret_vec
    
    

    
        

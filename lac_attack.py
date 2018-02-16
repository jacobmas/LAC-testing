#!/usr/bin/env sage
 
import sys
from sage.all import *

""" Comment: radix sort an array of integers """

def sample_binom(size=1):
    return reduce(lambda x,y: x+y,[randrange(0,2) - randrange(0,2) for i in range(size)])
 
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
    
    h_arr=[sample_func() - sample_min_val for i in range(max_hamming)]
    h_arr.extend([-sample_min_val for i in range(length-max_hamming)])

    rnd_arr = [(randrange(0,1 << 32) & ~sample_mask) ^ (h_arr[i] & sample_mask)    for i in range(length)]
    rnd_arr.sort()

    ret_vec=[(rnd_arr[i] & sample_mask) + sample_min_val for i in range(length)]
    return ret_vec

def sample_message(n):
    ret_lst = [randint(0,1) for i in range(n)]
    return ret_lst

def decode_message(S,q,noisy_pol):
    noisy_lift=noisy_pol.lift().list()
    ret_msg=[]
    for i in range(len(noisy_lift)):
        x=ZZ(noisy_lift[i]) % q
#        print('x={0}'.format(x))
        while x < 0:
            x=x+q
        while x >= q:
            x=x-q
        if q/4. <= x and x < 3.*q/4.:
            ret_msg.append(1)
        else:
            ret_msg.append(0)
    return ret_msg
    

 
def test_correctness(n=512, c2_vec_num=None, q=251, seed=None, hamming_add=0):
    """
    Testing LAC PQC submission correctness values
    n = such that ring is x^{n}+1 (not checking for power of 2  though) 
    q = modulus

    Returns number of positions where decrypted message differs from encrypted message
    ``hamming_add`` - for doing maximum hamming weight type testing
 
    """
    if c2_vec_num is None:
        c2_vec_num = n
       
    if seed is not None:
        from sage.misc.randstate import set_random_seed
        set_random_seed(seed)
 
 
    ZZ = IntegerRing()
    ZZ_q = IntegerModRing(q)
    
 
 
 
    R = ZZ_q['x'].quotient(cyclotomic_polynomial(2*n, 'x'), 'x')

#    print("R={0},{1}".format(R,type(R)))
    bar_q=R((q-1)/2)
    
    # Set up public key 
    a_pol=R.random_element()
    lac_samp=(sample_binom,-1,1) # shorten command
    s_pol=R(sample_noise(length=n,max_hamming=n,sample_elem=lac_samp))
    e_pol=R(sample_noise(length=n,max_hamming=n,sample_elem=lac_samp))

    b_pol=(a_pol*s_pol+e_pol)


    # Set up encryption
    message = sample_message(R.degree())

    my_hamm=n/2+hamming_add
    
    r_pol = R(sample_noise(length=n,max_hamming=my_hamm))
    e1_pol = R(sample_noise(length=n,max_hamming=my_hamm))
    e2_pol = R(sample_noise(length=n,max_hamming=my_hamm))

    c1_pol = a_pol * r_pol + e1_pol
    c2_pol = b_pol * r_pol + e2_pol + bar_q * R(message)


    # Set up decryption
    u_pol = R(c1_pol * s_pol)

    noisy_pol = R(c2_pol - u_pol)

    decoded_message = decode_message(R,q,noisy_pol)

    def comp_mess(x,y):
        return 1 if x!=y else 0
    mult_message= [comp_mess(message[i],decoded_message[i]) for i in range(len(decoded_message))]
    diff_weight=reduce(lambda x,y: x+y, mult_message)

    return diff_weight



def get_fail_probs(runs=100000,n=512,c2_vec_num=None,q=251, hamming_add=0):
    fail_map={}
    for i in range(n+1):
        fail_map[i]=0
    for i in range(runs):
        y=test_correctness(n,c2_vec_num,q=q,hamming_add=hamming_add)
        fail_map[y]+=1
        if i % 1000==0:
            print("i={0}".format(i))
        if i % 5000 == 0:
            for j in range(n+1):
                if fail_map[j]>0:
                    print("{0}: {1}".format(j,fail_map[j]))            
    for i in range(n+1):
        if fail_map[i]>0:
            print("{0}: {1}".format(i,fail_map[i]))

def compute_failure(n,t,delta):
    result=0.0
    for i in range(t+1,n+1):
        curr_val=binomial(n,i)*pow(1-delta,n-i)*pow(delta,i)
        result=result+curr_val
    return result
    


 

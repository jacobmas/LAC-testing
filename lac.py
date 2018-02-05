#!/usr/bin/env sage
 
import sys
from sage.all import *
 
def sample_single_error():
    t=randint(0,3)
    if t==0:
        return -1
    elif t<3:
        return 0
    return 1
 
def sample_noise(S):
    n=S.degree()
    ret_lst=[sample_single_error() for i in range(n)]
    return S(ret_lst)

def sample_message(n):
    ret_lst = [randint(0,1) for i in range(n)]
    return ret_lst

def decode_message(S,q,noisy_pol):
    noisy_lift=noisy_pol.lift().list()
    ret_msg=[]
    for i in range(len(noisy_lift)):
        x=ZZ(noisy_lift[i]) % q
#        print('x={0}'.format(x))
        if x < 0:
            x=x+q
        if q/4. <= x and x < 3.*q/4.:
            ret_msg.append(1)
        else:
            ret_msg.append(0)
    return ret_msg
    

 
def test_correctness(n=512, q=251, seed=None):
    """
    Testing LAC PQC submission correctness values
    n = such that ring is x^{n}+1 (not checking for power of 2 yet though) 
    q = modulus

    Returns number of positions where things differ
 
 
    """
       
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

    s_pol=sample_noise(R)
    e_pol=sample_noise(R)

    b_pol=(a_pol*s_pol+e_pol)


    # Set up encryption
    message = sample_message(R.degree())

    r_pol = sample_noise(R)
    e1_pol = sample_noise(R)
    e2_pol = sample_noise(R)

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
#    print(message)
#    print(decoded_message)
#    print(mult_message)
    #print('message={0}\nDecoded message={1}\n'.format(message,decoded_message))
    return diff_weight

def test_correctness2(n=512, q=251, seed=None):
    """
    Testing LAC PQC submission correctness values
    n = such that ring is x^{n}+1 (not checking for power of 2 yet though) 
    q = modulus

    Returns number of positions where things differ
 
 
    """
    from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
    from sage.matrix.constructor import identity_matrix, block_matrix
    from sage.matrix.matrix_space import MatrixSpace
    from sage.rings.integer_ring import IntegerRing
    from sage.modules.free_module_integer import IntegerLattice
       
    if seed is not None:
        from sage.misc.randstate import set_random_seed
        set_random_seed(seed)
 
 
    ZZ = IntegerRing()
    ZZ_q = IntegerModRing(q)
    
 

    #R=PolynomialRing(Integers(q),'x',implementation="FLINT")
 
    R = ZZ_q['x'].quotient(cyclotomic_polynomial(2*n, 'x'), 'x')
    bar_q=R((q-1)/2)
    
    # Set up public key 
    a_pol=R.random_element()
#    print("type(a_pol)={0}, type(a_pol.lift())={1}".format(type(a_pol),type(a_pol.lift())))

    s_pol=sample_noise(R)
    e_pol=sample_noise(R)

    b_pol=R(a_pol.lift()*s_pol.lift()+e_pol.lift())


    # Set up encryption
    message = sample_message(R.degree())

    r_pol = sample_noise(R)
    e1_pol = sample_noise(R)
    e2_pol = sample_noise(R)

    c1_pol = R(a_pol.lift() * r_pol.lift() + e1_pol.lift())
    c2_pol = R(b_pol.lift() * r_pol.lift() + e2_pol.lift()) + bar_q * R(message)


    # Set up decryption
    u_pol = R(c1_pol.lift() * s_pol.lift())

    noisy_pol = R(c2_pol - u_pol)

    decoded_message = decode_message(R,q,noisy_pol)

    def comp_mess(x,y):
        return 1 if x!=y else 0
    mult_message= [comp_mess(message[i],decoded_message[i]) for i in range(len(decoded_message))]
    diff_weight=reduce(lambda x,y: x+y, mult_message)
#    print(message)
#    print(decoded_message)
#    print(mult_message)
    #print('message={0}\nDecoded message={1}\n'.format(message,decoded_message))
    return diff_weight

def get_fail_probs(runs=100000,n=512,q=251):
    fail_map={}
    for i in range(513):
        fail_map[i]=0
    for i in range(runs):
        y=test_correctness(n,q)
        fail_map[y]+=1
        if i % 1000==0:
            print("i={0}".format(i))
        if i % 5000 == 0:
            for j in range(513):
                if fail_map[j]>0:
                    print("{0}: {1}".format(j,fail_map[j]))            
    for i in range(513):
        if fail_map[i]>0:
            print("{0}: {1}".format(i,fail_map[i]))
    


 

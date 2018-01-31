

def create_spter_vec(len,h_weight,seed=None):
    h_arr=[2 if i % 2 == 0 else 0 for i in range(0,h_weight)]
    h_arr.extend([1 for i in range(h_weight,len)])
    rnd_arr=[randrange(0,65536) for i in range(0,len)]
    ret_vec=[(~0x3 & rnd_arr[i]) ^ (h_arr[i] & 0x3) for i in range(0,len)]
    bob=sorted(ret_vec)
    ret_vec = [(bob[i] & 0x3) - 1 for i in range(0,len)]
    return ret_vec

def sample_D1():
    CDF_TABLE=[190,437,504,511]
    rnd=randrange(0,512)
    sample=0
    sign=randrange(0,1)
    for i in range(0,len(CDF_TABLE)):
        sample = sample + ((CDF_TABLE[i] - rnd) >> 15)

    sample = ((-sign) ^ sample) + sign
    return sample
    

def sample_a(len, q):
    a=[randrange(0,q) for i in range(0,len)]
    return a

def sample_e(len, q):
    e=[sample_D1() for i in range(0,len)]
    return e


def compute_lizard(len=1024,q=1024,h_weight=128):
    T=PolynomialRing(ZZ,'x')
    R=PolynomialRing(IntegerModRing(q),'x')
    f1=x**1024+1
    f2=x**4+1
    s=R(create_spter_vec(len,h_weight)) % f1
    a=R(sample_a(len,q)) % f1
    e=R(sample_e(len,q)) % f1
    b=(a*s+e) % f1

    e_prime = ((s % f2)*(a % f2) - (b % f2)) % f2
    eideal=e % f2
    print eideal
    print T(eideal.list()) % 2
    print e_prime
    

    

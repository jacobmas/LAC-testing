R=PolynomialRing(ZZ,'x')

def create_spter_vec(len,h_weight):
    h_arr=[2 if i % 2 == 0 else 0 for i in range(0,h_weight)]
    h_arr.extend([1 for i in range(h_weight,len)])
    rnd_arr=[randrange(0,65536) for i in range(0,len)]
    ret_vec=[(~0x3 & rnd_arr[i]) ^ (h_arr[i] & 0x3) for i in range(0,len)]
    bob=sorted(ret_vec)
    ret_vec = [(bob[i] & 0x3) - 1 for i in range(0,len)]
    return ret_vec

    

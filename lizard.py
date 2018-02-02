
def create_spter_vec(len,h_weight,seed=None):
    h_arr=[2 if i % 2 == 0 else 0 for i in range(0,h_weight)]
    h_arr.extend([1 for i in range(h_weight,len)])
    rnd_arr=[randrange(0,65536) for i in range(0,len)]
    ret_vec=[(~0x3 & rnd_arr[i]) ^ (h_arr[i] & 0x3) for i in range(0,len)]
    bob=sorted(ret_vec)
    ret_vec = [(bob[i] & 0x3) - 1 for i in range(0,len)]
    return ret_vec

def sample_D1(to_print=False):
    CDF_TABLE=[190,437,504,511]
    rnd=randrange(0,512)
    sample=0
    sign=randrange(0,2)
    if to_print:
        print('rnd={0}, sign={1}'.format(rnd,sign))

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

def get_odd_count(to_count):
    odd_ct=reduce(lambda x,y: x+y, list(map(lambda x: abs(x % 2),to_count)))
    return odd_ct
def compute_lizard(len=1024,sub_len=4,q=1024,h_weight=128,c=1):
    T=ZZ['x']
    R=PolynomialRing(ZZ,'x')
    I=R.ideal([2,(x-1**4)])

    S=QuotientRing(R,I*I)
    f1=R(x**len+1)
    f2=R(2*(x**sub_len-1))
    f3 = R(x**(2*sub_len)-2*x**sub_len+1)
    s=R(create_spter_vec(len,h_weight)) % f1
    print("f2={0},f3={1}".format(f2,f3))
    s=s % q
    a=R(sample_a(len,q)) % f1
    a = a % q
    e=R(sample_e(len,q)) % f1
    e = e % q
    b=(a*s+e) % f1
    b= b % q
    #print('a,b,s,e={0},{1},{2},{3}'.format(a,b,s,e))
    e_prime_ideal = (-(s % f3 % f2)*(a %  f3 % f2) + (b % f3 % f2)) %  f3 % f2
    e_prime_ideal = e_prime_ideal % q 
    eideal=e % f3 % f2
    eideal = eideal % q
    sideal=s  % f3 % f2
    sideal = sideal % q
#    emodtwo=R(((T(eideal.list())*f3) % 2).list())
#    smodtwo=R(((T(sideal.list())*f3) % 2).list())    
    ideal_diff=(eideal-e_prime_ideal) % f2 % f3
    ideal_diff = ideal_diff % q
    s_odd=get_odd_count(T(sideal).list())
    e_odd=get_odd_count(T(eideal).list())
    dictout={}
    dictout['s']=[s_odd,sub_len-s_odd]
    dictout['e']=[e_odd,sub_len-e_odd]
#    print("s: odd={0}, even={1}".format(s_odd,sub_len-s_odd))
#    print("e: odd={0}, even={1}".format(e_odd,sub_len-e_odd))

    print("ideal_diff={0}".format(ideal_diff))
    print ("e % f2 = {0}".format(eideal))
    print ("-(s % f2)*(a % f2) + (b % f2)={0}".format(e_prime_ideal))
    return dictout

def reducebygens(a,gens):
    for y in gens:
        a=a%y
    for y in gens[0:1]:
        a=a%y
    return a

def compute_lizard2(len=1024,sub_len=4,q=1024,h_weight=128,c=1):
    T=ZZ['x']
    R=PolynomialRing(ZZ,'x')
    I=R.ideal([2,(x**4-1)])
    I2=I*I
    S=QuotientRing(R,I2)
    print S
    f1=R(x**len+1)
    f2=R(2*(x**sub_len-1))
    f3 = R(x**(2*sub_len)-2*x**sub_len+1)
    s=R(create_spter_vec(len,h_weight)) % f1
    print("f2={0},f3={1}".format(f2,f3))
   # s=s % q
    a=R(sample_a(len,q)) % f1
    I2star=[S(f) for f in I2.gens()]
    print(I2.gens()[0])
    y=I2.gens()[0]
    
    print("a={0}".format(a))
    print("I2star={0}".format(I2star))
    Sa=S(a)
    a=reducebygens(a,I2.gens())
    print("Sa={0}".format(a))
#    \nS(a)={1}".format(a,a % y ))

   # a = a % q
    e=R(sample_e(len,q)) % f1
   # e = e % q
    b=(a*s+e) % f1
   # b= b % q
    #print('a,b,s,e={0},{1},{2},{3}'.format(a,b,s,e))
    e_prime_ideal = S(b-s*a)
    #e_prime_ideal = e_prime_ideal % q 
    eideal=S(e)
   # eideal = eideal % q
#    sideal=s  % f3 % f2
   # sideal = sideal % q
#    emodtwo=R(((T(eideal.list())*f3) % 2).list())
#    smodtwo=R(((T(sideal.list())*f3) % 2).list())    
    ideal_diff=(eideal-e_prime_ideal)
#    ideal_diff = ideal_diff % q
  #  s_odd=get_odd_count(T(sideal).list())
  #  e_odd=get_odd_count(T(eideal).list())
 #   dictout={}
 #   dictout['s']=[s_odd,sub_len-s_odd]
 #   dictout['e']=[e_odd,sub_len-e_odd]
#    print("s: odd={0}, even={1}".format(s_odd,sub_len-s_odd))
#    print("e: odd={0}, even={1}".format(e_odd,sub_len-e_odd))

   # print("ideal_diff={0}".format(ideal_diff))
   # print ("e,e % f2 = {0},{1}".format(e,eideal))
   # print ("-(s % f2)*(a % f2) + (b % f2)={0}".format(e_prime_ideal))
#    return dictout

def many_lizard(count=1000,sub_len=64):
    total_dict={'s': [0,0], 'e': [0,0]}
    for i in range(count):
        if i % 1000 == 0:
            print('i={0}'.format(i))
        ret=compute_lizard(sub_len=sub_len)
        for j in range(2):
            total_dict['s'][j]+=ret['s'][j]
            total_dict['e'][j]+=ret['e'][j]
    return total_dict

def sample_many_D1(num):
    dist={}
    for i in range(0,num):
        temp=sample_D1()
        if temp not in dist:
            dist[temp]=1
        else:
            dist[temp]=dist[temp]+1
    return dist
    

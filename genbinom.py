#import math
#import sys
#from sage.all import *

def binom(n,k):
    return int(math.factorial(n)/(math.factorial(k)*math.factorial(n-k)))
def do_stuff(coeffs,vardegs):
    result=1
    #print '{0},{1}'.format(coeffs,vardegs)
    for c in zip(coeffs,vardegs):
        #print '{0}, {1}'.format(c[0],c[1])
        result=result*(c[0]**c[1])
    #print 'result={0}'.format(result)
    return result
def gen_binom(p, num_var, deg,coeffs, the_vec):
    vardegs=[0 for x in range(num_var)]
    # Should assert that num vars is non-neg
    tot_deg=0
    curr_pos=num_var-1
    counter=0
    #print '{0}: {1}'.format(counter,vardegs)
    while 1:
        
        is_good=1
        ## for i in vardegs:
        ##     if i>=p:
        ##         is_good=0
        if is_good:
            print '{0}: {1}'.format(counter,vardegs)
            the_vec.append(do_stuff(coeffs,vardegs)) # append the result of taking the coeffs to the degrees specified by vardegs

        
        if vardegs[0]>=deg:
            break

        counter+=1
        if tot_deg<deg:
            curr_pos=num_var-1
        else:
            tot_deg-=vardegs[curr_pos]
            vardegs[curr_pos]=0
            curr_pos-=1
        vardegs[curr_pos]+=1
        tot_deg+=1
   
        #print '{0}: {1}'.format(counter,vardegs)
    result=binom(num_var+deg,deg)
    #print 'counter={0}, (var+deg)^C_(deg)={1}'.format(counter,result)

def test_smart_deg(p, q, deg):
    k=Mod(p,q).multiplicative_order()
    print 'k={0}'.format(k)
    F=FiniteField(p**k,'x')
    g=F.primitive_element() # The generator
    vec_list=[]
    b_list=[]
    for a in range(0,q):
        b_list.append(Mod(a,p))
        h=g**a
        #print h   
        coeffs=h.polynomial().list()
        while len(coeffs)<k:
            coeffs.append(0)
        #print '{0}, {1}, {2}'.format(h,coeffs,type(coeffs[0]))
        my_vec=[]
        gen_binom(p,k,deg,coeffs,my_vec)
        vec_list.append(my_vec)
            #if __name__ == "__main__":
    my_mat=matrix(GF(p),vec_list)
    my_vec=vector(GF(p),b_list)
    #print '\n\n'
#    print my_mat
    print 'Rank={0}, Columns={3}, q={1}, deg={2}'.format(my_mat.rank(),q,deg,my_mat.ncols())
    if my_mat.rank()==q:
        x=my_mat.solve_right(my_vec)
        print 'x={0}'.format(x)
        print 'my_vec={0}'.format(my_vec)
        print '{0}'.format(my_mat)
        print 'Success for degree={0}, p={1}, q={2}, (E+d-1)^C_(d-1)={3} (E+d)^C_d={4}'.format(deg,p,q,binom(k+deg-1,deg-1),binom(k+deg,deg))
        return 1
    else:
        return 0
def test_smart(p,q):
    deg=0
    result=0
    while not result:
        deg+=1
        result=test_smart_deg(p,q,deg)
    return
#    import sys
    # input modulus p, exponent k so q=p^k-1
#test_smart(int(sys.argv[1]),int(sys.argv[2]))

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
 
def my_gen_lattice2(n=4, q=11, seed=None,
                quotient=None, dual=False, ntl=False, lattice=False, GuessStuff=True):
    """
    This is a modification of the code for the gen_lattice function from Sage
 
    Randomness can be set either with ``seed``, or by using
    :func:`sage.misc.randstate.set_random_seed`.
 
    INPUT:
 
    - ``type`` -- one of the following strings
        - ``'cyclotomic'`` -- Special case of ideal. Allows for
          efficient processing proposed by [LM2006]_.
    - ``n`` -- Determinant size, primal:`det(L) = q^n`, dual:`det(L) = q^{m-n}`.
      For ideal lattices this is also the degree of the quotient polynomial.
    - ``m`` -- Lattice dimension, `L \subseteq Z^m`.
    - ``q`` -- Coefficient size, `q-Z^m \subseteq L`.
    - ``t`` -- BKZ Block Size
    - ``seed`` -- Randomness seed.
    - ``quotient`` -- For the type ideal, this determines the quotient
      polynomial. Ignored for all other types.
    - ``dual`` -- Set this flag if you want a basis for `q-dual(L)`, for example
      for Regev's LWE bases [Reg2005]_.
    - ``ntl`` -- Set this flag if you want the lattice basis in NTL readable
      format.
    - ``lattice`` -- Set this flag if you want a
      :class:`FreeModule_submodule_with_basis_integer` object instead
      of an integer matrix representing the basis.
 
    OUTPUT: ``B`` a unique size-reduced triangular (primal: lower_left,
      dual: lower_right) basis of row vectors for the lattice in question.
 
    EXAMPLES:
 
 
 
    Cyclotomic bases with n=2^k are SWIFFT bases::
 
        sage: sage.crypto.gen_lattice(type='cyclotomic', seed=42)
        [11  0  0  0  0  0  0  0]
        [ 0 11  0  0  0  0  0  0]
        [ 0  0 11  0  0  0  0  0]
        [ 0  0  0 11  0  0  0  0]
        [ 4 -2 -3 -3  1  0  0  0]
        [ 3  4 -2 -3  0  1  0  0]
        [ 3  3  4 -2  0  0  1  0]
        [ 2  3  3  4  0  0  0  1]
 
    Dual modular bases are related to Regev's famous public-key
    encryption [Reg2005]_::
 
        sage: sage.crypto.gen_lattice(type='modular', m=10, seed=42, dual=True)
        [ 0  0  0  0  0  0  0  0  0 11]
        [ 0  0  0  0  0  0  0  0 11  0]
        [ 0  0  0  0  0  0  0 11  0  0]
        [ 0  0  0  0  0  0 11  0  0  0]
        [ 0  0  0  0  0 11  0  0  0  0]
        [ 0  0  0  0 11  0  0  0  0  0]
        [ 0  0  0  1 -5 -2 -1  1 -3  5]
        [ 0  0  1  0 -3  4  1  4 -3 -2]
        [ 0  1  0  0 -4  5 -3  3  5  3]
        [ 1  0  0  0 -2 -1  4  2  5  4]
 
 
    """
    from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
    from sage.matrix.constructor import identity_matrix, block_matrix
    from sage.matrix.matrix_space import MatrixSpace
    from sage.rings.integer_ring import IntegerRing
    from sage.modules.free_module_integer import IntegerLattice
       
    if seed is not None:
        from sage.misc.randstate import set_random_seed
        set_random_seed(seed)
 
 
    m=n+1
    ZZ = IntegerRing()
    ZZ_q = IntegerModRing(q)
 
 
 
    from sage.arith.all import euler_phi
    from sage.misc.functional import cyclotomic_polynomial
 
    # we assume that n+1 <= min( euler_phi^{-1}(n) ) <= 2*n
    found = False
    for k in range(2*n,n,-1):
        if euler_phi(k) == n:
            found = True
            break
        if not found:
            raise ValueError("cyclotomic bases require that n "
                                 "is an image of Euler's totient function")
    R = ZZ_q['x'].quotient(cyclotomic_polynomial(2*n, 'x'), 'x')
    g=x**(n/2)+1
    T=ZZ_q['x'].quotient(x**(n/2)+1)
 
   
    a_pol=R.random_element()

    s_pol=sample_noise(R)
    e_pol=sample_noise(R)
 
    s_pol2=T((s_pol))
    e_pol2=T((e_pol))
    print("s={0},e={1}".format(T(s_pol),T(e_pol)))

    
    Z_mat=e_pol2.matrix().augment(s_pol2.matrix())
    Z_mattop=Z_mat[0:1].augment(matrix(1,1,[ZZ.one()*-1]))
  
  
    b_pol=(a_pol*s_pol+e_pol)
    print("s_pol={0}\ne_pol={1}".format((s_pol2).list(),(e_pol2).list()))
    # Does a linear mapping change the shortest vector size for the rest?/
    a_pol=a_pol#*x_pol
    b_pol=b_pol#*x_pol

    a_pol2 = T(a_pol.list())# % S(g.list())
    b_pol2 = T((b_pol).list())# % S(g.list())
#    print("a={0}\nb={1}".format(a_pol2,b_pol2))
    A=identity_matrix(ZZ_q,n/2)
    A=A.stack(a_pol2.matrix())
    
    
    b_prime=b_pol2.matrix()[0:1]
    b_prime=b_prime - 11*A[8:9]
    A=A.stack(b_pol2.matrix()[0:1])

    
    
#    print("X=\n{0}".format(X))

#    A = A.stack(identity_matrix(ZZ_q, n/2))
 
    print("A=\n{0}\n".format(A))
    # switch from representatives 0,...,(q-1) to (1-q)/2,....,(q-1)/2
    def minrepnegative(a):
        if abs(a-q) < abs(a): return (a-q)*-1
        else: return a*-1
    def minrep(a):
        if abs(a-q) < abs(a): return (a-q)
        else: return a
    A_prime = A[(n/2):(n+1)].lift().apply_map(minrep)
#    b_neg= A[(n):(n+1)].lift().apply_map(minrepnegative)
    Z_fixed=Z_mattop.lift().apply_map(minrep)
    print("Z_fixed={0}\n||Z_fixed||={1}".format(Z_fixed,float(Z_fixed[0].norm())))
    print('Z_fixed*A={0}\n\n'.format(Z_fixed*A))

    print("z_fixed[0].norm()={0}".format(float(Z_fixed[0].norm())))
#    B=block_matrix([[ZZ(q),ZZ.zero()],[A_neg,ZZ.one()]], subdivide=False)
#    B = block_matrix([[ZZ.one(), -A_prime.transpose()],
#                     [ZZ.zero(), ZZ(q)]], subdivide=False)
    B = block_matrix([[ZZ(q), ZZ.zero()], [-A_prime, ZZ.one()]], subdivide=False)
#    for i in range(m//2):
#        B.swap_rows(i,m-i-1)
    #    print("{0}\n".format(A_neg))
 #   B=block_matrix([[ZZ(q), ZZ.zero(),ZZ.zero()],[ZZ.one(),A_neg,ZZ.zero() ],[ZZ.zero(),b_neg,ZZ.one()]],
                     #  subdivide=False)
    #print("B=\n{0}".format(B))
    print("B*A=\n{0}\n\n".format(B*A))
    #print("A=\n{0}\n".format(A))
    def remap(x):
        return minrep((x*251)%251)
    BL=B.BKZ(block_size=n/2.)
    y=(BL.solve_left(Z_fixed))#.apply_map(remap))

#   print("y*B={0}".format(y*B))
    print("y:=B.solve_left(Z_fixed)={0}".format(y))
#    BL=B.BKZ(block_size=n/2.)
    print(BL[0])
    print("shortest norm={0}".format(float(BL[0].norm())))
#    L = IntegerLattice(B)
#    p
#    v=L.shortest_vector()
#    print("L.shortest_vector={0}, norm={1}".format(v,float(v.norm())))
    if ntl and lattice:
        raise ValueError("Cannot specify ntl=True and lattice=True ")
    if ntl:
        return B._ntl_()
    elif lattice:
        from sage.modules.free_module_integer import IntegerLattice
        return IntegerLattice(B)
    else:
        return B

#!/usr/bin/env sage
 
import sys
from sage.all import *
 

 
def do_factor(m=4, q=11):
    """

    Factor cyclotomic polynomial modulo q
    
 
    INPUT:
 
    - ``m`` -- factoring the mth cyclotomic polynomial = 
    - ``q`` -- modulus
 
    OUTPUT: factors of x
 

    """
    from sage.rings.finite_rings.integer_mod_ring import IntegerModRing
    from sage.matrix.constructor import identity_matrix, block_matrix
    from sage.matrix.matrix_space import MatrixSpace
    from sage.rings.integer_ring import IntegerRing
    from sage.modules.free_module_integer import IntegerLattice
       
 
    m=n+1
    ZZ = IntegerRing()
    ZZ_q = IntegerModRing(q)
 
 
 
    from sage.arith.all import euler_phi
    from sage.misc.functional import cyclotomic_polynomial
 
    R = ZZ_q['x']#.quotient(cyclotomic_polynomial(k, 'x'), 'x')
    f=cyclotomic_polynomial(m,'x')
    return f.factor()

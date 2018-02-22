class RingGramSchmidt:
    def __init__(self, the_basis):
        self.the_basis=the_basis # a matrix
        self.gram_schmidt=the_basis
        self.n=the_basis.ncols()
        self.coeffs = {}
        for j in range(0,self.n-1):
            for i in range(j+1,self.n):
                self.coeffs[(j,i)]=0
        self.update_gramschmidt()
            
    '''
    Update the Gram-Schmidt here
    '''
    def update_gramschmidt(self):
        for i in range(0,self.n-1):
            for j in range(0,i):
                # compute projections and stuff
                self.coeffs[(j,i)]=self.compute_gscoeff(self.the_basis.column(j),self.the_basis       
        
        
    

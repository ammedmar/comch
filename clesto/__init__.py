from collections import Counter
from itertools import combinations, product, chain, permutations

#_________________________________79_characters________________________________

def partitions(n, k, smallest_value=1, largest_value=None, ordered=False):
    '''n is the integer to partition and k is the length of partitions.
    It returns all k tuples of integers greater or equal to smallest_value 
    and less than or equal to largest_value that add to to n. 
    If ordered == True it returns all tuples if False it returns those 
    in non-decreassing order '''
    if largest_value == None:
        largest_value = n
    
    def unordered_partitions(n,k,l=smallest_value, m=largest_value):
        if k == 1:
            if l <= n <= m:
                yield (n,)
            return
        for i in range(l,m+1):
            for result in unordered_partitions(n-i,k-1,i,m):
                yield (i,)+result
                
    if ordered:
        return chain.from_iterable(set(permutations(p)) for p 
                                   in unordered_partitions(n,k))
    if not ordered:
        return unordered_partitions(n,k)

#_________________________________79_characters________________________________

class Module_element(Counter):
    '''Implements modules over the integers or the integers mod p'''

    torsion = None 

    def __init__(*args, **kwds):
        '''...'''
        # print('initializing as Module_element')
        self, *args = args
        super(Module_element, self).__init__(*args, **kwds)
        
        if not all( [type(v) is int for v in self.values()] ):
            raise TypeError('values must be integers')

        self.reduce_rep()

    def __str__(self):
        '''...'''
        if not self:
            return '0'
        else:
            answer = '' 
            for key, value in self.items():
                if value < -1:
                    answer += f'{value}{key}'
                elif value == -1:
                    answer += f'-{key}'
                elif value == 1:
                    answer += f'+{key}'
                elif value > 1:
                    answer += f'+{value}{key}'    
            if answer[0] == '+':
                answer = answer[1:]

            return answer

    def __add__(self, other):
        '''...'''
        answer = Module_element(self)
        answer.update(other)
        return type(self)( {k:v for k,v in answer.items() if v} )
    
    def __sub__(self, other):
        '''...'''
        answer = Module_element(self)
        answer.subtract(other)
        return type(self)( {k:-v for k,v in answer.items() if v} )

    def __mod__(self, p):
        '''...'''
        return type(self)( {k:v%p for k,v in self.items() if v % p} )
    
    def __neg__(self):
        '''...'''
        return type(self)( {k:-v for k,v in self.items() if v} )
    
    def __rmul__(self, c):
        '''...'''
        answer = Module_element()
        if type(c) is int:
            for key, value in self.items():
                answer[key] = value*c
            return type(self)(answer)
        else:
            raise TypeError(f"can't act by non-int of type {type(c)}")

    def __iadd__(self, other):
        '''...'''
        self.update(other)
        self.reduce_rep()
        return self
    
    def __isub__(self, other):
        '''...'''
        self.subtract(other)
        self.reduce_rep()
        return self
    
    def __imod__(self, p):
        '''...'''
        for key, value in self.items():
            self[key] = value % p
            
        return self

    def reduce_rep(self):
        '''chooses preferred representatives of values (if torsion 
        is specified) and deletes keys with 0 value'''
        # print('reducing as Module_element')
        if self.torsion:
            for key, value in self.items():
                self[key] = value % self.torsion

        zeros = [k for k,v in self.items() if not v]
        for key in zeros:
            del self[key]

# Module_element.torsion = 7
# a = Module_element({'a':6})
# print(a)

#_________________________________79_characters________________________________

class DGModule_element(Module_element):
    '''...'''
        
    def __str__(self):
        string = super().__str__()
        return string.replace(', ', ',')

    def reduce_rep(self):
        '''deletes degenerate keys and reduces as Module_element'''
        # print('reducing as DGModule_element')
        if not all([isinstance(x,tuple) for x in self.keys()]):
            raise TypeError('keys must be tuples')

        for simplex,v in self.items():
            if not type(simplex) is tuple:
                raise TypeError('keys must be tuples')
            
            for i in range(len(simplex)-1):
                if simplex[i] == simplex[i+1]:
                    self[simplex] = 0
                
        super().reduce_rep()
    
    def boundary(self):
        '''...'''
        bdry = DGModule_element()
        for spx, coeff in self.items():
            for i in range(len(spx)):
                i_face = tuple(spx[:i]+spx[i+1:])
                i_coeff = coeff*((-1)**i)
                bdry += DGModule_element({i_face: i_coeff})
        return bdry

# Module_element.torsion = 4
# x = DGModule_element({(1,1):3})
# print(x)

# #_________________________________79_characters________________________________

class Cyclic_Module_element(Module_element):
    '''Modeling elements in Z[C] or Z/nZ[C_n] where C is the 
    infinite cyclic group and C_n the n-torsion cyclic group 
    generated by abd element denoted by a'''
    
    psi_dict = {}

    def __str__(self):
        '''...'''
        if not self:
            return '0'
        else:
            answer = '' 
            for exponent, coefficient in self.items():
                if coefficient != 1 or exponent == 0:
                    answer += f'+{coefficient}a^{exponent}'
                else:
                    answer += f'+a^{exponent}'    
            if answer[0] == '+':
                answer = answer[1:]

            return answer.replace('a^0','').replace('a^1','a')
            
    def __mul__(self, other):
        '''...'''
        answer = Counter()
        for k1,v1 in self.items():
            for k2,v2 in other.items():
                answer[k1+k2] += v1*v2
        
        return Cyclic_Module_element(answer)
        
    def __call__(self, other):
        '''...'''
        if isinstance(other, Cyclic_Module_element):
            return self*other
        if isinstance(other, Cyclic_DGModule_element):
            answer = Counter()
            for k,v1 in self.items():
                for x,v2 in other.items():
                    y = tuple(k+i for i in x)
                    answer[y] += v1*v2
            return Cyclic_DGModule_element(answer)

    def reduce_rep(self):
        '''in place mod p reduction of the keys'''
        # print('reducing as Cyclic_Module_element')
        if not all([isinstance(k,int) for k in self.keys()]):
            raise TypeError('keys must be integers')
        if self.torsion:    
            aux = list(self.items())
            self.clear()
            for k,v in aux:
                self[k%self.torsion] += v
        
        super().reduce_rep()

    def psi(self, n):
        '''considering self as an element in W_n we compute 
        its image in E(Z/p[C_p]). We use that the coefficients
        in psi[n] are equal to 1 for each basis element'''

        try:
            answer = Cyclic_DGModule_element()
            for k1 in (Cyclic_Module_element.psi_dict[n]).keys():
                for k2,v2 in self.items():
                    to_add = Cyclic_DGModule_element({tuple(k2+i for i in k1): v2})
                    answer += to_add
            return answer

        except KeyError:
            Cyclic_Module_element._compute_psi(n)
            return self.psi(n)
    
    def transpo_element():
        return Cyclic_Module_element({1:1, 0:-1})

    def norm_element():
        if Cyclic_Module_element.torsion:
            return Cyclic_Module_element({i:1 for i in 
                                   range(Cyclic_Module_element.torsion)}) 

    @staticmethod
    def all_elements():
        '''...'''
        p = Cyclic_Module_element.torsion
        if p:
            return ( Cyclic_Module_element( dict(zip(range(p), coeffs)) ) 
                       for coeffs in product(range(p),repeat=p) )

    @classmethod
    def _compute_psi(self, i):
        '''...'''

        try:
            previous_psi = Cyclic_Module_element.psi_dict[i-1]

            operators = {0: Cyclic_Module_element.norm_element(),
                         1: Cyclic_Module_element.transpo_element()}

            op_psi = operators[i%2](previous_psi)
            
            Cyclic_Module_element.psi_dict[i] = ( 
                             Cyclic_DGModule_element(
                                {(0,) + k:v for k,v in op_psi.items()}) )
        
        except KeyError:
            if i == 0:
                Cyclic_Module_element.psi_dict[0] = \
                        Cyclic_DGModule_element({(0,):1})
            else:            
                Cyclic_Module_element._compute_psi(i-1)
                Cyclic_Module_element._compute_psi(i)

# Module_element.torsion = 3
# print(Cyclic_Module_element({3:4, 4:1}))

#_________________________________79_characters________________________________

class Cyclic_DGModule_element(DGModule_element):
    '''...'''

    def __str__(self):
        '''...'''
        s = super().__str__()
        s = s.replace(', ', ',')
        return s.replace('(','a^(')

    def reduce_rep(self):
        '''reduces mod p the keys and values and deletes keys with 0 value 
        or which are degenerate'''
        # print('reducing as Cyclic_DGModule_element')
        if self.torsion:    
            aux = list(self.items())
            self.clear()
            for x, v in aux:
                y = tuple(i % self.torsion for i in x)
                self[y] += v
        super().reduce_rep()

    def phi(self):
        p = self.torsion
        answer = Counter()
        for k, v in self.items():
            x = []
            for i in k:
                x.append(tuple(j%p + 1 for j in range(i, p+i)))
            answer[tuple(x)] = v
            
        return Symmetric_DGModule_element(answer)

# Module_element.torsion = 3
# x = Cyclic_Module_element({0:1})
# print('after x')
# print(str(x.psi(3)))
# print(Cyclic_Module_element.psi_dict)

# #_________________________________79_characters________________________________

class Symmetric_DGModule_element(DGModule_element):
    '''...'''
    
    def reduce_rep(self):
        # print('reducing as Symmetric_DGModule_element')
        for spx in self.keys():
            if not ( {len(sigma) for sigma in spx} 
                == {Module_element.torsion} ):
                raise TypeError(f'length of all tuples in {spx} must '
                                 + f'be {Module_element.torsion}')
        super().reduce_rep()

    def table_reduction(self):
        '''given a set of basis element in the Barratt-Eccles operad, it returns 
        the set of surjections in its image via the table reduction morphism'''
        
        answer = Counter()
        for bar_ecc_element, value in self.items():
            d, a = len(bar_ecc_element)-1, max(bar_ecc_element[0]) #dim and arity
            for pi in partitions(d+a, d+1, ordered=True): 
                surjection, removed = [], []
                degenerate = False
                for idx, i in enumerate(pi):
                    filtered =  [i for i in bar_ecc_element[idx] 
                                         if i not in removed]
                    if idx > 0 and surjection[-1] == filtered[0]:
                        degenerate = True
                        break    
                    if i > 1:
                        removed += filtered[:i-1]
                    surjection += filtered[:i]

                if not degenerate:
                    answer += Counter({tuple(surjection):value})
     
        return DGModule_element(answer)


# Module_element.torsion = 3
# x = Symmetric_DGModule_element({((1,2,3),(1,3,2)):2})
# print(x.table_reduction() + x.table_reduction())


# Module_element.torsion = 3
# x = Cyclic_DGModule_element({(0, 1, 2): 1, (0, 2, 3): 1})
# print(x.phi().table_reduction())



# Not needed but aesthetically appealing
# class Symmetric_Module_element(Z_p_Module_element):
#     '''...'''

#     def __init__(*args, **kwds):
#         '''...'''
#         self, *args = args
#         super(Z_pS_p_element, self).__init__(*args, **kwds)
#         p = Z_p_Module_element.prime
#         if bool(self) and (
#             not {len(k) for k in self.keys()} == {p} or \
#             not {frozenset(k) for k in self.keys()} \
#             == {frozenset(range(1,p+1))} ):
#                 raise TypeError(f'keys must be permutations \
#                                   of {tuple(range(1,p+1))}')
#         self.reduce_rep()

#     def __add__(self, other):
#         '''...'''
#         self.check_prime(other)
#         return Z_pS_p_element( super().__add__(other) )
    
#     def __sub__(self, other):
#         '''...'''
#         self.check_prime(other) 
#         return Z_pS_p_element( super().__sub__(other) )
    
#     def __neg__(self):
#         '''...'''
#         return Z_pS_p_element( {k:-v for k,v in self.items() if v} )
            
#     def __mul__(self, other):
#         '''...'''
#         self.check_prime(other)
#         answer = Z_pS_p_element()
#         for k1,v1 in self.items():
#             for k2,v2 in other.items():
#                 answer[tuple(k1[i-1] for i in k2)] += v1*v2
#         answer.reduce_rep()
#         return answer
    
#     def __call__(self, other):
#         '''...'''
#         if isinstance(other, Z_pS_p_element):
#             return self*other
#         if isinstance(other, EZ_pS_p_element):
#             self.check_prime(other) 
#             answer = Counter()
#             for k1,v1 in self.items():
#                 for k2,v2 in other.items():
#                     answer[tuple(tuple(k1[i-1] for i in x) for x in k2)] = v1*v2
#             return EZ_pS_p_element(answer)
    
#     def __repr__(self):
#         '''...'''
#         s = super().__repr__()
#         if self:
#             return s.replace('})',f'}}, p={self.prime})')
#         if not self:
#             return s.replace('()', f'({{}}, p={self.prime})')
    
#     def __str__(self):
#         '''...'''
#         self.reduce_rep()
#         if not self:
#             return '0'
#         else:
#             s = super().__str__()
#             return s.replace(', ',',')
    
#     @staticmethod
#     def all_elements():
#         '''...'''
#         p = Z_p_Module_element.prime
#         return ( Z_pC_p_element( dict(zip(range(p), coeffs)) ) 
#                    for coeffs in product(range(p),repeat=p) )


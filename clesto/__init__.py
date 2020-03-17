from collections import Counter
from itertools import combinations, product, chain, permutations
from math import floor, factorial

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
    """
    Counter with arithmetic improvements to handle (modular) integer values.

    Class constructed to model free module elements over the ring of 
    (modular) integers.

    Attributes
    ----------
    torsion : int or None 
        Chooses the underlying ring R. 
        An int n sets R = Z/nZ whereas None sets R = Z

    """

    torsion = None 

    def __init__(*args, **kwds):
        # print('initializing as Module_element')
        self, *args = args
        super(Module_element, self).__init__(*args, **kwds)
        
        if not all( [type(v) is int for v in self.values()] ):
            raise TypeError('values must be integers')

        self._reduce_rep()

    def __str__(self):
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
        '''The sum of two free module elements.

        >>> Module_element({'a':1, 'b':2}) + Module_element({'a':1})
        Module_element({'a':2, 'b':2})

        '''
        answer = Module_element(self)
        answer.update(other)
        return type(self)( {k:v for k,v in answer.items() if v} )
    
    def __sub__(self, other):
        '''The substraction of two free module elements.

        >>> Module_element({'a':1, 'b':2}) - Module_element({'a':1})
        Module_element({'b':2})

        '''
        answer = Module_element(self)
        answer.subtract(other)
        return type(self)( {k:-v for k,v in answer.items() if v} )

    def __neg__(self):
        '''The additive inverse of a free module element.

        >>> -Module_element({'a':1, 'b':2})
        Module_element({'a':-1, 'b':-22})

        '''
        return type(self)( {k:-v for k,v in self.items() if v} )
    
    def __rmul__(self, c):
        '''The scaling by c of a free module element.

        >>> 3*Module_element({'a':1, 'b':2})
        Module_element({'a':3, 'b':6})

        '''
        answer = Module_element()
        if type(c) is int:
            for key, value in self.items():
                answer[key] = value*c
            return type(self)(answer)
        else:
            raise TypeError(f"can't act by non-int of type {type(c)}")

    def __mod__(self, n):
        '''The reduction mod n of the coefficients of a free module element.

        >>> Module_element({'a':3, 'b':2}) % 2
        Module_element({'a':1})

        '''
        return type(self)( {k:v%n for k,v in self.items() if v % n} )
    
    def __iadd__(self, other):
        '''The in place addition of two free module elements.

        >>> x = Module_element({'a':1, 'b':2})
        >>> x += Module_element({'a':3, 'b':6})
        Module_element({'a':4, 'b':8})

        '''
        self.update(other)
        self._reduce_rep()
        return self
    
    def __isub__(self, other):
        '''The in place addition of two free module elements.

        >>> x = Module_element({'a':1, 'b':2})
        >>> x -= Module_element({'a':3, 'b':6})
        Module_element({'a':-2, 'b':-4})

        '''
        self.subtract(other)
        self._reduce_rep()
        return self
    
    def __imod__(self, n):
        '''The in place reduction mod n of the values of a free module element.

        >>> x = Module_element({'a':1, 'b':2})
        >>> x += Module_element({'a':3, 'b':6})
        Module_element({'a':4, 'b':8})

        '''
        for key, value in self.items():
            self[key] = value % p
            
        return self

    def _reduce_rep(self):
        '''The preferred representative of the free module element.

        It reduces all values mod n if torsion is n and removes 
        key:value pairs with value = 0.
        
        >>> Module_element({'a':1, 'b':2, 'c':0})
        Module_element({'a':1, 'b':2})
         
        '''
        # print('reducing as Module_element')
        if self.torsion:
            for key, value in self.items():
                self[key] = value % self.torsion

        zeros = [k for k,v in self.items() if not v]
        for key in zeros:
            del self[key]

    def set_torsion(self, n):
        '''...'''
        self.torsion = n
        self._reduce_rep()
        return self

#_________________________________79_characters________________________________

class Cyclic_Module_element(Module_element):
    '''Modeling elements in Z[C] or Z/nZ[C_n] where C is the 
    infinite cyclic group and C_n the n-torsion cyclic group 
    generated by abd element denoted by a'''
    
    order = None

    psi_dict = {}

    def __str__(self):
        '''...'''
        if not self:
            return '0'
        else:
            answer = '' 
            for exponent, coefficient in self.items():
                if coefficient != 1 or exponent == 0:
                    answer += f'+{coefficient}q^{exponent}'
                else:
                    answer += f'+q^{exponent}'    
            if answer[0] == '+':
                answer = answer[1:]

            return answer.replace('q^0','').replace('q^1','q')
            
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

    def _reduce_rep(self):
        '''in place mod p reduction of the keys'''
        # print('reducing as Cyclic_Module_element')
        if not all([isinstance(k,int) for k in self.keys()]):
            raise TypeError('keys must be integers')
        if self.order:    
            aux = list(self.items())
            self.clear()
            for k,v in aux:
                self[k%self.order] += v
        
        super()._reduce_rep()

    def set_order(self, r):
        self.order = r
        self._reduce_rep()
        return(self)

    @staticmethod    
    def transposition_element(n):
        '''...'''
        if not(isinstance(n,int) and n>0):
            raise TypeError('order must be a positive integer')
        else:
            T = Cyclic_Module_element({1:1, 0:-1})
            T.set_order(n)
            T.set_torsion(n)
            return T

    @staticmethod
    def norm_element(n):
        '''...'''
        if not(isinstance(n,int) and n>0):
            raise TypeError('order must be a positive integer')
        else:
            N = Cyclic_Module_element({i:1 for i in range(n)})
            N.set_order(n)
            N.set_torsion(n)
            return N

    def psi(self, d):
        '''...'''
        # recursive function to find psi(e_d)
        def _psi_on_generator(d):
            if d == 0:
                return Cyclic_DGModule_element({(0,):1})
            else:
                operators = {
                    0: Cyclic_Module_element.norm_element(self.order),
                    1: Cyclic_Module_element.transposition_element(self.order)}

                op_psi = operators[d%2](_psi_on_generator(d-1))
            
                previous = \
                    Cyclic_DGModule_element({(0,)+k:v for k,v in op_psi.items()})
                previous.set_torsion(self.order)
                previous.set_order(self.order)
                
                return previous

        # using linearity of psi knowing psi(e_d)
        answer = Cyclic_DGModule_element()
        answer.set_torsion(self.order)
        answer.set_order(self.order)
        for k1 in _psi_on_generator(d).keys():
            for k2,v2 in self.items():
                to_add = Cyclic_DGModule_element({tuple(k2+i for i in k1): v2})
                to_add.set_torsion(self.order)
                to_add.set_order(self.order)
                answer += to_add
        return answer


#_________________________________79_characters________________________________

class Symmetric_Module_element(Module_element):
    '''...'''

    def __str__(self):
        '''...'''
        if not self:
            return '0'
        else:
            s = super().__str__()
            return s.replace(', ',',')

    def __mul__(self, other):
        '''...'''
        answer = Symmetric_Module_element()
        for k1,v1 in self.items():
            for k2,v2 in other.items():
                answer[tuple(k1[i-1] for i in k2)] += v1*v2
        return answer
    
    def __call__(self, other):
        '''...'''
        if isinstance(other, Symmetric_Module_element):
            return self*other

        if isinstance(other, Barratt_Eccles_element):
            answer = Counter()
            for k1,v1 in self.items():
                for k2,v2 in other.items():
                    answer[tuple(tuple(k1[i-1] for i in x) for x in k2)] = v1*v2
            return Barratt_Eccles_element(answer)

        if isinstance(other, Surjection_element):
            answer = Counter()
            for k1,v1 in self.items():
                for k2,v2 in other.items():
                    answer[tuple(k1[i-1] for i in k2)] = v1*v2
            return Surjection_element(answer)

    def _reduce_rep(self):
        '''...'''
        # print('reducing as Symmetric_Module_element')
        if any([set(k) != set(range(1,len(k)+1)) for k in self.keys()]):
            raise TypeError('keys must be permutations of (1,2,...,r)') 

        super()._reduce_rep()

    def compose(self, *others):
        '''...'''
        pass

    @classmethod
    def all_elements(r):
        pass

#_________________________________79_characters________________________________

class DGModule_element(Module_element):
    '''...'''
        
    def __str__(self):
        string = super().__str__()
        return string.replace(', ', ',')

    def _reduce_rep(self):
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
                
        super()._reduce_rep()
    
    def boundary(self):
        '''...'''
        bdry = type(self)()
        for spx, coeff in self.items():
            for i in range(len(spx)):
                i_face = tuple(spx[:i]+spx[i+1:])
                i_coeff = coeff*((-1)**i)
                bdry += type(self)({i_face: i_coeff})
        return bdry

#_________________________________79_characters________________________________

class Cyclic_DGModule_element(DGModule_element):
    '''...'''

    order = None

    def __str__(self):
        '''...'''
        s = super().__str__()
        s = s.replace(', ', ',')
        return s.replace('(','q^(')

    def _reduce_rep(self):
        '''reduces mod p the keys and values and deletes keys with 0 value 
        or which are degenerate'''

        # print('reducing as Cyclic_DGModule_element')
        if self.order:    
            aux = list(self.items())
            self.clear()
            for x, v in aux:
                y = tuple(i % self.order for i in x)
                self[y] += v

        super()._reduce_rep()

    def phi(self):
        if not self.order:
            raise('Function phi not define for infinite order elements')

        r = self.order
        answer = Counter()
        for k, v in self.items():
            x = []
            for i in k:
                x.append(tuple(j%r + 1 for j in range(i, r+i)))
            answer[tuple(x)] = v
            
        return Barratt_Eccles_element(answer)

    def set_order(self, r):
        '''...'''
        self.order = r
        self._reduce_rep()
        return self

#_________________________________79_characters________________________________

class Barratt_Eccles_element(DGModule_element):
    '''...'''
    
    def compose(self, *others):
        '''...'''
        pass

    def table_reduction(self):
        '''given a set of basis element in the Barratt_Eccles operad, it returns 
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
     
        return Surjection_element(answer)

#_________________________________79_characters________________________________

class Surjection_element(DGModule_element):
    '''...'''
    
    def compose(self, *others):
        '''...'''
        pass

    def _reduce_rep(self):
        '''...'''

        zeros = [k for k in self.keys() if set(k) != set(range(1,max(k)+1))]
        for k in zeros:
            del self[k]

        super()._reduce_rep()

#_________________________________79_characters________________________________

class Eilenberg_Zilber_element(Module_element):
    '''...'''
    def __str__(self):
        '''...'''
        string = ''
        for multiop, coeff in self.items():
            if coeff != 1:
                string += str(coeff)
            string += '('
            for op in multiop:
                if not op:
                    d = 'id'
                else:
                    d = f'd_{"d_".join(str(i) for i in op)}'

                string += d+')x('
            string = string[:-2]+' + '
        return string[:-3]

    def _reduce_rep(self):
        '''...'''
        aux = list(self.items())
        self.clear()
        for multiop, coeff in aux:
            ordered_multiop = tuple()
            for op in multiop:
                # check input
                if not(isinstance(multiop,tuple) and isinstance(op,tuple) and
                       all([isinstance(i, int) for i in op])):
                    raise TypeError('keys must be tuple of tuple of int')
                # order input
                deg_maps  = Eilenberg_Zilber_element._face_maps_sort(op)
                ordered_multiop += (deg_maps,)
            self[ordered_multiop] += coeff

        super()._reduce_rep()

    @staticmethod
    def _face_maps_sort(face_maps):
        '''puts the face maps in canonical order d < ... < d using the
        simplicial identity d_i d_j = d_j d_{i+1} if i >= j'''

        face_maps = list(face_maps)
        for index in range(1, len(face_maps)):

            currentvalue = face_maps[index]
            position = index

            while (position > 0 and
                   face_maps[position-1] >= currentvalue):

                face_maps[position] = face_maps[position-1]+1
                position = position-1

            face_maps[position] = currentvalue

        return tuple(face_maps)
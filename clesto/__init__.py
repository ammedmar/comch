from collections import Counter
from itertools import combinations, product, chain, permutations
from math import floor, factorial

# TODO: 
# fix EZ_element
# change psi, phi, and table_reduction to to_... methods

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
    default_torsion : int or string 'free' 
        Chooses the underlying ring R. 
        An int n sets R = Z/nZ whereas 'free' sets R = Z

    """

    default_torsion = 'free' 

    def __init__(self, data=None, torsion=None):
        # print('initializing as Module_element')

        # checking input data: dict with int values
        if data:
            if not ( isinstance(data, dict) and
                     all((type(v) is int for v in data.values()))
                ):
                raise TypeError('input must be dict with int values')

        # checking input torsion: positive int or 'free'
        if not torsion is None:
            if not( isinstance(torsion, int) and torsion > 0 
                    or torsion == 'free'
                ):
                raise TypeError("torsion must be a positive int or 'free'") 
        
        # setting torsion
        n = torsion if torsion else type(self).default_torsion
        setattr(self, 'torsion', n)

        # initialize element
        super(Module_element, self).__init__(data)

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
        self.compare_attributes(other)
        answer = type(self)(self).copy_attrs_from(self)
        answer.update(other)
        answer._reduce_rep()
        return answer

    def __sub__(self, other):
        '''The substraction of two free module elements.

        >>> Module_element({'a':1, 'b':2}) - Module_element({'a':1})
        Module_element({'b':2})

        '''
        self.compare_attributes(other)
        answer = type(self)(self).copy_attrs_from(self)
        answer.subtract(other)
        answer._reduce_rep()
        return answer
    
    def __rmul__(self, c):
        '''The scaling by c of a free module element.

        >>> 3*Module_element({'a':1, 'b':2})
        Module_element({'a':3, 'b':6})

        '''
        if not isinstance(c,int):
            raise TypeError(f"can't act by non-int of type {type(c)}")

        scaled = {k:c*v for k,v in self.items()}
        answer = type(self)(scaled).copy_attrs_from(self)
        return answer

    def __neg__(self):
        '''The additive inverse of a free module element.

        >>> -Module_element({'a':1, 'b':2})
        Module_element({'a':-1, 'b':-22})

        '''
        return self.__rmul__(-1)
    
    def __iadd__(self, other):
        '''The in place addition of two free module elements.

        >>> x = Module_element({'a':1, 'b':2})
        >>> x += Module_element({'a':3, 'b':6})
        Module_element({'a':4, 'b':8})

        '''
        self.compare_attributes(other)
        self.update(other)
        self._reduce_rep()
        return self
    
    def __isub__(self, other):
        '''The in place addition of two free module elements.

        >>> x = Module_element({'a':1, 'b':2})
        >>> x -= Module_element({'a':3, 'b':6})
        Module_element({'a':-2, 'b':-4})

        '''
        self.compare_attributes(other)
        self.subtract(other)
        self._reduce_rep()
        return self

    def _reduce_rep(self):
        '''The preferred representative of the free module element.

        It reduces all values mod n if torsion is n and removes 
        key:value pairs with value = 0.
        
        >>> Module_element({'a':1, 'b':2, 'c':0})
        Module_element({'a':1, 'b':2})
         
        '''
        # print('reducing as Module_element')

        # reducing coefficients mod torsion    
        if self.torsion != 'free':
            for key, value in self.items():
                self[key] = value % self.torsion

        # removing key:value pairs with value = 0
        zeros = [k for k,v in self.items() if not v]
        for key in zeros:
            del self[key]

    def set_torsion(self, n):
        '''...'''
        setattr(self, 'torsion', n)
        self._reduce_rep()
        return self

    def compare_attributes(self, other):
        '''...'''
        if self.__dict__ != other.__dict__:
            raise AttributeError('not the same attributes')

    def copy_attributes_to(self, other):
        '''...'''
        for attr, value in self.__dict__.items():
            setattr(other, attr, value)
        other._reduce_rep()
        return(other)

    def copy_attrs_from(self, other):
        '''...'''
        for attr, value in other.__dict__.items():
            setattr(self, attr, value)
        self._reduce_rep()
        return(self)

#_________________________________79_characters________________________________

class Cyclic_Module_element(Module_element):
    '''Modeling elements in Z[C] or Z/nZ[C_n] where C is the 
    infinite cyclic group and C_n the n-torsion cyclic group 
    generated by abd element denoted by a'''
    
    default_order = 'infinite'

    def __init__(self, data=None, torsion=None, order=None):
        # print("initializing as Cyclic_Module_element")

        # checking input data: dict with int keys
        if data:
            if not( isinstance(data, dict) and
                    all((type(k) is int for k in data.keys()))
                ):
                raise TypeError('data type must be dict with int keys')

        # checking input order: positive int or 'infinite'
        if not order is None:
            if not( isinstance(order, int) and order > 0 
                    or order != 'infinite'
                ):
                raise TypeError("order must be a positive int or 'infinite'")

        # setting order
        n = order if order else type(self).default_order
        setattr(self, 'order', n)

        # initializing element
        super(Cyclic_Module_element, self).__init__(data=data, torsion=torsion)

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
        self.compare_attributes(other)
        answer = type(self)().copy_attrs_from(self)
        for k1,v1 in self.items():
            for k2,v2 in other.items():
                answer[k1+k2] += v1*v2
        answer._reduce_rep()
        return answer
        
    def __call__(self, other):
        '''...'''
        if isinstance(other, Cyclic_Module_element):
            return self.__mul__(other)

        if isinstance(other, Cyclic_DGModule_element):
            self.compare_attributes(other)
            answer = type(other)()
            for attr, value in other.__dict__.items():
                setattr(answer,attr,value)
            for k,v1 in self.items():
                for x,v2 in other.items():
                    y = tuple(k+i for i in x)
                    answer[y] += v1*v2
            answer._reduce_rep()
            return answer

    def _reduce_rep(self):
        '''in place mod p reduction of the keys'''
        # print('reducing as Cyclic_Module_element')

        # reducing keys mod order
        if self.order != 'infinite':    
            aux = list(self.items())
            self.clear()
            for k,v in aux:
                self[k%self.order] += v
        
        super()._reduce_rep()

    def set_order(self, r):
        setattr(self, order, r)
        self._reduce_rep()
        return(self)

    @staticmethod    
    def transposition_element(n):
        '''...'''
        if not(isinstance(n,int) and n>0):
            raise TypeError('order must be a positive integer')
        else:
            return Cyclic_Module_element({1:1, 0:-1}, torsion=n, order=n)
            
    @staticmethod
    def norm_element(n):
        '''...'''
        if not(isinstance(n,int) and n>0):
            raise TypeError('order must be a positive integer')
        else:
            return Cyclic_Module_element( {i:1 for i in range(n)}, 
                                           torsion=n, order=n )

    def psi(self, d):
        '''...''' 
        # recursive function to find psi(e_d)
        def _psi_on_generator(d):
            if d == 0:
                return Cyclic_DGModule_element({(0,):1}, torsion=self.torsion,
                                                         order=self.order)
            else:
                operators = {
                    0: Cyclic_Module_element.norm_element(self.order),
                    1: Cyclic_Module_element.transposition_element(self.order)}

                op_psi = operators[d%2](_psi_on_generator(d-1))
            
                return Cyclic_DGModule_element(
                            {(0,)+k:v for k,v in op_psi.items()},
                            torsion=self.torsion, order=self.order)

        # using linearity of psi knowing psi(e_d)
        answer = Cyclic_DGModule_element().copy_attrs_from(self)
        for k1 in _psi_on_generator(d).keys():
            for k2,v2 in self.items():
                to_add = Cyclic_DGModule_element({tuple(k2+i for i in k1): v2},
                                        torsion=self.torsion, order=self.order)
                answer += to_add
        return answer

#_________________________________79_characters________________________________

class Symmetric_Module_element(Module_element):
    '''...'''

    def __init__(self, data=None, torsion=None):
        #print("initializing as Symmetric_Module_element")

        # checking input data: dict with tuple of int keys
        if data:
            if not( isinstance(data, dict) and
                    all(isinstance(perm, tuple) for perm in data.keys()) and
                    all(isinstance(i, int) for i in 
                    chain.from_iterable(data.keys()))
                ):
                raise TypeError('data type must be dict '+
                                'with tuple of int keys')
            
            if any((set(k) != set(range(1,len(k)+1)) for k in data.keys())):
                raise TypeError('keys must be permutations of (1,2,...,r)')

        # setting attribute arity
        arity = None # arity of 0 element
        if data:
            arities = set(max(k) for k in data.keys())
            if len(arities) != 1:
                raise TypeError('keys must have equal arity')  
            else:
                arity = arities.pop()
        setattr(self, 'arity', arity)

        # initialize element
        super(Symmetric_Module_element, self).__init__( data=data,
                                                        torsion=torsion )        

    def __str__(self):
        '''...'''
        if not self:
            return '0'
        else:
            s = super().__str__()
            return s.replace(', ',',')

    def __mul__(self, other):
        '''...'''
        self.compare_attributes(other)
        answer = type(other)().copy_attrs_from(self)
        for k1,v1 in self.items():
            for k2,v2 in other.items():
                answer[tuple(k1[i-1] for i in k2)] += v1*v2
        answer._reduce_rep()
        return answer

    def __call__(self, other):
        '''...'''
        if isinstance(other, Symmetric_Module_element):
            return self*other

        if isinstance(other, Barratt_Eccles_element):
            self.compare_attributes(other)
            answer = type(other)().copy_attrs_from(other)
            for k1,v1 in self.items():
                for k2,v2 in other.items():
                    answer[tuple(tuple(k1[i-1] for i in x) for x in k2)] = v1*v2
            answer._reduce_rep()
            return answer

        if isinstance(other, Surjection_element):
            self.compare_attributes(other)
            answer = type(other)().copy_attrs_from(other)
            for k1,v1 in self.items():
                for k2,v2 in other.items():
                    answer[tuple(k1[i-1] for i in k2)] = v1*v2
            answer._reduce_rep()
            return answer
            
    def compose(self, *others):
        '''...'''
        if self == Symmetric_Module_element():
            # composing with 0 is 0
            return Symmetric_Module_element()

        if len(others) == 2 and isinstance(others[1], int): 
            # partial composition
            answer = Symmetric_Module_element()
            for attr, value in self.__dict__.items():
                setattr(answer, attr, value)

            k = others[1]
            other = others[0]
            for perm1, coeff1 in self.items():
                for perm2, coeff2 in other.items():
                    s = len(perm2)-1
                    to_insert = tuple(i+k-1 for i in perm2)
                    at = perm1.index(k)
                    shifted = tuple(map(lambda i: i+s if i>k else i, perm1))
                    inserted = shifted[:at] + to_insert + shifted[at+1:]
                    new_coeff = coeff1*coeff2
                    to_add = Symmetric_Module_element({inserted:new_coeff})
                    for attr, value in self.__dict__.items():
                        setattr(to_add, attr, value)
                    answer += to_add 
            return answer

        else: 
            # total composition
            if len(others) != self.arity:
                raise TypeError('the number of arguments must equal the arity '
                                +f'of self, but {len(others)} != {self.arity}')

            answer = self
            for idx, other in reversed(list(enumerate(others))):
                answer = answer.compose(other, idx+1)
            return answer
            
    @staticmethod
    def all_elements(r):
        pass

#_________________________________79_characters________________________________

class DGModule_element(Module_element):
    '''...'''
    
    def __init__(self, data=None, torsion=None):
        #print('initializing as DGModule_element')

        # checking input data: dict with tuple keys
        if data:
            if not all((isinstance(x, tuple) for x in data.keys())):
                raise ValueError('data type must be dict with tuple keys')

        # initializing element
        super(DGModule_element, self).__init__(data=data, torsion=torsion)

    def __str__(self):
        string = super().__str__()
        return string.replace(', ', ',')

    def _reduce_rep(self):
        '''deletes degenerate keys and reduces as Module_element'''
        # print('reducing as DGModule_element')

        # removes degenerate simplices
        for simplex, v in self.items():            
            for i in range(len(simplex)-1):
                if simplex[i] == simplex[i+1]:
                    self[simplex] = 0
                
        super()._reduce_rep()
    
    def boundary(self):
        '''...'''
        bdry = type(self)().copy_attrs_from(self)
        for spx, coeff in self.items():
            for i in range(len(spx)):
                i_term = {tuple(spx[:i]+spx[i+1:]) : coeff*((-1)**i)}
                to_add = type(self)(i_term).copy_attrs_from(bdry)
                bdry += to_add
        bdry._reduce_rep()
        return bdry

#_________________________________79_characters________________________________

class Cyclic_DGModule_element(DGModule_element):
    '''...'''

    default_order = 'infinite'

    def __init__(self, data=None, torsion=None, order=None):
        # print("initializing as Cyclic_DGModule_element")

        # checking input data: dict with tuple of int keys
        if data:
            if not( isinstance(data, dict) and
                    all(( isinstance(i, int) for i in 
                    chain.from_iterable(data.keys()) ))
                ):
                raise ValueError('data type must be dict'+
                                 'with tuple of int keys')
        # setting order
        n = order if order else type(self).default_order
        setattr(self, 'order', n)

        #initializing element
        super(Cyclic_DGModule_element, self).__init__(data=data, 
                                                      torsion=torsion)
    def __str__(self):
        '''...'''
        s = super().__str__()
        s = s.replace(', ', ',')
        return s.replace('(','q^(')

    def _reduce_rep(self):
        '''reduces mod p the keys and values and deletes keys with 0 value 
        or which are degenerate'''
        # print('reducing as Cyclic_DGModule_element')
        
        # reducing keys mod order 
        if self.order != 'infinite':
            aux = list(self.items())
            self.clear()
            for x, v in aux:
                y = tuple(i % self.order for i in x)
                self[y] += v

        super()._reduce_rep()

    def phi(self):
        '''from Cyclic_DGModule to Barrat_Eccles'''

        if self.order == 'infinite':
            raise('phi is not define for elements of infinite order')

        r = self.order
        answer = {}
        for k, v in self.items():
            x = []
            for i in k:
                x.append(tuple(j%r + 1 for j in range(i, r+i)))
            answer[tuple(x)] = v

        return Barratt_Eccles_element(answer, torsion=self.torsion)

    def set_order(self, r):
        '''...'''
        setattr(self, order, r)
        self._reduce_rep()
        return self

#_________________________________79_characters________________________________

class Barratt_Eccles_element(DGModule_element):
    '''...'''
    def __init__(self, data=None, torsion=None):
        # print("initializing as Symmetric_Module_element")

        # checking input data: dict with tuple of tuple of int keys
        if data:
            if not( isinstance(data, dict) and
                    all(isinstance(x, tuple) for x in data.keys()) and
                    all(isinstance(perm, tuple) for perm in 
                    chain.from_iterable(data.keys())) and
                    all(isinstance(i, int) for i in 
                    chain.from_iterable(chain.from_iterable(data.keys())))
                ):
                raise TypeError('data type must be dict ' +
                                'with tuple of tuple of int keys')
            
            if any((set(perm) != set(range(1,len(perm)+1)) for perm in
                    chain.from_iterable(data.keys()))):
                raise TypeError('keys must tuples of ' +
                                'permutations of (1,2,...,r)')

        # setting attribute arity 
        arity = None
        if data:
            arities = set()
            for k in data.keys():
                arities_in_k = set()
                for perm in k:
                    if len(perm) != max(perm):
                        raise ValueError(f'the term {perm} is not a permutation')
                    arities_in_k.add(max(perm))
                if len(arities_in_k) != 1:
                    raise ValueError(f'the key {k} is not well formed')
                arities |= arities_in_k # in place union
            if len(arities) != 1:
                raise ValueError('keys must have the same arity')
            arity = arities.pop()
        setattr(self, 'arity', arity)

        # initializing element
        super(Barratt_Eccles_element, self).__init__(data=data, 
                                                     torsion=torsion)
    def domain(self):
        '''returns the integer r such that self is a linear combination 
           of vector of permutations in Sigma_r'''
           
        perm_vect = list(self.keys())[0]
        return len(perm_vect[0])
    
    def _int_to_bin(n):
        '''returns a list representing n in binary'''
        
        if n == 0: return [0]
        
        result = []
        while n > 0:
            result.insert(0, n % 2)
            n = n // 2
        return result
        
    def _paths(p, q):
        '''returns as a list all increasing paths from (0,0) to (p,q)'''
        
        if (p,q) == (0,0): return [((0,0),)]
        
        answer = list()
        if p > 0:
            west = Barratt_Eccles_element._paths(p-1, q)
            for path in west:
                answer.append(path + ((p,q),))
            
        if q > 0:
            south = Barratt_Eccles_element._paths(p, q-1)
            for path in south:
                answer.append(path + ((p,q),))
              
        return answer

    def _sgn_of_path(path):
        '''...'''
        
        segments = range(1,len(path))
        horizontal_segments = []
        vertical_segments = []
        for i in segments:
            vertex1 = path[i-1]
            vertex2 = path[i]
            if vertex2[0] > vertex1[0]:
                horizontal_segments.append(i)
            else:
                vertical_segments.append(i)
     
        ordered_segments = horizontal_segments + vertical_segments
        
        # find the permutation that transforms segments to orderedSegments
        permutation = {}
        for seg in segments:
            for j in range(1,len(ordered_segments)+1):
                if seg == ordered_segments[j-1]:
                    permutation[seg] = j
       
        # compute the signature of the permutation
        sgn = 1
        for i in range(1,len(segments)+1):
            for j in range(i+1,len(segments)+1):
                diff = permutation[j] - permutation[i]
                sgn = diff//abs(diff)
        return sgn
    
    def compose(self, *others):
        '''...'''
        
        if len(others) == 2 and isinstance(others[1], int):
            # partial composition
            answer = Barratt_Eccles_element()
            k = others[1]
            other = others[0]
            for perm_vect1, coeff1 in self.items():
                for perm_vect2, coeff2 in other.items():
                    comp = Barratt_Eccles_element()
                    d = len(perm_vect1)-1
                    e = len(perm_vect2)-1
                    for path in Barratt_Eccles_element._paths(d,e):
                        new_perm_vect = ()
                        for i,j in path:
                            perm1 = Symmetric_Module_element({perm_vect1[i]:1})
                            perm2 = Symmetric_Module_element({perm_vect2[j]:1})
                            partial_comp = perm1.compose(perm2, k)
                            new_perm_vect += (list(partial_comp.keys())[0],)
                        sgn = Barratt_Eccles_element._sgn_of_path(path)
                        comp += Barratt_Eccles_element({new_perm_vect:sgn})
                    answer += coeff1*coeff2*comp
            return answer
        else:
            if not len(others) == self.domain():
                raise TypeError('the number of arguments must be equal to ' + 
                                'the domain of self')
                                
            answer = self
            for idx, other in reversed(list(enumerate(others))):
                answer = answer.compose(other, idx+1)
            return answer

    def table_reduction(self):
        '''given a set of basis element in the Barratt_Eccles operad, it returns 
        the set of surjections in its image via the table reduction morphism'''

        answer = Surjection_element(torsion=self.torsion)
        setattr(answer, 'arity', self.arity)

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
                    answer += Surjection_element({tuple(surjection):value}, 
                                                  torsion=self.torsion)
        answer._reduce_rep()
        return answer

#_________________________________79_characters________________________________

class Surjection_element(DGModule_element):
    '''...'''
    
    def __init__(self, data=None, torsion=None):
        '''...'''

        # checking input data: dict with tuple of int keys
        if data:
            if not( isinstance(data, dict) and
                    all(isinstance(surj, tuple) for surj in data.keys()) and
                    all(isinstance(i, int) for i in
                    chain.from_iterable(data.keys()))
                ):
                raise TypeError('data type must be dict ' +
                                'with tuple of int keys')

        # checking input and setting attribute arity
        arity = None
        if data:
            data_copy = dict(data)
            arities = set()
            for surj in data_copy.keys():
                if set(surj) == set(range(1, max(surj)+1)):
                    arities.add(max(surj))
                else:
                    del data[surj] # degenerate surjection
            if len(arities) != 1:
                raise ValueError('keys must have the same arity')
            arity = arities.pop()
        setattr(self, 'arity', arity)

        # initializing element
        super(Surjection_element, self).__init__(data=data, torsion=torsion)

    def compose(self, *others):
        '''...'''
        pass

    def _reduce_rep(self):
        '''...'''

        zeros = (k for k in self.keys() if 
                 set(k) != set(range(1, self.arity+1)))
        for k in zeros:
            del self[k]

        super()._reduce_rep()

    def to_Eilenber_Zilber_element(self, n):
        '''I forgot how this function works. 
        Bad documentation comes back to bite'''

        def _new_term(term, num_to_append, pos_to_append, k):
            tuple_to_replace = term['seq'][pos_to_append]+(num_to_append,)

            return {'pos': term['pos']+k, 
                    'seq': term['seq'][:pos_to_append] 
                            + (tuple_to_replace,)
                            + term['seq'][pos_to_append+1:]}

        def _get_sign(term):
            '''TBW: given a sequence of "step" of faces of (0,...,n) returns the associate
            sign following the Berger-Fresse convention'''
            return 1
        
        answer = Eilenberg_Zilber_element().copy_attrs_from(self)
        for surj, coeff in self.items():
            after = [{'pos': 0, 
                      'seq': ((),)*(surj[0]-1) + ((0,),) + ((),)*(max(surj)-surj[0])}]
            
            for i in range(n+len(surj)-1):
                before = after[:]
                after = []
                for term in before:
                    if term['pos'] < len(surj)-1:
                        num_to_append = term['seq'][surj[term['pos']]-1][-1]
                        pos_to_append = surj[term['pos']+1]-1

                        empty = bool(term['seq'][pos_to_append])
                        if not empty or num_to_append != term['seq'][pos_to_append][-1]:
                            after.append(_new_term(term, num_to_append, pos_to_append,1))

                    if term['seq'][surj[term['pos']]-1][-1] < n:
                        num_to_append = term['seq'][surj[term['pos']]-1][-1] + 1
                        pos_to_append = surj[term['pos']]-1

                        after.append(_new_term(term, num_to_append, pos_to_append,0))

            for term in after:
                sign = _get_sign(term['seq'])
                multioperator = tuple()
                for seq in term['seq']:
                    multioperator += ( tuple(i for i in range(n+1) if i not in seq),)
                answer += Eilenberg_Zilber_element(
                          {multioperator:sign*coeff}).copy_attrs_from(answer)

        return answer

#_________________________________79_characters________________________________

class Eilenberg_Zilber_element(Module_element):
    '''...'''

    def __init__(self, data=None, torsion=None):
        '''...'''
        
        # checking input and setting attribute arity
        arity = None # arity of the 0 element
        if data:
            # cheking input: dict of tuple of tuple of int
            if not( all((isinstance(multiop, tuple)) for multiop in 
                        data.keys()) and
                    all((isinstance(op, tuple) for op in 
                        chain.from_iterable(data.keys()))) and
                    all((isinstance(i, int) for i in 
                        chain.from_iterable(chain.from_iterable(data.keys()))))
                    ):
                raise TypeError('keys must be tuple of tuple of int')  
        
            # setting arity
            arities = set(len(multiop) for multiop in data.keys())
            if len(arities) != 1:
                raise TypeError('keys must have same arity')
            else:
                arity = arities.pop()        

        setattr(self, 'arity', arity)

        # initialize object
        super(Eilenberg_Zilber_element, self).__init__(data=data, torsion=torsion)

    def __str__(self):
        '''...'''
        if not self:
            return '0'

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
        # print('reducing as Eilenberg_Zilber_element')

        # order face maps in increasing value
        self_data = dict(self)
        self.clear()
        for multiop, coeff in self_data.items():            
            ordered_multiop = tuple()
            for op in multiop:
                ordered_op  = Eilenberg_Zilber_element._face_maps_sort(op)
                ordered_multiop += (ordered_op,)
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

#_________________________________79_characters________________________________

class Power_operation(object):
    '''...'''
    def __init__(self, p, s, d, bockstein=False):
        coeff = (-1)**(floor(d/2)+s)
        if d/2 - floor(d/2):
            coeff *= factorial((p-1)/2)

        b = int(bockstein)
    
        e = Cyclic_Module_element({0:coeff}, torsion=p, order=p)

        self.cyclic = e.psi((2*s-d)*(p-1)-b)
        self.barrat_eccles = self.cyclic.phi()
        self.surjection = self.barrat_eccles.table_reduction()

import re 
from Util import _NumberCollectUtil

LABEL="LABEL"

class Operator(object):
    
    def set_members():
        self.arguments={}














class Vertex(object):

    def _set_members(self):
        self.arguments={LABEL:None}
        self._primes=()

    def __init__ (self,ilabel):
        self._set_members()
        self.set_label(ilabel)
        self._get_prefactor()
        if not self._is_empty():
            self._label=self._strip_primes(self._label)
        self._register()

    ##@param ilabel string with the new name of the vertex (with ticks)
    def set_label(self,ilabel):
        """ sets the label stripped of primes"""
        self._label=ilabel

    def _get_prefactor(self):
        """extracts the prefactor from an operator ans stores it in fac and fac_inv
        """
        helper= _NumberCollectUtil()
        label_list=[self._label]
        self.fac,self.fac_inv=helper.pop_numbers(label_list,1,1)
        if len(label_list)>1:
            raise Exception("number extraction split operator in half:"+self._label)
        print label_list
        self._label=label_list[0] if len(label_list)==1 else ""
        
    def _is_empty(self):
        return self._label==""

    def _strip_primes(self,ilabel):
        """ strips primes from the label and stores their positions"""
        matches=re.finditer("'",ilabel)
        list_=[]
        for match in matches:
            list_.append(match.start())
        self._primes=tuple(list_)
        return re.sub("'","",ilabel)

    def get_primes(self):
        return self._primes
        
    ##@param other other Vertex this is to be compared with
    #@return returns True if the other Vetex has the same Name and number of tics(>=1). 
    def same(self,other):
        """Tests, if the other belongs to the same operator"""
        return ( other._label == self._label ) and\
            len(self._primes) != 0  and\
            (len(other._primes) == len(self._primes))
 
    ##@param other other Vertex this is to be compared with
    #@return returns True if the other Vetex has the same name and position of tics
    def identical(self,other):
        """Test, if the other means an identical Vertex

        Note: if two vertices don't have tics. they can match identical eventhough they don't match in same.
        """ 
        return self._identical2(other)
        
    def _identical(self,other):
        """Test if two vertices are identical.
        

        Two vertices have to be from the same operator to match by this function.
        """
        return  ( other._label == self._label )\
            and len(self._primes) != 0 and\
            other._primes == self._primes

    def _identical2(self,other):
        """Test if two vertices are identical.

        May test true even if same() tests false.
        """
        return  ( other._label == self._label ) and\
            other._primes == self._primes

    ##@param new_label string with the name of the new operator
    def replace(self,new_label):
        self.arguments["LABEL"]=new_label

    def _register(self):
        pass

    def __str__(self):
        return self._label

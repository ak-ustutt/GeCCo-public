import re 


LABEL="LABEL"

class Operator(object):
    
    def set_members():
        self.arguments={}














class Vertex(object):

    def _set_members(self):
        self.arguments={LABEL:None}
        self._tics=()

    def __init__ (self,ilabel):
        self._set_members()
        self.set_label(ilabel)
        self._register()

    ##@param ilabel string with the new name of the vertex (with ticks)
    def set_label(self,ilabel):
        """ sets the label stripped of tics and memorizes the tick muster"""
        self.arguments[LABEL]=self._strip_tics(ilabel)
        
    def _strip_tics(self,ilabel):
        """ strips tics from the label and stores their positions"""
        matches=re.finditer("'",ilabel)
        list_=[]
        for match in matches:
            list_.append(match.start())
        self._tics=tuple(list_)
        return re.sub("'","",ilabel)

    ##@param other other Vertex this is to be compared with
    #@return returns True if the other Vetex has the same Name and number of tics(>=1). 
    def same(self,other):
        """Tests, if the other belongs to the same operator"""
        return ( other.arguments[LABEL] == self.arguments[LABEL] ) and len(self._tics != 0 ) and (len(other._tics) == len(self._tics))
 
    ##@param other other Vertex this is to be compared with
    #@return returns True if the other Vetex has the same Name and number of tics(>=1). 
    def identical(self,other):
        """Test, if the other means an identical Vertex

        Note: if two vertices don't have tics. they can match identical eventhough they don't match in same.
        """ 
        return self._identical1(other) or self._identical2(other)

    def _identical(self,other):
        """Test if two vertices are identical.
        

        Two vertices have to be from the same operator to match by this function.
        """
        return  ( other.arguments[LABEL] == self.arguments[LABEL] ) and len(self._tics != 0 ) and other._tics == self._tics

    def _identical2(self,other):
        """Test if two vertices are identical.

        May test true even if same() tests false.
        """
        return  ( other.arguments[LABEL] == self.arguments[LABEL] ) and other._tics == self._tics

    ##@param new_label string with the name of the new operator
    def replace(self,new_label):
        self.arguments["LABEL"]=new_label

    def _register(self):
        pass

    def __str__(self):
        return self.arguments["LABEL"]

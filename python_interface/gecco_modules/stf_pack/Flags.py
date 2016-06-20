class Flag(list):
    _id=0
    
    collection=None
    @classmethod
    def _make_id(cls):
        cls._id+=1
        return cls._id
        
    def __init__(self):
        list.__init__(self,[self._make_id()])

    def __contains__(self, other):
        if isinstance(other,Flag) and len(other) == 1 :
            return list.__contains__(self, other[0])
        else:
            raise TypeError("in <Flag> requires single Flag on the left, not {typ}".format(typ=type(other)))

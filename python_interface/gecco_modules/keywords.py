from python_interface.gecco_modules import default_keywords as dk
from python_interface import gecco_interface as gi 


ContextError = dk.Context

class KeywordTree(object):
    def __init__(self, input_=gi.GeCCo_Input(), registry=dk.RegistryHandler() ):
        self._input = input_
        self._registry = registry

        
    def is_set(self, context, name):
        return self._input.is_set(context+"." +name)

    def get_value(self, context, name):
        if self.is_set(context, name):
            strvalue = self._input.get(context+"."+name)
            typ = self._registry.get_type(context,name)
            if typ == bool:
                def bool_conv(val):
                    if val == "T":
                        return True
                    elif val == "F":
                        return False
                    else:
                        raise ValueError
                typ = bool_conv
            if isinstance(strvalue, list):
                return [typ(elem) for elem in strvalue ]
            else:
                return typ(elem)
        else:
            return self._registry.get_value(context, name)

keyword_tree = KeywordTree()

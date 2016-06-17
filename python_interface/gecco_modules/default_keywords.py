import re #regular expression support
import xml.etree.ElementTree as ET #xml and DOM functionality
import os # to find the set_keywords file




keyword_file=os.getenv("GECCO_DIR")+"/data/keyword_registry.xml"




_type_dict={8:"str",
            2:"int",
            1:"log",
            4:"rl8",}


_tree=ET.parse(keyword_file)



def _find_element_core(context,root):
    """gets an Element  context
    """
    for item in root.getchildren():
        #iterates through all children of root.
        regexp_context_extr=re.compile("^"+item.get("name")+"(?:\.|$)"+"(?P<rest>.+)?")
        match=regexp_context_extr.match(context)
        #Looks if context is of the form "item_name.<rest>"
        ##todo It may be faster and more straigthforward to test by string slicing 
        #and extract in the end with only one regexp 
        if match and match.group("rest"):
            #there was a match and 
            return _find_parent_core(match.group("rest"),item)
        elif match:
            #there was a match but all context is resolved
            return item
    raise UnknownContextError(context)
    
def _find_element  (context,root):   
    """wrapper for _find_element_core
    """
    try:
        return _find_parent_core(context,root)
    except UnknownContextError as ex :
        print "default_keywords.py"
        print "no context found for: " ,match.group("keyword")
        print match.group("context")+"\n at level:"+ex.context




    

def find_by_context(context,name,tree=_tree):
    full_context=context+"."+name
    root=tree.getroot()

    # finds actually the designated element
    return _find_parent(full_context,root)

def convert_to_bool(string):
    """takes a string of xml bools and converts it to python booleans"""
    if re.match("true",string,flags=re.I):
        return True
    elif re.match("false",string,flags=re.I):
        return False
    else :
        print "convert_to_bool:"
        print "unconvertible string: "+string
        exit
    
def convert_to_float(string):
    """Fortran uses d as exponent separator, this function replaces this by e"""
    try:
        return float(re.sub("[dD]","e",string))
    except:
        print "convert_to_float:"
        print "unconvertible string: "+string
        exit

def get_value_by_context(context,name,tree=_tree):
    elem= find_by_context(context,name,default_keys)
    type_=_type_dict[int(elem.get("type","8"))]
    val=elem.get("value",None)
    if val is None:
        raise Exception("No default value set")

    if  type_=="str":
        return str(val)
        #by xml definition elem.text is saved as string
        #but better make this stable
    elif  type_=="int":
        return int(val)
    elif  type_=="log":
        return convert_to_bool(val)
    elif  type_=="rl8":
        return convert_to_float(val)
    else:
        raise Exception("Unknown type: "+type_)







        

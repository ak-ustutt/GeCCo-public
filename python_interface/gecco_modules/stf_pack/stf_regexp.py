#Regular expressions for the string_to_form package
import re


##@param context ="regexp" "!regexp" (inverted RegExp) or "tuple"
#@return separators either as string (uncompiled regexp) or as tuple (e.g: "[de]","[^de]", ("d","e"))
def _pow_sep(context):
    """ returns the possible separators of the decimal power
    
    Returns the separators either as string (uncompiled regexp) or as tuple
    """
    pow_sep=('d','e')
    if context=="regexp":
        return "["+"".join(pow_sep)+"]"
    elif context =="!regexp":
        return "[^"+"".join(pow_sep)+"]"
    elif context =="tuple":
        return pow_sep
    else:
        raise NotImplementedError("context "+str(context)+"not implemented")

## @param context one of "start", "anywhere" or "bare"
#@return returns a string (uncompiled regular expression)
def _num_regexp_total(context="start"):
    """just a function to store and explain the regexp used to scan for numbers"""
    d=_pow_sep("regexp")
    regexp1="-?[0-9]*"+    "\.[0-9]+"  +"(?:"+d+"-?[0-9]+)?"
    regexp2="-?[0-9]+"+ "(?:\.[0-9]*)?"+"(?:"+d+"-?[0-9]+)?"
    #as described in the documentation of re 
    # "-?" matches zero or one "-"
    # "[0-9]*" matches zero or more digits "[0-9]+" matches one or more digits
        #for internationality this might be changed to "\d"
    # "\." matches a "." backslash for escaping
    # "(?:"  and ")" groups an expression combined with the tailing "?" 
        #this matches zero or one times this group of expressions      
    # more digit matching 
    # more "(?:" ")" and "?"
    # d returns a string that is the regexp for the power separator so "[de]" by default
    # matching anoter group of positive or negative integers
    # the connection by "|" let us have the longer match by either RegExp1 or RegExp2
    # "^" matches the beginning of a string
    # getting it all in Parentheses tells re.split to also return the matched expression
    # So this regexp matches every of the following number formats: 1 , 1. 1.234 , .234
        # all cases can have a preceeding minus or a following exponent(which in turn can 
        #be negative)
    if context == "start" :
        return "(^"+regexp1+"|^"+regexp2+")"
    elif context == "anywhere" :
        return "("+regexp1+"|"+regexp2+")"
    elif context == "bare" :
        return regexp1+"|"+regexp2
    elif context == "bare/start" :
        return "^"+regexp1+"|"+regexp2
    else:
        raise Exception("_num_regexp_total: unkown context")
    

def _OP_regexp_total(context="bare"):
    """a function to store and explain the regexp used to scan formatch OPs"""
    total="(?:" + "[^-+<>\[\],()*]*" + "(?:^+)?" + ")"
    total=".*"
    #An OP (without numerical prefix) is any number of signs that do not include -+<>[],()* 
    # followed by one or zero "^+" and this muster repeated at least zero times. 
    # I wrote into the documentation that ^+ is escaped, so ^+ is escaped.
    # this also matches an empty match (so the numerical regexp also match purely numerical entities)
    # e.g.: -5.63 is matched with an empty operator part
    if context == "bare" :
        return total
######################################################################
#regular expressions


#regexp used to scan for numbers at the start of a string (e.g an operator)
#any match by _NumberRegExp is (apart from _pow_sep) avalid float 
stf_number_regexp=re.compile(_num_regexp_total())

stf_number_extract_regexp=re.compile("(?P<number>"+ _num_regexp_total(context="bare/start") +")"
                                     +"(?P<operator>"+ _OP_regexp_total(context="bare")+")")


#matches anything that begins with a  minus
#second regexp like numbers, but has to be treated as -1 
stf_negative_regexp=re.compile("(^-)")

stf_negative_extract_regexp=re.compile("(^-)(?P<operator>" 
                                  + _OP_regexp_total(context="bare") 
                                  + ")")

#ignore case, but if you ever type bCh, I will hate you
stf_bch_regexp=re.compile("^BCH",flags=re.I)



# matches any division of numbers
stf_div_regexp=re.compile("^(?P<nominator>"
                          + _num_regexp_total(context="bare") +
                          ")/(?P<denominator>"
                          + _num_regexp_total(context="bare") + ")"  )


stf_div_extract_regexp=re.compile("^(?P<nominator>"
                                  + _num_regexp_total(context="bare")
                                  +")/(?P<denominator>"
                                  + _num_regexp_total(context="bare")
                                  +")(?P<operator>" 
                                  + _OP_regexp_total(context="bare") 
                                  +")" )

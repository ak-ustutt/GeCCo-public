from python_interface.gecco_interface import *


ntest=000
def _print_str(string):
    """ Wrapper for PRINT """
    PRINT({'STRING':string})

def comment(string):
    """Easy way to print comments. No formatting yet"""
    _print_str(string)

def heading(string):
    """Easy way to print headings for a target. No formatting yet"""
    _print_str(string)

def notice(arg_name,arg_val="",printlvl=0): 
    """ For notifications of used input values

    will print a notification of the form 
    arg_name  : arg_val (type of arg_val)
    @param arg_name argument name
    @param arg_val argument value 
    """
    _print_str(arg_name+" :  "+str(arg_val)+"("+str(type(arg_val))+")")   

def mark(string,printlvl=0):
    """To set marks for the tests"""
    _print_str("Mark: "+string)

#Debug functions
def debug_MEL(label,only_this=False,nthresh=100,info_only=False):
    """function to debug, prints Me-lists"""
    global ntest
    if ( ntest>=nthresh or only_this):
        PRINT({'STRING':'This is ' + label})
        PRINT_MEL_INFO({'LIST':label})
        if (not info_only):
            PRINT_MEL({'LIST':label})

def debug_FORM(label,only_this=False,nthresh=50,mode='SHORT', output=""):
    global ntest
    interm_dict={'LABEL':label,'MODE':mode}
    if ( output != ""):
        interm_dict['OUTPUT']=output 
    if ( ntest >= nthresh or only_this):
        PRINT({'STRING':'This is the formula ' + label})
        PRINT_FORMULA(interm_dict)



from python_interface.gecco_interface import *


if keywords.is_keyword_set('general.print'):
    ntest = int(keywords.get('general.print'))
else:
    ntest = 0

def _print_str(string):
    """ Wrapper for PRINT """
    PRINT({'STRING':string})

def comment(string):
    """Easy way to print comments. No formatting yet"""
    _print_str("--> " + string)

def heading(string):
    """Easy way to print headings for a target."""
    _print_str("===== " + string + " =====")

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
    if (ntest >= printlvl):
        _print_str("Mark: "+string)

#Debug functions
def debug_MEL(label,only_this=False,nthresh=100,info_only=False):
    """Debugs the ME-lists

    label              The label for the MEL
    @param only_this   If True, prints information for any print level. Default = False
    @param nthresh     Prints the information if print level is higher than nthresh. Default = 100
    @param info_only   Prints only informations about the ME, but not the ME itself. Default = False
    """
    global ntest
    if ( ntest>=nthresh or only_this):
        PRINT({'STRING':'This is ' + label})
        PRINT_MEL_INFO({'LIST':label})
        if (not info_only):
            PRINT_MEL({'LIST':label})

def debug_FORM(label,only_this=False,nthresh=50,mode='SHORT', output=""):
    """Debugs the Formula

    label       The label for the formula
    @param only_this   If True, prints information for any print level. Default = False
    @param nthresh     Prints the information if print level is higher than nthresh. Default = 50
    @param mode        The MODE argument for printing formula. Default = SHORT
    @param output      The OUTPUT argument for printing formula. Default = ""
    """
    global ntest
    interm_dict={'LABEL':label,'MODE':mode}
    if ( output != ""):
        interm_dict['OUTPUT']=output 
    if ( ntest >= nthresh or only_this):
        PRINT({'STRING':'This is the formula ' + label})
        PRINT_FORMULA(interm_dict)


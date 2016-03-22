from gecco_interface import * # This is an extension to the gecco interface. 
                              # it relies on the functions provided there 

import re #support for Regular expressions
import itertools as it # some iterators
#import traceback as tb # to print exceptions

from stf_exceptions import * #custom made exceptions for this package
from stf_regexp import * # The regular expressions for numbers ... 
from stf_string_expansion import InputString,BracketRep  # what I actually need to expand the string 


from ArneUtil import combine_dicts  # a function to combine dictionaries 

##@module A module to convert string representations of formulas to formula objects
#
#Allows to create formula objects and to set them as multiple EXPAND_OP_PRODUCT from a string
#\author Arne Bargholz
#\date 21.05.2015
#\version 1.09 tested, but without updated documentation


#*********************************************************
#General idea: 
#String expansion: 
#    The expansion of a string into a Formula object is handled in two steps
#    In the first step the object hierachy is build from the string. 
#    This is done by taking the string and passing it to a lower tier object, which is stored 
#    in a member list (self.content)
#        The lowest tier objects (_OPRep) contain strings of the operator names.
#    In the second step the extract method demands extraction of all 
#    lower tier objects and then returns a sum of operator products (a list of lists)
#Generation of prefactors
#    The Bracket object (_ExpBracket) 
#     (which at this point contains only a list of strings representing OP product)
#     scans every string with an regular expression(that only fits on leading numbers)
#     Every found number is combined into a single factor that is later given as FAC argument
#Generation of Idx_SV-lists
#     The Bracket object (_ExpBracket) walks over all OPs. 
#     Every OP without a single quote in it is deemed a standalone vertex 
#     and get its own number in the idx_sv list.
#     those with single quotes are compared to each other. (see specs for details, )
## \todo don't forget to change this comment when the  implementation changes

#General considerations: 
#     this module works not on a naive string but on a object termed _InputString()
#     differences (presign of minus) is ignored in the extraction step 
#     (considered part of an OP and invokes most times an addition) and 
#     collected in the number collection as a prefactor of -1

#*****************************************************************
#Naming conventions:
#    Python naming conventions are generally followed. although 
#    I forgot to make several methods of hidden objects hidden
#    The Objects involved in the string unpacking are named "_"+ShortenedOperation+"Rep"
#       Rep for representation
#    The Util classes collect methods for a certain functionality. 
#      they should be considered abstract classes
#    every class has a set_member function that initializes its permanent variables.
#      Why? Because  like my variables initialized
#    There is a distinction between quantenmechanical operators (refferered to as OP) 
#    and the operatrions that connect them (termed operations) 

#What to find where? 
# everything is found in the following order:
# class declarations 
# Exception declarations
# global functions for functionality
# global objects
# class definions for _InputString
# class definitions concerning the disassembling and extraction of the string
# class definitions for the formula object
# class definitions and functions for the user interface
# class definiton of the Test-object


#Possible arguments for EXPAND_OP_PRODUCT:
#I could make a dictionary, mapping the actual keywords to placeholder I define 
# Nah, it's good like this, effectively the list positions are my placeholders

_g_argument_names=[
              'LABEL', #0
              'OP_RES',#1
              'OPERATORS',#2
              'IDX_SV', #3
              'NEW',    #4
              'TITLE', #5
              'FAC',    #6
              'FAC_INV',#7 
              'CONNECT',#8
              'AVOID',  #9
              'LABEL_DESCR',#10
              'FIX_VTX',#11
              'BLK_MIN',#12
              'BLK_MAX',#13
              'N_CONNECT',#14
              'N_AVOID',#15
              'INPROJ',#16
              'N_INPROJ',#17
              'N_DESCR',#18
]




_g_required_args=[LABEL,OP_RES,OPERATORS,IDX_SV]
###############################################
#Classes

##########################################################################################




class _NumberCollectUtil(object):
    """Utility class to derive a factor from operator prefactors"""

#    ##
#    #@param parts output of a re.split function
#    #@param OP_no current OP
#    #@return returns 0 if no (zero) splits occured
#    #         returns 1 if a split occurred
#    def process_parts(self,parts,OP_no):
#        """processes the output of a re.split operation"""#
#
#        if len(parts) not in (1,3):
#            raise BadSplitError("Operator"+str(self.operators[OP_no])+
#                                "was split into "+str(len(parts))+" parts")
#        elif len(parts) == 1:
#            return 0
#        else:
#            self.operators[OP_no]=parts[2]
#            if stf_number_extract_reg_exp.match(parts[1]):
#                self.factor=self.factor*float(re.sub(_pow_sep("reg_exp"),"e",parts[1])) 
#            elif _special_reg_exp.match(parts[1]):
#                self.factor=-self.factor
#            return 1

     
#    ##
#    #@param regexp a compiled regular expression
#    #@param OP_no number of current operator
#    #@return returns the return of self.process parts
#
#    def split(self,regexp,OP_no):
#        """ Splits by a regular expression """
#        #
#        return self.process_parts(regexp.split(self.operators[OP_no]),OP_no)
#        #if there was a match:
#        #[0] is an empty string 
#            #because our strings match the beginning. read documentation of re
#        #[1] contains our match 
#        #[2] contains the rest (at least an empty string)
#        # if there is no match:
#        #[0] contains our string

    ##
    #@param string str to be matched
    #@param nominator
    #@param denominator
    #@return tuple (x,y,z): x new nominator, y new denominator, 
    #z nonnumber part of string.

    def _eval_numbers(self, string, nominator=1, denominator=1):
        """converts the number matched by the regexps to floats

 
        """
        residue=""
        match_div=stf_div_extract_regexp.match(string)
        match_num=stf_number_extract_regexp.match(string)
        match_neg=stf_negative_extract_regexp.match(string)
        if match_div is not None :
            nominator=nominator*float(match_div.group("nominator") )
            denominator=denominator*float(match_div.group("denominator") )
            residue=match_div.group("operator")
        elif match_num is not None :
            nominator=nominator*float(match_num.group("number"))
            residue=match_num.group("operator")
        elif match_neg is not None :
            nominator=-nominator
            residue=match_neg.group("operator")
        else:
            residue=string
        nominator,denominator=self._normalize(nominator,denominator)
        return nominator,denominator,residue


    def _pop_numbers(self):
        """Collects the numbers from operators
        
        Find all numerical prefixes,
        multiply the numbers, eliminate the numbers 
        from the elements and delete those elements that become empty
        
        """
        try:
            OPs=self._OPs
        except KeyError:
            print "pop_numbers: 'OPERATORS' not set"
            raise 
        #identical to self.arguments.get('FAC',1)
        nominator=self.arguments.get(FAC,1)
        #identical to self.arguments.get('FAC_INV',1)
        denominator=self.arguments.get(FAC_INV,1)

        #indexing the list from the end to avoid messing with indices during deletion
        i=-len(OPs) # 0
        while i<0:
            i+=1

            a,b,c=self._eval_numbers(OPs[i])
            nominator*=a
            denominator*=b
            OPs[i]=c

            if OPs[i] == '' or OPs[i] is None :
                del OPs[i]
                # OPs is an alias not a copy for self.arguments[OPERATORS]
                # so its not required to copy it back to arguments[OPERATORS]
        self.arguments[FAC]=nominator
        self.arguments[FAC_INV]=denominator
    

    def _normalize(self,nominator,denominator):
        """function called to normalize nominator and denominator

        Does nothing right now.
        """
        return nominator,denominator


class _IDXUtil(object):
    """Collection of methods to create the idx_sv list"""
    
    def _set_idx_sv(self):
        """function to generate a idx_sv list

        """
        #generator that returns a row of increasing numbers
        counter=it.count(1)
        #Two lists: 
        #the first contains the operators with ticks ("'") 
        #the second one their respective idx_indices.
        special_OPs=[]
        special_OP_numbers=[]
        idx_sv=[]
        try:
            OPs=self._OPs
        except KeyError:
            print "set_idx_sv: 'OPERATORS' not set"
            raise 
        for OP in OPs:
            if OP.find("'")>=0:
                if OP.replace("'","") in special_OPs:
                    index=special_OPs.index(OP.replace("'",""))
                    idx_sv.append(special_OP_numbers[index])
                else:
                    special_OPs.append(OP.replace("'",""))
                    special_OP_numbers.append(counter.next())
                    idx_sv.append(special_OP_numbers[-1])
            else:
                idx_sv.append(counter.next())
        self.arguments[IDX_SV]=idx_sv


class _OPProduct(_IDXUtil,_NumberCollectUtil):
    """Handles a sequence of multiplicatively joined operators"""

    def _set_members(self):
        #self.arguments[OPERATORS] contains the OPs as they appear 
        #in EXPAND_OP_PRODUCT(without ticks).
        #self.OPs contains the original copy
        self.arguments={OPERATORS:[],FAC:1,FAC_INV:1,IDX_SV:[]}
        self._OPs=[]

    def __init__(self,OPs):
        self._set_members()
        self._OPs=OPs
        self._pop_numbers()
        self._set_idx_sv()
        self._process_operators()

    def _strip_ticks_from_OPs(self):
        """returns the OPs without single quotes"""
        return  map(lambda x: re.sub("'","",x),self._OPs)

    def _process_operators(self):
        """Put OPs into the arguments dict 

        Puts strips the single quotes from the OPs and puts them into the arguments
        dictionary.
        """
        self.arguments[OPERATORS]=self._strip_ticks_from_OPs()

    def __str__(self):
        return str(self.arguments[FAC]) + "/" + str(self.arguments[FAC_INV])\
            + "*<"+"*".join(self.arguments[OPERATORS])+">"

    def _combine_dicts(self,bracket_dict={},other_dict={}):
        """Combines the arguments dictionary with another dict

        The foreign dictionary takes precedence, 
        except in for the FAC and FAC_INV keys. their arguments are multiplicated 
        """
        if other_dict=={}:
            other_dict=self.arguments
        nominator=other_dict[FAC]*bracket_dict[FAC]
        denominator=other_dict[FAC_INV]*bracket_dict[FAC_INV]
        new_dict=combine_dicts(bracket_dict,other_dict)
        return combine_dicts({FAC:nominator,FAC_INV:denominator}, new_dict)

    def _set_rule(self,bracket_dict={},settle=True):
        """Invokes EXPAND_OP_PRODUCT"""
        ges_dic=self._combine_dicts(bracket_dict,self.arguments)
        for argument in _g_required_args:
            if argument not in ges_dic:
                raise MissingArgumentError(argument+\
                        "is missing for EXPAND_OP_PRODUCT of "+str(self))
        if settle:
            EXPAND_OP_PRODUCT(ges_dic)
            return [None]
        else:
            return ges_dic

    def _show(self,form_dic={}):
        return self._set_rule(form_dic,settle=False)

#****************************************************************************#
#A Bracket is a collection of operator products

class _Bracket(_NumberCollectUtil):
    """Handles information, that are specific to every bracket"""

    def set_members(self):
        self.arguments={FAC:1,FAC_INV:1}
        self.OP_products=[]
        #self.content normally is a BracketRep
        self.content=""

    def __init__(self, string):
        self.set_members()
        self.process_string(string)
        self.extract()

    def process_string(self, string):
        string.set_start()
        try:
            string.goto("<")
        except ValueError:
            raise UnexpectedEndOfStringError(string.display_error)
        else : 
            nominator,denominator,residue=self._eval_numbers(string.substring)
            if residue not in ("",None,"*"):
                raise UnknownEntity("there is something non numerical"\
                                   "before this bracket"+string.display_error)
            self.arguments[FAC]=nominator
            self.arguments[FAC_INV]=denominator
        #If bracket has "|" as delimiters, ignore everything before "|"
        which,where=string.find_first( ("|",">") )
        if which == 0 :
            try:
                string.goto("|")            
            except ValueError:
                raise UnexpectedEndOfStringError(string.display_error)
        self.content=BracketRep(string)        
        try:
            string.goto(">")
        except ValueError:
            raise UnexpectedEndOfStringError(string.display_error)


    def extract(self):
        self.OP_products=map( _OPProduct, self.content.extract())

    def __str__(self):
        return str(self.arguments[FAC]) + "/" + str(self.arguments[FAC_INV])+"("\
            +"\n+ ".join(map(str, self.OP_products))+")"

    def _combine_dicts(self,form_dict={}):
        """Combines the arguments dictionary with another dict

        The foreign dictionary takes precedence, 
        except in for the FAC and FAC_INV keys. their arguments are multiplicated 
        """
        return combine_dicts(form_dict,self.arguments)


    def set_rule(self,form_dic={},settle=True):
        """Invokes set_rule"""
        ges_dic=self._combine_dicts(form_dic)
        new=ges_dic.get(NEW,False)
        ret_list=[]
        for prod in self.OP_products:
            ret_list.append( prod._set_rule(ges_dic,settle=settle) )
            ges_dic[NEW]=False
        if settle:
            return [None]
        else:
            return ret_list

    def show(self,form_dic={}):
        return self.set_rule(form_dic,settle=False)

##########################################################################################
#the surface class Formula for handling formula on the python surface

class _Formula(object):
    """Base class for all User interfaces"""


    def set_members(self):
        self.content=[]
        self.arguments={NEW:True, FIX_VTX:True}
        self.string=InputString("")
        self.new=True
        self.OP_res=""
        self.label=""

    def _combine_dicts(self,special_dict={}):        
        """ Creates a dictionary of the member variables 
 
        The dictionary is the basis for all EXPORT_OP_PRODUCT dicts of this formula
        
        """
        return combine_dicts(special_dict,self.arguments)

    def __init__(self):
        self.set_members()

    def preprocess_string(self,string):
        """ prepares the input for processing
        
        Removes whitespaces and
        saves the input as _InputString in self.string
        """
        #remove all white spaces
        string= re.sub(" ", "",string)
        self.string=InputString(string)



    def process_string(self,string):
        while True:
            self.content.append(_Bracket(self.string))
            try:
                string.goto("+")
                string.next()
            except ValueError:
                break

    def set_rule(self,special_dic={},settle=True):
        """Invokes set_rule of the brackets"""
        ges_dic=self._combine_dicts(special_dic)
        new=ges_dic.get(NEW,True)
        ret_list=[]
        for bracket in self.content:
            ret_list+=bracket.set_rule(ges_dic,settle=settle)
            ges_dic[NEW]=False
        if settle:
            return None
        else:
            return ret_list

    def show(self,special_dic={}):
        """Returns the dictionaries that arewould be given to EXPAND_OP_PRODUCT"""
        return self.set_rule(special_dic=special_dic,settle=False)

    def str_short(self):
        """Returns the original input as string"""
        return str(self.string)
        
    def __str__(self):
        return self.strings[0]+"\n+ ".join(self.strings[1:])

    def strings(self):
        """Returns the formula as list of strings representing strings 
        

        """
        beginning= [self.label+":"+self.OP_res+"="]
        body=[str(bracket) for bracket in self.content]
        return beginning+body

    def extract(self):
        """ triggert die extraction in the deeper layers"""
        for bracket in self.content:
            bracket.extract()

    def append_(self,other):
        """Appends another _Formula object.
        
        With tailing underline so other classes can implement .append and access this method.
        """
        self.content=self.content+other.content




    

#User interfaces (classes and functions)
#
class Formula(_Formula):
    """Class exposed to user"""
    def extract_label(self):
        """Extracts the formula label from the input string"""
        self.string.set_start()
        try: 
            self.string.goto(":")<0
        except ValueError:
            raise MissingLabelError("Please precede your formula by <my_label>: ")
        self.arguments[LABEL]=self.string.substring
        if len(self.arguments[LABEL])==0:
            raise MissingLabelError("This formula label is empty")
        self.string.next()
        
    def extract_OP_res(self):
        """Extracts the result operators label from the input string"""
        self.string.set_start()
        try:
            self.string.goto("=")
        except:
            raise MissingOP_ResError("Please precede your formula by"\
                                     "'LABEL:OP_RES=' ")
        self.arguments[OP_RES]=self.string.substring
        if len(self.arguments[OP_RES])==0:
            raise MissingOP_ResError("This result Operator is empty")
        self.string.next()
        
    def preprocess_string(self,string):
        """Preprocesses the input
        
        Removes any whitespaces (including tabs)
        Extracts Formula Label and OP_res
        Puts input in self.string (as _InputString)
        """
        #remove all (Unicode) white space characters
        string= re.sub("\s", "",string)
        #set string as _InputString
        self.string=InputString(string)
        self.extract_label()
        self.extract_OP_res()

    def __init__(self,string):
        """The constructor"""
        self.set_members()     
        self.preprocess_string(string)
        self.process_string(self.string)
        self.extract()

    def append(self,other):
        if isinstance(other, _Formula):
            self.append_(other)

        #Are there Problems with unicode strings?
        elif isinstance(other, basestring):
            interm=_Formula()
            interm.preprocess_string(other)
            interm.process_string(interm.string)
            interm.extract()
            self.append_(interm)

        else:
            raise NotImplemented


def make_dicts(string):
    thisformula=Formula1(string)
    thisformula.show()

def make_formula(string):
    thisformula=Formula1(string)
    thisformula.set_rules()








class Test(object):
    def test_OP_product(verbose=True):
        if verbose:
            print("testing _OPProduct")
        try:
            i=_OPProduct(["C0","H0","T1"])
        except( Exception):
            print("There was an error at creation of _OPProduct")
            tb.print_exc()
            raise
        inp={LABEL:"LAG", 
             IDX_SV:[1,2,3], 
             OP_RES:"L", 
             FAC:1, 
             FAC_INV:1 }
        exp_res= {"OP_RES": 'L',
                  'OPERATORS': ['C0', 'H0', 'T1'], 
                  'LABEL': 'LAG', 
                  'FAC_INV': 1, 
                  'FAC': 1, 
                  'IDX_SV': [1, 2, 3]}
        if ( i._show(inp)!=  exp_res ):
            print("there was a problem with _OPProduct._show()"
                  +"\nprobably a problem in the dictionary generation"
                  +"\nresult: " +str(i.show(inp))
                  +"\n expected res:" +str(exp_res))
            raise Exception

        inp={LABEL:"LAG", 
             OP_RES:"L", 
             FAC:1, 
             FAC_INV:1 }
        if i._show(inp) != exp_res :
            print("the idx_sv generation does not yield correct results"
                   +"\nresult:" + str(i.show(inp))
                   +"\n expected res:" + str(exp_res))
            raise Exception
        
        def try_num(op_list,inp,exp_res):
            
            try:
                i=_OPProduct(op_list)
            except( Exception):
                print("There was an error at creation of _OPProduct")
                tb.print_exc()  
                raise

            if  i._show(inp) != exp_res :
                print( "error"
                       +"\nresult:" + str(i.show(inp))
                       +"\n expected res:" + str(exp_res))
                raise Exception
        op_list=["2C0"]
        inp={LABEL:"LAG", 
             OP_RES:"L", 
             FAC:1, 
             FAC_INV:1 }
        exp_res={'OP_RES': 'L', 
                 'OPERATORS': ['C0'], 
                 'LABEL': 'LAG', 
                 'FAC_INV': 1, 
                 'FAC': 2.0, 
                 'IDX_SV': [1]}
        try_num(op_list,inp,exp_res)        

        op_list=["2.4C0"]
        inp={LABEL:"LAG", 
             OP_RES:"L", 
             FAC:1, 
             FAC_INV:1 }
        exp_res={'OP_RES': 'L', 
                 'OPERATORS': ['C0'], 
                 'LABEL': 'LAG', 
                 'FAC_INV': 1, 
                 'FAC': 2.4, 
                 'IDX_SV': [1]}
        try_num(op_list,inp,exp_res)

        op_list=[".4C0"]
        inp={LABEL:"LAG", 
             OP_RES:"L", 
             FAC:1, 
             FAC_INV:1 }
        exp_res={'OP_RES': 'L', 
                 'OPERATORS': ['C0'], 
                 'LABEL': 'LAG', 
                 'FAC_INV': 1, 
                 'FAC': 0.4, 
                 'IDX_SV': [1]}
        try_num(op_list,inp,exp_res)
        
        
        op_list=["-2.4C0"]
        inp={LABEL:"LAG", 
             OP_RES:"L", 
             FAC:1, 
             FAC_INV:1 }
        exp_res={'OP_RES': 'L', 
                 'OPERATORS': ['C0'], 
                 'LABEL': 'LAG', 
                 'FAC_INV': 1, 
                 'FAC': -2.4, 
                 'IDX_SV': [1]}
        try_num(op_list,inp,exp_res)

        op_list=["-.4C0"]
        inp={LABEL:"LAG", 
             OP_RES:"L", 
             FAC:1, 
             FAC_INV:1 }
        exp_res={'OP_RES': 'L', 
                 'OPERATORS': ['C0'], 
                 'LABEL': 'LAG', 
                 'FAC_INV': 1, 
                 'FAC': -0.4, 
                 'IDX_SV': [1]}
        try_num(op_list,inp,exp_res)     

        op_list=["-2C0"]
        inp={LABEL:"LAG", 
             OP_RES:"L", 
             FAC:1, 
             FAC_INV:1 }
        exp_res={'OP_RES': 'L', 
                 'OPERATORS': ['C0'], 
                 'LABEL': 'LAG', 
                 'FAC_INV': 1, 
                 'FAC': -2.0, 
                 'IDX_SV': [1]}
        try_num(op_list,inp,exp_res)

        op_list=["-4e-001C0"]
        inp={LABEL:"LAG", 
             OP_RES:"L", 
             FAC:1, 
             FAC_INV:1 }
        exp_res={'OP_RES': 'L', 
                 'OPERATORS': ['C0'], 
                 'LABEL': 'LAG', 
                 'FAC_INV': 1, 
                 'FAC': -0.4, 
                 'IDX_SV': [1]}
        try_num(op_list,inp,exp_res) 

        op_list=["1/24C0"]
        inp={LABEL:"LAG", 
             OP_RES:"L", 
             FAC:1, 
             FAC_INV:1 }
        exp_res={'OP_RES': 'L', 
                 'OPERATORS': ['C0'], 
                 'LABEL': 'LAG', 
                 'FAC_INV': 24.0, 
                 'FAC': 1, 
                 'IDX_SV': [1]}
        try_num(op_list,inp,exp_res) 


        op_list=["-C0"]
        inp={LABEL:"LAG", 
             OP_RES:"L", 
             FAC:1, 
             FAC_INV:1 }
        exp_res={'OP_RES': 'L', 
                 'OPERATORS': ['C0'], 
                 'LABEL': 'LAG', 
                 'FAC_INV': 1, 
                 'FAC': -1.0, 
                 'IDX_SV': [1]}
        try_num(op_list,inp,exp_res)

        #There are many more combinations to test for numbers.
        # I am omitting many of them.
        #

        
        op_list=["1/24C0","-H","10.543T1"]              
        inp={LABEL:"LAG", 
             OP_RES:"L", 
             FAC:1, 
             FAC_INV:1 }
        exp_res={'OP_RES': 'L', 
                 'OPERATORS': ['C0','H','T1'], 
                 'LABEL': 'LAG', 
                 'FAC_INV': 24.0, 
                 'FAC': -10.543, 
                 'IDX_SV': [1,2,3]}
        try_num(op_list,inp,exp_res)

        if verbose:
            print("number extraction seems to be working")

        #try_num also works for tests of sv_idx
        #

        
        op_list=["GAM'","H","GAM'"]              
        inp={LABEL:"LAG", 
             OP_RES:"L", 
             FAC:1, 
             FAC_INV:1 }
        exp_res={'OP_RES': 'L', 
                 'OPERATORS': ['GAM','H','GAM'], 
                 'LABEL': 'LAG', 
                 'FAC_INV': 1, 
                 'FAC': 1, 
                 'IDX_SV': [1,2,1]}
        try_num(op_list,inp,exp_res)

        
        op_list=["GAM1'","H","GAM2'"]              
        inp={LABEL:"LAG", 
             OP_RES:"L", 
             FAC:1, 
             FAC_INV:1 }
        exp_res={'OP_RES': 'L', 
                 'OPERATORS': ['GAM1','H','GAM2'], 
                 'LABEL': 'LAG', 
                 'FAC_INV': 1, 
                 'FAC': 1, 
                 'IDX_SV': [1,2,3]}
        try_num(op_list,inp,exp_res)



        op_list=["GAM'","H","GA'M"]              
        inp={LABEL:"LAG", 
             OP_RES:"L", 
             FAC:1, 
             FAC_INV:1 }
        exp_res={'OP_RES': 'L', 
                 'OPERATORS': ['GAM','H','GAM'], 
                 'LABEL': 'LAG', 
                 'FAC_INV': 1, 
                 'FAC': 1, 
                 'IDX_SV': [1,2,1]}
        try_num(op_list,inp,exp_res)


        if verbose:
            print("creation of index list seems to be working")

        #testing __str__
        #
        #reusing last _OPProduct
        i=_OPProduct(op_list)
        try:
            assert str(i) == '1/1*<GAM*H*GAM>'
        except Exception:
            print("Warning _OPProduct.__str__ seems not to be working"
                  +" according to specifications"
                  +"\nresult:"+str(i)) 
        else:
            if verbose:
                print("all test of _OPProduct were successful")


    def test_Bracket(verbose=True):
        def factory(string):
            return _Bracket(InputString(string))

        def try_(string,inp,exp_res):
            try:
                i=factory(string)
            except( Exception):
                print("There was an error at creation of _OPProduct from input"+str(string))
                tb.print_exc()  
                raise

            if  i.show(inp) != exp_res :
                print( "error"
                       +"\nresult:" + str(i.show(inp))
                       +"\n expected res:" + str(exp_res))
                raise Exception

        string="<H>"
        inp={"LABEL":"L", 
             "OP_RES":"L", 
             "NEW":True}

        exp_res=[{'OP_RES': 'L', 'OPERATORS': ['H'], 'LABEL': 'L', 'FAC_INV': 1, 'NEW': True, 'FAC': 1, 'IDX_SV': [1]}]
        try_(string,inp,exp_res)


        string="2.4*<H>"
        inp={"LABEL":"L", 
             "OP_RES":"L", 
             "NEW":True}

        exp_res=[    {'OP_RES': 'L', 
                      'OPERATORS': ['H'], 
                      'LABEL': 'L', 
                      'FAC_INV': 1, 
                      'NEW': True, 
                      'FAC': 2.4, 
                      'IDX_SV': [1]
                  }
        ]
        try_(string,inp,exp_res)

        string="2/7<H>"
        inp={"LABEL":"L", 
             "OP_RES":"L", 
             "NEW":True}

        exp_res=[    {'OP_RES': 'L', 
                      'OPERATORS': ['H'], 
                      'LABEL': 'L', 
                      'FAC_INV': 7.0, 
                      'NEW': True, 
                      'FAC': 2.0, 
                      'IDX_SV': [1]
                  }
        ]
        try_(string,inp,exp_res)


        string="2.5<H>"
        inp={"LABEL":"L", 
             "OP_RES":"L", 
             "NEW":True}

        exp_res=[    {'OP_RES': 'L', 
                      'OPERATORS': ['H'], 
                      'LABEL': 'L', 
                      'FAC_INV': 1, 
                      'NEW': True, 
                      'FAC': 2.5, 
                      'IDX_SV': [1]
                  }
        ]
        try_(string,inp,exp_res)

        
        string="<C0*D|H|C0+-*[[[>"
        inp={"LABEL":"L", 
             "OP_RES":"L", 
             "NEW":True}

        exp_res=[    {'OP_RES': 'L', 
                      'OPERATORS': ['H'], 
                      'LABEL': 'L', 
                      'FAC_INV': 1, 
                      'NEW': True, 
                      'FAC': 1, 
                      'IDX_SV': [1]
                  }
        ]
        try_(string,inp,exp_res)
        
        i=factory(string)
        try:
            assert str(i) == '1/1(1/1*<H>)'
        except Exception:
            print("Warning _OPProduct.__str__ seems not to be working")
        else:
            if verbose:
                print("all test of _Bracket were successful")


    def test_formula(verbose=True):
        inp={"LABEL":"F_LAG","OP_RES":"L"}
        exp_res=[{'OP_RES': 'L', 'FIX_VTX': True, 'OPERATORS': ['H'], 'LABEL': 'F_LAG', 'FAC_INV': 1, 'FAC': 1, 'NEW': True, 'IDX_SV': [1]},
{'OP_RES': 'L', 'FIX_VTX': True, 'OPERATORS': ['H', 'T1'], 'LABEL': 'F_LAG', 'FAC_INV': 1, 'FAC': 0.5, 'NEW': False, 'IDX_SV': [1, 2]},
{'OP_RES': 'L', 'FIX_VTX': True, 'OPERATORS': ['T1', 'H'], 'LABEL': 'F_LAG', 'FAC_INV': 1, 'FAC': -0.5, 'NEW': False, 'IDX_SV': [1, 2]},
{'OP_RES': 'L', 'FIX_VTX': True, 'OPERATORS': ['H', 'T1', 'T1'], 'LABEL': 'F_LAG', 'FAC_INV': 6.0, 'FAC': 1.0, 'NEW': False, 'IDX_SV': [1, 2, 3]},
{'OP_RES': 'L', 'FIX_VTX': True, 'OPERATORS': ['T1', 'H', 'T1'], 'LABEL': 'F_LAG', 'FAC_INV': 6.0, 'FAC': -1.0, 'NEW': False, 'IDX_SV': [1, 2, 3]},
{'OP_RES': 'L', 'FIX_VTX': True, 'OPERATORS': ['T1', 'H', 'T1'], 'LABEL': 'F_LAG', 'FAC_INV': 6.0, 'FAC': -1.0, 'NEW': False, 'IDX_SV': [1, 2, 3]},
{'OP_RES': 'L', 'FIX_VTX': True, 'OPERATORS': ['T1', 'T1', 'H'], 'LABEL': 'F_LAG', 'FAC_INV': 6.0, 'FAC': 1.0, 'NEW': False, 'IDX_SV': [1, 2, 3]}]


        i=_Formula()
        i.preprocess_string("<H>+0.5<[H,T1]>+1/6*<[[H,T1],T1]>")
        i.process_string(i.string)
        assert i.show(inp) == exp_res

        print("all tests of _Formula were successful")       

        #currently only one test for _Formula. 
    def test_Formula(verbose=True):
        print("Formula is currently not tested")
    

    def test_all(verbose=True):
        try:
            LABEL
        except NameError: 
            for i in _g :
                exec(i +"='"+i+"'" )
	test_OP_product(verbose)
        test_Bracket(verbose)
        test_formula(verbose)
        test_Formula(verbose)


"""  A module to convert string representations of formulas to formula objects

Allows to create formula objects and to set them as multiple EXPAND_OP_PRODUCT from a string


\author Arne Bargholz
\date 17.06.2015
\version 1.19 tested,
"""
from gecco_interface import * # This is an extension to the gecco interface. 
                              # it relies on the functions provided there 
import copy # to produce deepcopys of objects
import re #support for Regular expressions
import itertools as it # some iterators

from stf_exceptions import * #custom made exceptions for this package
from stf_regexp import * # The regular expressions for numbers ... 
from stf_string_expansion import InputString,BracketRep  # what I actually need to expand the string 


from Util import combine_dicts,_IDXUtil,_NumberCollectUtil,remove_whites




#*********************************************************
#General idea: 
#Generation of prefactors
#     The generation of prefactors is now done in operators.py
#Generation of Idx_SV-lists
#     The Bracket object (_ExpBracket) walks over all OPs. 
#     Every OP without a single quote in it is deemed a standalone vertex 
#     and get its own number in the idx_sv list.
#     those with single quotes are compared to each other. (see specs for details, )
#Append functionality.
#     Either strings or other formula objects can be appended. If another formula is 
#     appended, only the body part is appended, LABEL and OP_RES stay active
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
#    every class has a set_member function that initializes its permanent variables.
#      Why? Because  like my variables initialized
#    There is a distinction between quantenmechanical operators(vertices actually)
#              (refferered to as OP) 
#    and the operations that connect them (termed operations) 

#What to find where? 
# everything is found in the following order:
# environment preparation for unit tests
# global functions  and objects
# class definitions(each preceeded by base classes and Util classes):
#    _OPProduct
#    _Bracket
#    Formula and GenForm
# functions of the user interface
# class definiton of the Test-object





#for unit testing ensures that all functionality (besides set_rule) can be tested even when gecco_interface is not imported
_g_argument_names=[
              'LABEL', 
              'OP_RES',
              'OPERATORS',
              'IDX_SV', 
              'NEW',    
              'TITLE', 
              'FAC',    
              'FAC_INV',
              'CONNECT',
              'AVOID', 
              'LABEL_DESCR',
              'FIX_VTX',
              'BLK_MIN',
              'BLK_MAX',
              'N_CONNECT',
              'N_AVOID',
              'INPROJ',
              'N_INPROJ',
              'N_DESCR',
]

_g_required_functions=["EXPAND_OP_PRODUCT"]


class _Logger(object):
    def __init__(self):
        self.events=[]

    def log(self, function_name, argument_dict):
        """logs the function name and the arguments"""
        self.events.append((function_name, argument_dict))
        
    def clear(self):
        self.events=[]

#used to detect if gecco_interface is loaded otherwise testmode is assumed
try:
    LABEL
except NameError:
    for i in _g_argument_names:
        exec(i +"='"+i+"'" )
    def quit_error(msg):
        print msg

    logger=_Logger()
    for name in _g_required_functions:
        exec("def {name}(dictionary):\n"\
             "    logger.log('{name}', dictionary)"\
             .format(name=name))

try:
    keywords #keywords is defined in gecco_interface
except NameError:
    #gecco_interface not loaded make mock up
    for i in _g_argument_names:
        exec(i +"='"+i+"'" )
    def quit_error(msg):
        print msg

    logger=_Logger()
    for name in _g_required_functions:
        exec("def {name}(dictionary):\n"\
             "    logger.log('{name}', dictionary)".\
             format(name=name))
else:
    if( keywords.is_keyword_set("method.unit_test") ):
        logger=_Logger()
        #gecco_interface loaded but test mode
        for name in _g_required_functions:
             exec("def {name}(dictionary):\n"\
                 "    logger.log('{name}', dictionary)".\
                 format(name=name))

        
_g_required_args=[LABEL,OP_RES,OPERATORS,IDX_SV]





#----------------------------------------------------------------------------------------
#Classes
#----------------------------------------------------------------------------------------






#----------------------------------------------------------------------------------------
#
#----------------------------------------------------------------------------------------

    




class _OPProduct(object):
    """Handles a sequence of multiplicatively joined Vertices"""

    def _set_members(self):
        #self.arguments[OPERATORS] contains the OPs as they appear 
        #in EXPAND_OP_PRODUCT(without ticks).
        #self.OPs contains the original copy
        self.arguments={OPERATORS:[],FAC:1,FAC_INV:1,IDX_SV:[]}
        self._OPs=[]

#----------------------------------------------------------------------
#   public methods and hooks 
    ##@param OPs a list of Vertex objects 
    def __init__(self,OPs):
        self._set_members()
        self.arguments[FAC],\
            self.arguments[FAC_INV]=self._extract_prefactors(OPs,(1,1))
        OPs=self._eliminate_empty_OPs(OPs)
        self.arguments[IDX_SV]=self._generate_idx_sv(OPs)
        self._OPs=OPs

    ##@return a string representation of the object
    def __str__(self):
        ret=""
        if self.arguments[FAC_INV]!=1:
            ret=str(self.arguments[FAC])+"/"+str(self.arguments[FAC_INV])+"*"
        elif self.arguments[FAC]!=1:
            ret=str(self.arguments[FAC])+"*"
        return ret + "*".join(
            self._OPs_to_string(
                self._OPs)
        )

    ##@param other an _OP_PRODUCT
    def compare(self,other):
        """compares if self is a product of the same vertices as other"""
        #assert isinstance(other,type(self))
        self._compare_operators(self._OPs,other._OPs)


    #@raises MissingArgumentError if an argument is missing
    #@param bracket_dic dic given by _Bracket
    #@param settle parameter that differentiates show and set modus
    def set_rule(self,bracket_dict={},settle=True):
        """Builds argument dict for and invoces EXPAND_OP_PRODUCT
        """
        arg_dict=copy.copy(self.arguments)
        arg_dict[OPERATORS]=self._OPs_to_string(self._OPs)
        ges_dic=self._combine_dicts(arg_dict,bracket_dict)
        self._check_required_args(_g_required_args,ges_dic)
        if settle:
            EXPAND_OP_PRODUCT(ges_dic)
            return [None]
        else:
            return ges_dic

    ##@param args list of arguments, to be checked
    #@param dic dic to be checked
    #@raises MissingArgumentError if an argument is missing
    def _check_required_args(self, args, dic):
        for argument in args:
            if argument not in dic:
                raise MissingArgumentError(argument+\
                        "is missing for EXPAND_OP_PRODUCT of "+str(self))

    ##@param other a tuple of the form (fac,fac_inv) with the preset factors
    def add_prefac(self,other):
        """ adds anothers fac and fac_inv to selfs"""
        facs=other.arguments[FAC],other.arguments[FAC_INV]
        facs_own=self.arguments[FAC],self.arguments[FAC_INV]
        self.arguments[FAC],self.arguments[FAC_INV]=self._add_ratios(facs_own,facs)

#-------------------------------------------------------------------------
# data transformation methods, independent of self
        
    ##@param facs1 a tuple of the form (fac,fac_inv)
    #@param facs2 a tuple of the form (fac,fac_inv) 
    def _add_rations(self, facs1, facs2):
        """adds two ratios """
        fac1,fac1_inv=facs1
        fac2,fac2_inv=facs2
        ##\TODO normalize result of addition
        return fac1*fac2_inv+fac2*fac1_inv, fac1_inv*fac2_inv
        
    ##@param[in] OPs a list of OPs
    #@return a list of integers denoting the supervertices a vertex at a given index belongs to
    def _generate_idx_sv(self,OPs):
        helper=_IDXUtil()
        return helper.generate_idx_sv(OPs)


    ##@param[in] OPs a list of OPs
    #@return a list of OPs none of which is empty
    def _eliminate_empty_OPs(self, OPs):
        """eliminates the empty OPs"""
        return [OP for OP in OPs if not OP._is_empty()]

        
    ##@param OPs list of OPs
    #@param facs a tuple of the form (fac,fac_inv) with the preset factors
    #@return  a tuple of the form (fac,fac_inv) with the combined factors of all OPs
    def _extract_prefactors(self,OPs,facs):
        """extracts the the prefactor and the inverse prefactor from the operator sequence"""
        return reduce(lambda x,y : (x[0]*y.fac,x[1]*y.fac_inv),
                      OPs,
                      facs )
        
    ##@param OPs list of OPs
    ##@return list of OP-names
    def _OPs_to_string(self, OPs):
        """converts all OPs to their Name(without primes)"""
        #removal of primes is done in Vertex
        return [str(x) for x in OPs]

    ##@param OPs1 list of OPs
    #@param OPs2 list of OPs
    #@return True if all OPs in list one test as identical to their
    # respective counterparts in list two
    def _compare_operators(self,OPs1,OPs2):
        """tests if the OPs in two lists are identical"""
        if len(OPs1)!=len(OPs2):
            return False
        else:
            return reduce(lambda x,y : x and y,
                          [OPs1[i].identical(OPs2[i]) for i  in xrange(len(OPs1))])


    ##@param bracket_dict,other_dict two dictionaries containing entries for FAC and FAC_INV
    #
    def _combine_dicts(self,other_dict,bracket_dict={}):
        """Combines the arguments dictionary with another dict

        The foreign dictionary takes precedence, 
        except in for the FAC and FAC_INV keys. their arguments are multiplicated 
        """
        nominator=other_dict[FAC]*bracket_dict[FAC]
        denominator=other_dict[FAC_INV]*bracket_dict[FAC_INV]
        new_dict=combine_dicts(bracket_dict,other_dict)
        return combine_dicts({FAC:nominator,FAC_INV:denominator}, new_dict)



#****************************************************************************#
#A Bracket is a collection of operator products

#declaring them

        ##\TODO move all possible code from _Bracket to _AbstractBracket
class _AbstractBracket(_NumberCollectUtil):
    def _set_members(self):
        self.arguments={FAC:1,FAC_INV:1}
        self.OP_products=[]
        #self.content normally is a BracketRep
        self.content=""
        
class _AltBracket(_AbstractBracket):
    def __init__(self):
        self._set_members()

class _Bracket(_AbstractBracket):
    """Handles information, that are specific to every bracket"""


    ##@param string an _InputString object currently pointing to a bracket opening
    def __init__(self, string):
        self._set_members()
        self.process_string(string)
        self.extract()
        self.sum_OP_Products()

    def sum_OP_Products(self,other=_AltBracket()):
        OPProd_list=self.OP_products+other.OP_products
        self.OP_products=self._sum_prod_list(OPProd_list)
                
    #Functions of the process - extract mechanic
    # see stf_string_expansion for details
    def process_string(self, string):
        string,\
            self.arguments[FAC],\
            self.arguments[FAC_INV]=self._process_prelude(string)
        string=self._process_bra(string)
        string, self.content=self._process_operators(string)
        string=self._process_ket(string)

    def extract(self):
        self.OP_products= [_OPProduct(x) for x in  self.content.extract()]

    ##@param string an _InputString pointing to the start of the numerical
    #prefactor or start of the bracket
    def _process_prelude(self,string):
        """processes everything that is before the actual start of a string"""
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
        return string,nominator,denominator

    ##@param string an _InputString pointing to the start of the numerical        
    def _process_bra(self,string):
        #If bracket has "|" as delimiters, ignore everything before "|"
        which,where=string.find_first( ("|",">") )
        if which == 0 :
            try:
                string.goto("|")            
            except ValueError:
                raise UnexpectedEndOfStringError(string.display_error)
        return string

    def _process_operators(self, string):
        return string,BracketRep(string)

    def _process_ket(self, string):
        try:
            string.goto(">")
        except ValueError:
            raise UnexpectedEndOfStringError(string.display_error)
        return string


    #other functions
    def set_rule(self,form_dic={},settle=True):
        """Invokes set_rule of contained objects

        Combines own arguments with formula lvl arguments (from form_dic)
        and invokes EXPAND_OP_PRODUCT
        """
        ges_dic=self._combine_dicts(form_dic)
        new=ges_dic.get(NEW,False)
        ret_list=[]
        for prod in self.OP_products:
            ret_list.append( prod.set_rule(ges_dic,settle=settle) )
            ges_dic[NEW]=False
        if settle:
            return [None]
        else:
            return ret_list

    def show(self,form_dic={}):
        return self.set_rule(form_dic,settle=False)

    def __str__(self):        
        ret=""
        if self.arguments[FAC_INV]!=1:
            ret=str(self.arguments[FAC])+"/"+str(self.arguments[FAC_INV])+"*"
        elif self.arguments[FAC]!=1:
            ret=str(self.arguments[FAC])+"*"
        return ret + "<\n" + "\n+".join(map(str, self.OP_products))+"\n>"

    ##@param OPProd_list List of OPProducts
    def _sum_prod_list(self, OPProd_list):
        """sums the OPProducts of a list
        """
        ran =xrange(-len(OPProd_list),0,1) 
        for i in ran:
            for j in xrange(i,0,1):
                if OPProd_list[i].compare(
                        OPProd_list[j]):
                    OPProd_list[j].add_prefac(OPProd_list[i])
                    del OPProd_list[i]
                    break
        return OPProd_list

    def _combine_dicts(self,form_dict={}):
        """Combines the arguments dictionary with another dict

        The foreign dictionary takes precedence, 
        except in for the FAC and FAC_INV keys. their arguments are multiplicated 
        """
        return combine_dicts(form_dict,self.arguments)



#************************************************************************
#the surface class Formula for handling formula on the user surface

class _FormulaStringRepUtil(object):
    """Collection of string representation methods 
    
    Collection of methods that represent the Formula as a string or list of strings
    """
    def show(self,special_dic={}):
        """Returns the dictionaries that would be given to EXPAND_OP_PRODUCT"""
        return self.set_rule(special_dic=special_dic,settle=False)

    def str_short(self):
        """Returns the original input as string"""
        return str(self._string)
        
    def __str__(self):
        return self._to_string()

    def _to_string(self):
        return self._strings()[0]+"\n+ ".join(self._strings()[1:])

    def _strings(self):
        """Returns the formula as list of strings representing strings """
        head= [self.arguments[LABEL]+":"+self.arguments[OP_RES]+"="]
        body=[str(bracket) for bracket in self._content]
        return head+body
    

class _Formula( _FormulaStringRepUtil):
    """Base class for all User interfaces"""

    def set_members(self):
        """initialze members
        """
        self._content=[]
        self.arguments={LABEL:None,OP_RES:None,NEW:True, FIX_VTX:True}
        self._string=InputString("")


    def preprocess_string(self,string):
        """ prepares the input for processing
        
        Removes whitespaces and converts string to InputString
        @param string string to be preprocessed
        @return preprocessed string 
        """
        #remove all white spaces
        string=remove_whites(string)
        return InputString(string)


    #Functions of the process - extract mechanic
    # see stf_string_expansion for details
    def process_string(self,string):
        """processes the string to a hierachy of operation objects

        @param string string to process (should be preprocessed)
        @return None 
        """
        while True:
            self._content.append(_Bracket(string))
            try:
                string.goto("+")
                string.next()
            except ValueError:
                break

    def extract(self):
        """triggers the extraction of the deeper layers"""
        for bracket in self._content:
            bracket.extract()

    #
    #Methods to invoke EXPAND_OP_PRODUCT
    ##@param special_dic dictionary of touples that may overwrite the self.arguments of this instance (default: empty)
    #@param settle for internal use only
    #@param cleanup 
    #@return None (the dictionaries fed to EXPORT_OP_PRODUCT if settle is False)
    def set_rule(self,special_dic={},settle=True,flags=None):
        """Invokes set_rule of the brackets
        
        """
        ges_dic=self._combine_dicts(special_dic)
        ret_list=[]
        for bracket in self._content:
            ret_list+=bracket.set_rule(ges_dic,settle=settle)
            ges_dic[NEW]=False
        if settle:
            return None
        else:
            return ret_list

    def _combine_dicts(self,special_dict={}):        
        """ Creates a dictionary of the member variables 
 
        @return dictionary that is the basis for all EXPORT_OP_PRODUCT dicts of this formula
        """
        return combine_dicts(special_dict,self.arguments)

    #Append functionality
    def _append(self,other):
        """Appends another formula body either from _Formula or from a string.

        @parameter other either a formula (only the body is appended) 
        or a  string representing a formula body """
        if isinstance(other, _Formula):
            self._content+=other._content
            self._string.extend(str(other._string))
        elif isinstance(other, basestring):
            interm=_Formula()
            interm.set_members()
            #relying on the parent class, not GenFormula, 
            #  to make Formula independent from other interface classes 
            string=interm.preprocess_string(other)
            interm.process_string(string)
            interm.extract()
            self._append(interm)
        elif isinstance(other,list):
            for member in other:
                self._append(member)
        else:
            raise NotImplementedError

    def append(self,other):
        """Appending a string or a formula

        @parameter other either a formula  
        or a  string representing a formula body """
        try:
            self._append(other)
        except ExpFormError as ex:
            quit_error(ex.msg)


    #the "remove" functionality
    ##\todo implemente remove functionality

    #operator overload
    def __add__(self,other):
        """Overload the addition to append formula or string to a deepcopy of self

        @param other see this_class.append(other)
        @return the extended formula
        """
        interm=copy.deepcopy(self)
        interm.append(other)
        return interm

    def __iadd__(self,other):
        """Overload increment '+=' with appending
        
        @param other see this_class.append(other)
        @return self with 
        """
        self.append(other)
        return self




# interfaces (classes and functions)
#
class GenForm(_Formula):
    """ Interface for generic formulas"""
    ##@param label LABEL of the formula (can be None)
    # @param OP_res OP_RES of this formula (can be None)
    # @param body body of the formula
    def __init__(self,label=None,OP_res=None,body=None):
        """The constructor

        """
        self.set_members()  
        self._set_values(label,OP_res,body)

    def _set_values(self,label,OP_res,body):
        """Deal with the arguments of __init__"""
        if label is not None:
            self.arguments[LABEL]=label
        if OP_res is not None:
            self.arguments[OP_res]=OP_res
        if body is not None:
            self._string=self.preprocess_string(body)
            self.process_string(self._string)

        

    
class Formula(_Formula):
    """Interface for complete formulas """
    ##  @param string string of the form 'LABEL:OP_RES=body'
    def __init__(self,string):
        """The constructor
        
        """
        self.set_members()     
        self._string=self._preprocess_string(string)
        self.process_string(self._string)
        self.extract()
        
    def _extract_label(self,string):
        """Extracts the formula label from the input string"""
        string.set_start()
        try: 
            string.goto(":")<0
        except ValueError:
            raise MissingLabelError("Please precede your formula by LABEL: ")
        self.arguments[LABEL]=string.substring
        if len(self.arguments[LABEL])==0:
            raise MissingLabelError("This formula label is empty")
        string.next()
        return string
        
    def _extract_OP_res(self,string):
        """Extracts the result operator label from the input string"""
        string.set_start()
        try:
            string.goto("=")
        except:
            raise MissingOP_ResError("Please precede your formula by"\
                                     "'LABEL:OP_RES=' ")
        self.arguments[OP_RES]=string.substring
        if len(self.arguments[OP_RES])==0:
            raise MissingOP_ResError("This result Operator is empty")
        string.next()
        return string

        ##@param string string to be preprocessed
        #@return preprocessed string
    def _preprocess_string(self,string):
        """Preprocesses the input
        
        Removes any whitespaces (including tabs)
        Extracts Formula Label and OP_res
        converts string to InputString
        @param string string to be preprocessed
        @return preprocessed string
        """
        string=remove_whites(string) 
        return self._extract_OP_res(self._extract_label(InputString(string)))
        
def make_formula(string):
    Formula(string).set_rule()

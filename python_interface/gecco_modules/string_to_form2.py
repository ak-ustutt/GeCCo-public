from gecco_interface import * # This is an extension to the gecco interface. 
                              # it relies on the functions provided there 
import copy # to produce deepcopys of objects
import re #support for Regular expressions
import itertools as it # some iterators
from traceback import format_exc
#mport traceback as tb # to print exceptions

from stf_exceptions import * #custom made exceptions for this package
from stf_regexp import * # The regular expressions for numbers ... 
from stf_string_expansion import InputString,BracketRep  # what I actually need to expand the string 


from ArneUtil import combine_dicts  # a function to combine dictionaries 

##@module A module to convert string representations of formulas to formula objects
#
#Allows to create formula objects and to set them as multiple EXPAND_OP_PRODUCT from a string
#\author Arne Bargholz
#\date 17.06.2015
#\version 1.19 tested, but without updated documentation


#*********************************************************
#General idea: 
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
#    The Util classes collect methods for a certain functionality. 
#      they should be considered abstract classes
#    every class has a set_member function that initializes its permanent variables.
#      Why? Because  like my variables initialized
#    There is a distinction between quantenmechanical operators (refferered to as OP) 
#    and the operatrions that connect them (termed operations) 

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


#used to detect if gecco_interface is loaded otherwise testmode is assumed
try:
    LABEL
except NameError: 
    for i in _g_argument_names:
        exec(i +"='"+i+"'" )
    def quit_error(msg):
        print msg
        print(format_exc())
_g_required_args=[LABEL,OP_RES,OPERATORS,IDX_SV]




def _remove_whites(string):    
    """removes all (Unicode) white space characters"""
    return re.sub("\s", "",string,flags=re.U)


#----------------------------------------------------------------------------------------
#Classes
#----------------------------------------------------------------------------------------




class _NumberCollectUtil(object):
    """Utility class to derive a factor from operator prefactors"""

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
            raise KeyError("set_idx_sv: 'OPERATORS' not set")
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
        ret=""
        if self.arguments[FAC_INV]!=1:
            ret=str(self.arguments[FAC])+"/"+str(self.arguments[FAC_INV])+"*"
        elif self.arguments[FAC]!=1:
            ret=str(self.arguments[FAC])+"*"
        return ret + "*".join(self.arguments[OPERATORS])

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
#        print(string.display_error)
        self.set_members()
        print "debug: "+string.display_error
        print "debug: Bracket entered"
        self.process_string(string)
        print "debug: Bracket entered"
        self.extract()

    #Functions of the process - extract mechanic
    # see stf_string_expansion for details
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
        print(self.content.extract())
        self.OP_products=map( _OPProduct, self.content.extract())

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
            ret_list.append( prod._set_rule(ges_dic,settle=settle) )
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

    def _combine_dicts(self,form_dict={}):
        """Combines the arguments dictionary with another dict

        The foreign dictionary takes precedence, 
        except in for the FAC and FAC_INV keys. their arguments are multiplicated 
        """
        return combine_dicts(form_dict,self.arguments)

#    def scan 
    
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
        self._content=[]
        self.arguments={LABEL:None,OP_RES:None,NEW:True, FIX_VTX:True}
        self._string=InputString("")

    def __init__(self):
        self.set_members()

    def preprocess_string(self,string):
        """ prepares the input for processing
        
        Removes whitespaces and converts string to InputString
        @param string string to be preprocessed
        @return preprocessed string 
        """
        #remove all white spaces
        string=_remove_whites(string)
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
            print "debug: nth bracket used"
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
    def set_rule(self,special_dic={},settle=True):
        """Invokes set_rule of the brackets
        
        @param special_dic dictionary of touples that may overwrite the self.arguments of this instance (default: empty)
        @param settle for internal use only
        @return None (the combined dicts if settle is False)
        """
        ges_dic=self._combine_dicts(special_dic)
        new=ges_dic.get(NEW,True)
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
        """
        if isinstance(other, _Formula):
            print "it's a formular'"
            self._content+=other._content
            self._string.extend(str(other._string))
        elif isinstance(other, basestring):
            print "debug: it's a string'"
            interm=_Formula()
            print "debug: formula build"
            #relying on the parent class, not GenFormula, 
            #  to make Formula independent from other interface classes 
            string=interm.preprocess_string(other)
            print "preprocessed"
            print string.display_error
            interm.process_string(string)
            print "processed"
            interm.extract()
            self._append(interm)
            print "debug: formula build"
        elif isinstance(other,list):
            for member in other:
                self._append(member)
        else:
            raise NotImplementedError

    def append(self,other):
        """Appending a string or a formula

        @parameter other either a formula (only the body is appended) 
        or a  string representing a formula body """
        print "starting to append"
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

    def __init__(self,label=None,OP_res=None,body=None):
        """The constructor

        @param label LABEL of the formula (can be None)
        @param OP_res OP_RES of this formula (can be None)
        @param body body of the formula"""
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

    def __init__(self,string):
        """The constructor
        
        @param string string of the form 'LABEL:OP_RES=body' """
        #        try:
        self.set_members()     
        self._string=self._preprocess_string(string)
        self.process_string(self._string)
        self.extract()
#        except Exception as ex:
#            quit_error("string_to_form.py:"+ex.msg)
#            raise ex
        
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
        
    def _preprocess_string(self,string):
        """Preprocesses the input
        
        Removes any whitespaces (including tabs)
        Extracts Formula Label and OP_res
        converts string to InputString
        @param string string to be preprocessed
        @return preprocessed string
        """
        string=_remove_whites(string) 
        return self._extract_OP_res(self._extract_label(InputString(string)))
        
        




def make_formula(string):
    thisformula=Formula(string)
    thisformula.set_rule()








class Test(object):
    def test_OP_product(self,verbose=True):
        if verbose:
            print("testing _OPProduct")
        try:
            i=_OPProduct(["C0","H0","T1"])
        except( Exception):
            print "There was an error at creation of _OPProduct" 
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
            print "there was a problem with _OPProduct._show()"\
                  "\nprobably a problem in the dictionary generation"\
                  "\n result: " +str(i.show(inp))+\
                  "\n expected res:" +str(exp_res)
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


    def test_Bracket(self,verbose=True):
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
        

        string="<C0^+*(H+1/2V)C0>"
        i=factory(string)
        try:
            assert str(i) == '<H>'
        except Exception:
            print("Warning _OPProduct.__str__ seems not to be working")
        else:
            if verbose:
                print("all test of _Bracket were successful")


    def test_formula(self,verbose=True):
        inp={"LABEL":"F_LAG","OP_RES":"L"}
        exp_res=[{'OP_RES': 'L', 'FIX_VTX': True, 'OPERATORS': ['H'], 'LABEL': 'F_LAG', 'FAC_INV': 1, 'FAC': 1, 'NEW': True, 'IDX_SV': [1]},
{'OP_RES': 'L', 'FIX_VTX': True, 'OPERATORS': ['H', 'T1'], 'LABEL': 'F_LAG', 'FAC_INV': 1, 'FAC': 0.5, 'NEW': False, 'IDX_SV': [1, 2]},
{'OP_RES': 'L', 'FIX_VTX': True, 'OPERATORS': ['T1', 'H'], 'LABEL': 'F_LAG', 'FAC_INV': 1, 'FAC': -0.5, 'NEW': False, 'IDX_SV': [1, 2]},
{'OP_RES': 'L', 'FIX_VTX': True, 'OPERATORS': ['H', 'T1', 'T1'], 'LABEL': 'F_LAG', 'FAC_INV': 6.0, 'FAC': 1.0, 'NEW': False, 'IDX_SV': [1, 2, 3]},
{'OP_RES': 'L', 'FIX_VTX': True, 'OPERATORS': ['T1', 'H', 'T1'], 'LABEL': 'F_LAG', 'FAC_INV': 6.0, 'FAC': -1.0, 'NEW': False, 'IDX_SV': [1, 2, 3]},
{'OP_RES': 'L', 'FIX_VTX': True, 'OPERATORS': ['T1', 'H', 'T1'], 'LABEL': 'F_LAG', 'FAC_INV': 6.0, 'FAC': -1.0, 'NEW': False, 'IDX_SV': [1, 2, 3]},
{'OP_RES': 'L', 'FIX_VTX': True, 'OPERATORS': ['T1', 'T1', 'H'], 'LABEL': 'F_LAG', 'FAC_INV': 6.0, 'FAC': 1.0, 'NEW': False, 'IDX_SV': [1, 2, 3]}]


        i=_Formula()
        i._string=i.preprocess_string("<H>+0.5<[H,T1]>+1/6*<[[H,T1],T1]>")
        i.process_string(i._string)
        assert i.show(inp) == exp_res

        print("all tests of _Formula were successful")       

        #currently only one test for _Formula. 
    def test_Formula(self,verbose=True):
        print("Formula is currently not tested")
    

    def test_all(self,verbose=True):
        #if arguments are not defined as variables, define them
	self.test_OP_product(verbose)
        self.test_Bracket(verbose)
        self.test_formula(verbose)
        self.test_Formula(verbose)


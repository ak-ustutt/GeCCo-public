import re #support for Regular expressions
import itertools as it # some iterators
import traceback as tb

from gecco_interface import *

##@module A module to convert string representations of formulas to formula objects
#
#Allows to create formula objects and to set them as multiple EXPAND_OP_PRODUCT from a string
#\author Arne Bargholz
#\date 21.05.2015
#\vaersion 1.0


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



#Declare classes
#I have circular dependencies (e.g. parenteses as summands in a sum in parentheses ... )
#Instead of declaring some class I just declare all Rep and base classes
    #better to remember which classes are declare

class _StringConvUtil(object):
    pass
class _ExpFormBaseClass(_StringConvUtil):
    pass


class _SinglePartOperation( _ExpFormBaseClass):
    pass
class _DoublePartOperation( _ExpFormBaseClass):
    pass
class _MultiPartOperation( _ExpFormBaseClass):
    pass


class _InitOperation(_ExpFormBaseClass):
    pass
class _UInitOperation(_ExpFormBaseClass):
    pass


class _OPRep(_InitOperation,_SinglePartOperation):
    pass
class _OPRepAlt(_OPRep):
    pass

class _AddRep(_UInitOperation,_MultiPartOperation):
    pass

class _ComRep( _InitOperation,_DoublePartOperation):
    pass

class _MultRep( _UInitOperation,_DoublePartOperation):
    pass
class _MultRepAlt(_MultRep):
    pass

class _ParRep( _InitOperation,_SinglePartOperation):
    pass

class _BracketRep( _InitOperation,_SinglePartOperation):
    pass



class _BracketUtil( object):
    pass

#******************************************************************************************
#Exceptions
class ExpFormError(Exception):
    #base class for all exceptions special to this module
    pass

class DevError(ExpFormError):
    #class for all Errors in Develoment
    pass

class NotImplementedError(DevError):
    pass

class InputError(ExpFormError):
    # For all bugs on the other side of the monitor
    pass

class IllegalCharError(InputError):
    pass

class IllegalFirstCharError(IllegalCharError):
    pass

class MissingCommaError(IllegalCharError):
    pass

class TooManyCommaError(IllegalCharError):
    pass


#not implemented
#class UnclosedCommutatorError(IllegalCharError):
#    pass

#class UnclosedBracketError(IllegalCharError):
#    pass

#class UnclosedParenthesisError(IllegalCharError):
#    pass

class UnexpectedEndOfStringError(InputError):
    pass

class BadSplitError(InputError):
    pass

class MissingArgumentError(InputError):
    pass

class MissingLabelError(MissingArgumentError):
    pass

class MissingOP_ResError(MissingArgumentError):
    pass


#######################################################################
#global functions

def _pow_sep(context):
    """ returns the possible separators of the decimal power
    
    Returns the separators either as RegExp string "[...]" or as tuple
    @param context ="reg_exp" "!reg_exp" (inverted RegExp) or "tuple"
    """
    pow_sep=('d','e')
    if context=="reg_exp":
        return "["+"".join(pow_sep)+"]"
    elif context =="!reg_exp":
        return "[^"+"".join(pow_sep)+"]"
    elif context =="tuple":
        return pow_sep
    else:
        raise NotImplementedError("context "+str(context)+"not implemented")

def _num_reg_exp_total(context="start"):
    """just a function to store and explain the total regexp used to scan for numbers"""
    d=_pow_sep("reg_exp")
    reg_exp1="-?[0-9]*"+    "\.[0-9]+"  +"(?:"+d+"-?[0-9]+)?"
    reg_exp2="-?[0-9]+"+ "(?:\.[0-9]*)?"+"(?:"+d+"-?[0-9]+)?"
    #as described in the re libary
    # "-?" matches zero or one "-"
    # "[0-9]*" matches zero or more digits "[0-9]+" matches one or more digits
        #for internationality this might be changed to "\d"
    # "\." matches a "." backslash for escaping
    # "(?:"  and ")" groups an expression combined with the tailing "?" 
        #this matches zero or one times this group of expressions      
    # more digit matching 
    # more "(?:" ")" and "?"
    # d returns a string that is the reg_exp for the power separator so "[de]" by default
    # matching anoter group of positive or negative integers
    # the connection by "|" let us have the longer match by either RegExp1 or RegExp2
    # "^" matches the beginning of a string
    # getting it all in Parentheses tells re.split to also return the matched expression
    # So this regexp matches every of the following number formats: 1 , 1. 1.234 , .234
        # all cases can have a preceeding minus or a following exponent(which in turn can 
        #be negative)
    if context=="start":
        return "(^"+reg_exp1+"|^"+reg_exp2+")"
    elif context=="anywhere":
        return "("+reg_exp1+"|"+reg_exp2+")"
    elif context=="bare":
        return reg_exp1+"|"+reg_exp2
    else:
        raise Exception("_num_reg_exp_total: unkown context")
    


def _combine_dicts(first,second):
    """concatenates two dictionaries

    if both have an identical key. the first dictionary takes precedence
    """

    def _combine_dicts_core(alpha,beta,x):
        if x in alpha:
            return alpha.get(x)
        else:
            return beta.get(x)
    return {x: _combine_dicts_core(first,second,x)  for x in set(first).union(second)}


######################################################################
#global objects 

_number_reg_exp=re.compile(_num_reg_exp_total())
#regexpused to scan for numbers at the start of a string (e.g an operator)
#any match bby _NumberRegExp is (apart from _pow_sep) avalid float 

_special_reg_exp=re.compile("(^-)")
#matches anything that begins with a  minus
#second regexp like numbers, but has to be treated as -1 

_special_bch_reg_exp=re.compile("BCH",flags=re.I)
#ignore case, but if you ever type bCh, I will hate you


_div_reg_exp=re.compile("("+_num_reg_exp_total(context="bare")+
                           " */ *"
                           +_num_reg_exp_total(context="bare")+ ")" )
# matches any division of numbers
# 


###############################################
#Classes

class _InputStringBase(object):

    def __init__(self,string):
        """The constructor"""
        self.set_member()
        self.string=string

    #functions that change the internal state
    def set_start(self,start=0):
        """Sets self.begin to this position
        
        Later a substring can be extracted beginning at this postition
        """
        self.begin=self.i+start

    def goto(self,character):
        """Jumps to character
        
        This method jumps to the next occurence of character.
        this method does not jump to a specific index.
        """
        self.i=self.string.find(character,self.i)
        return self.i
        ## \todo  something. I forgot. The world is fading

    def reset(self,i=0,begin=0,toggle=False):
        """ Resets the counters

        For testing
        """
        self.begin=begin
        self.i=i
        self.toggle_advance_on_minus=toggle

    # Overwrite of builtins
    def __len__(self):
        """allows the use of len(string)"""
        return len(self.string)

    def __str__(self):
        """makes object printable"""
        return self.string

    #display functions (as property getters)
    @property
    def display_error(self): 
        """returns the string and a vizualization of the current character
        
        will get called by the Exception handlers"""
        return "\n"+self.string+"\n"+"-"*self.i+"^"


    @property
    def has_not_ended(self):
        """tests if end of string is reached"""
        return self.i<len(self.string)
        #for while string.has_not_ended

    @property
    def substring_number_match(self):
        """test if current position is in the number part of an OP"""
        #If current position is within the number part of an Operator starting at self.begin
        #  return True
        #  else return False
        for reg_exp in [_number_reg_exp,_special_reg_exp]:
            match=reg_exp.match(self.string[self.begin:])
            if match and  match.end() >= (self.i-self.begin):
                return True
        return False


    @property
    def substring(self):
        return self.string[self.begin:self.i]

class _HiddenMinusUtil(_InputStringBase):
    """Utility class with methods that hide the minus from objects

    The true_xxx methods are also defined here
    """
    def set_member(self):
        """sets the variables of the class that store information"""
        self.string=''
        self.i=0 #"pointer" to current char
        self.begin=0  
        #"pointer" to beginn of current string. is updated at beginn of an object
        self.toggle_advance_on_minus=False
        #next() is called when an object has processed a character
        # minus signs must be processed by two objects (_AddRep and _OPRep)
        # therefor only every second something should advance over a minus the 
        # next() really will advance

    def next(self):
        """advances to next char
        
        Sets pointer to next char3
        """
        if self.true_this_char != "-":
            self.i +=1
        else:
            if self.toggle_advance_on_minus:
                self.i+=1
            self.toggle_advance_on_minus=not self.toggle_advance_on_minus

    #display only functions as property getters
    @property
    def this_char(self):
        """returns the current character hides minus"""
        if self.string[self.i] != "-":
            return self.string[self.i]
        else:
            return "+"
        #it's basically the same

    @property
    def last_char(self): 
        """returns the previous character hides minus"""
        return self.string[self.i-1]
        #if this is a minus you ended your operator ended on(was a) minus

    @property
    def next_char(self): 
        """returns the next character hides minus"""
        #notice: self.next_char only returns something.
        # self.next() increases the pointer and returns nothing
        if self.string[self.i+1] != "-" :
            return self.string[self.i+1]
        else:
            return "+"

    #and the implementations that don't hide "-"
    @property 
    def true_this_char(self):
        """returns the current character"""
        return self.string[self.i]

    @property
    def true_last_char(self): 
        """returns the previous character"""
        return self.string[self.i-1]

    @property
    def true_next_char(self): 
        """returns the next character"""
        return self.string[self.i+1]


class _InputString(_HiddenMinusUtil):
    """class to handle the string that is converted to the formula"""

    def __init__(self,string):
        """The constructor"""
        self.set_member()
        self.string=string

class _InputStringAlt(_InputString):
    """advanced constructor for testing"""
    def __init__(self,string,i=0,begin=0,toggle=False):
        self.set_member()
        self.string=string
        self.reset(i,begin,toggle)

#########################
#Classes for the string to expanded (without addition) 

class _StringConvUtil(object):
    """Store all functions helping to convert the string to an object hierachy

    
    """
    @staticmethod
    def check_illegal_char(string,tuple,context='not first'):
        """Tests if the current char is forbidden at this point
        
        raises an exception if the current character is in tuple. 
        @param tuple tuple should be tuple or list
        """
        if string.true_this_char in tuple:
            if context=="first":
                raise IllegalFirstCharError(string.true_this_char,
                                            " may not come directly after ",
                                            string.true_last_char,
                                            string.display_error,)
            ##\todo implement a switch for different objects
            else :
                raise IllegalCharError("Illegal character detected:",
                                       string.true_this_char,
                                       string.display_error,)

    @staticmethod
    def check_close_char(string,tuple):
        """Tests if the current char closes the object."""
        #
        if string.this_char in tuple:
            return True
        else:
            return False 

    def _create_next_operation(self,string):
        """Creates and returns a lower tier object"""
        ## \todo find a better name, that doesn't imply case is a loop

        #intended order of use:
        #self.check_illegal_char( ... )
        #if self.check_close_char(...)
        #    close
        #elif "special cases"
        #    ...
        #else:
        #    self._create_next_operation(string)
        if string.this_char=="[":
            return _ComRep(string)
                    
        elif string.this_char=="(":
            return _ParRep(string)
    
        elif string.this_char=="*":
            return _MultRep(string,self.content[self.part])

        elif string.this_char=="(":
            return _AddRep(string,self.content[self.part])
        else:
            return _OPRep(string)




class _ExpFormBaseClass(_StringConvUtil): 
    """
    Abstract class for operator containing objects
    This class has only prototypes of functions
    """

    def set_members(self):
        # sets content and part 
        self.content=[]
        #contains a list of objects joined by the operation the class represents
        self.part=0
        #pointer to current object in list   
        raise NotImplementedError()

    def process_string(self,string):
        #creates content from string
        raise NotImplementedError()
    
    def close(self):
        #method that is called, when the string processing finishes
        #part of an old architecture. does nothing now.
        raise NotImplementedError()

    def extract(self):
        #returns a sum of Operator products (a list of lists)
        raise NotImplementedError()

class _SinglePartOperation( _ExpFormBaseClass):
    def set_members(self):
        self.content=[""]
        self.part=0

class _DoublePartOperation( _ExpFormBaseClass):
    def set_members(self):
        self.content=["",""]
        self.part=0

class _MultiPartOperation( _ExpFormBaseClass):
    def set_members(self):
        self.content=[""]  
        # For _AddRep : the first object is stored by __init__ 
        # all other objects are appended
        self.part=0


class _InitOperation(_ExpFormBaseClass):
    """Abstract class for all operations with an initializing command
    
    These operations are "(" "[" and operators
    """
    def __init__(self,string):
        self.set_members()
        self.process_string(string)

class _UInitOperation(_ExpFormBaseClass):
    """ Abstract class for all operations without an initializing command

    These operations are:"+" "*"
    """
    def __init__(self,obj,string):
        self.set_members()
        self.content[self.part]=obj
        self.part=1
        self.process_string(string)



class _OPRep(_InitOperation,_SinglePartOperation):
    """Class for operators """

    def process_string(self,string):
        string.set_start()
        self.check_illegal_char(string, ("<",">","+",""),"first") 
        ##\todo think about which chars are illegal at the beginning of objects 
        while string.has_not_ended:

            self.check_illegal_char(string,("<"))            
            if ( string.true_this_char == "+" and string.last_char == "^" ):
                #Ignore ^+ combinations
                string.next()
            elif ( string.true_this_char=="-" and string.substring_number_match ): 
                string.next()
                #ignore "-" that belong to a number -1.0d-4
                #                                  -^----^
                # thanks to toggle minus advance this will be done twice on the second minus
            elif self.check_close_char(string,("[",",","]","(",")","*","+",">","-")):
                self.content=[string.substring]
                self.close()
                return 0
            else:
                string.next()
        raise UnexpectedEndOfStringError(string.display_error)

    def close(self):
        return 
    
    def extract(self):
        return [self.content]
        #Since self.content is also a list (of one object)
        # this returns a list of lists like any other Rep object, allowing recursion.

class _OPRepAlt(_OPRep):
    """Alternative class for operators,
 
    when you only have to convert a str to an operator
    """
    def __init__(self,str_):
        self.set_members()
        self.content=[str(str_)]


class _AddRep(_UInitOperation,_MultiPartOperation):
    """Class to store all summands of an addition"""


    def process_string(self,string):
        self.part-=1
        while string.has_not_ended:
            if string.this_char=="+" and string.true_next_char =="-" :
                # can occur with input: T1+-0.5T2 
                string.next()
            elif string.this_char=="+":
                string.next()
                self.check_illegal_char(string, ("]",")",",","*","+")) 
                self.part +=1
                #string.i can point neither to * nor to + so part may point to out of scope
                self.content.append(self._create_next_operation(string))


            elif self.check_close_char(string,(",","]",")",">")):
                self.close()
                return 

            elif string.this_char=="*": 
                print self.part,self.content
                self.content[self.part]=_MultRep(self.content[self.part],string)
      
            else:
                #implicit multipication T1+(T2+T3)T4
                # after handling the parethesis the _AddRep sees the T of T4 
                self.content[self.part]=_MultRep(self.content[self.part],string)
        raise UnexpectedEndOfStringError(string.display_error)

    def close(self):
        return 0    

    def extract(self):
        #self.content contains a list of objects. upon extraction those return a list of 
        # additively joined lists of mutiplicatively joined operators
        # to extract those lists have simply to be concatenated
        return reduce(lambda x,y: x+y.extract() ,self.content, [])


        
class _ComRep( _InitOperation,_DoublePartOperation):
    """Class to store the commutator parts"""
#    counter=0
#    @classfunction
#    def increase_counter(cls):
#        cls.counter+=1

    def process_string(self,string):                
        string.next() # String pointed to "["

        self.check_illegal_char(string, ("<",">","*","","+"),"first")
        self.content[self.part]=self._create_next_operation(string)
        while string.i <len(string):
            self.check_illegal_char(string,("<",")")) 
            
            if string.this_char==",":
                if self.part==0: # If not, it is a [,, combination
                    self.part=1
                    string.next()
                    self.check_illegal_char(string,("<",")","",",","]"))
                    self.content[self.part]=self._create_next_operation(string)
                else:
                    raise TooManyCommasError("There may only be one ',' ",
                                             "in any commutator",
                                             string.display_error)

            elif self.check_close_char(string,("]")):
                if self.part==1: # If not, it is a [] combination
                    self.close()
                    string.next()
                    return 
                else:
                    raise MissingCommaError("All commutators are of the form [...,...]",
                                            string.display_error)

            elif string.this_char=="+": 
                self.content[self.part] = _AddRep(self.content[self.part],string)

            else:
                #implicit multipication T1+(T2+T3)T4
                # after handling the parethesis the AddRep sees the T of T4 
                self.content[self.part]=_MultRep(self.content[self.part],string)
        raise UnexpectedEndOfStringError(string.display_error)

    def close(self):
        return 0

    def extract(self):
        #Transform [O1,O2] -> O1*O2 + ("-1"*(O2*O1)) then extract those and concatenate them 
        #(see addition)
        interm1 = _MultRepAlt(self.content[0],self.content[1])
        interm2 = _MultRepAlt(_OPRepAlt("-1"),_MultRepAlt(self.content[1],self.content[0]))
        return interm1.extract()+interm2.extract()
        #Yeah, i could do this without intermediates. Do you want me to push that in one line?



class _MultRep( _UInitOperation,_DoublePartOperation):
    """Class to store two Factors of a multiplication
    
    _MultRep only handles two factors at a time to allow for easy extraction
    """

    def process_string(self,string,obj=''):
        if string.this_char=="*": 
            #not necessarily true since (F+V)(T1+T2) also triggers a Multiplication
            string.next()

        self.check_illegal_char(string, (")","]","*",",","","+") )
        
        self.content[self.part]=self._create_next_operation(string)

        while string.has_not_ended:
            self.check_illegal_char(string,("<")) 

            if self.check_close_char(string,(",","]",")",">","+")):
                self.close()
                return 

            else:
                self.content[self.part]=_MultRep(self.content[self.part],string)

        raise UnexpectedEndOfStringError(string.display_error)

    def close(self):
        if self.part != 1:
            raise ToManyPartsError("this is a _MultRep object with "+str(self.part) +"parts" )
        return 

    def extract(self):
        return [x+y for x in self.content[0].extract() for y in self.content[1].extract()]


class _MultRepAlt(_MultRep):
    """ When you already have two objects an only need to multiply them. """
    #look at _ComRep
    def __init__(self,obj1,obj2):
        self.set_members()
        self.content=[obj1,obj2]
        self.part=1

        
class _ParRep( _InitOperation,_DoublePartOperation):
    """Class to store the content of a parenthesis"""

    def process_string(self,string):
        string.next()
        self.check_illegal_char(string, ("<",">","*","]",")","","+"),"first") 
        self.content[self.part]=self._create_next_operation(string)

        while string.has_not_ended:
            self.check_illegal_char(string,("<",",","]")) 

            if self.check_close_char(string,(")")):
                self.close()
                string.next()
                return 

            elif string.this_char=="+": 
                self.content[self.part] = _AddRep(self.content[self.part],string)

            else: 
                self.content[self.part] = _MultRep(self.content[self.part],string)

        raise UnexpectedEndOfStringError(string.display_error)

    def close(self):
        return 0

    def extract(self):
        return self.content[0].extract()



class _BracketRep(_InitOperation,_SinglePartOperation):
 
    def process_string(self,string):
        string.next()
        self.check_illegal_char(string,("<",",","]",")","")) #
        
        self.content[self.part]=self._create_next_operation(string)

        while string.has_not_ended:
            self.check_illegal_char(string,("<",",","]",")","")) # 

            if string.this_char=="+": 
                self.content[self.part]=_AddRep(self.content[self.part],string)

            elif self.check_close_char(string,(">")):
                self.close()
                string.next()
                return 

            else: 
                self.content[self.part]=_MultRep(self.content[self.part],string)

        raise UnexpectedEndOfStringError(string.display_error)
        
    def close(self):
        return  

    def extract(self):
        return self.content[0].extract()


##########################################################################################
#the surface classes: _ExpBraketRep and _Formula 
# _ExpBracketRep is for already 
class _NumberCollectUtil(object):
    """Utility class to derive a factor for bracket from operator prefactors"""

    def process_parts(self,parts,OP_no):
        """processes the output of a re.split operation"""
        ##@return returns 0 if no (zero) splits occured
        #         returns 1 if a split occurred
        if len(parts) not in (1,3):
            raise BadSplitError("Operator"+str(self.operators[OP_no])+
                                "was split into "+str(len(parts))+" parts")
        elif len(parts) == 1:
            return 0
        else:
            self.operators[OP_no]=parts[2]
            if _number_reg_exp.match(parts[1]):
                self.factor=self.factor*float(re.sub(_pow_sep("reg_exp"),"e",parts[1])) 
            elif _special_reg_exp.match(parts[1]):
                self.factor=-self.factor
            return 1

    def split(self,reg_exp,OP_no):
        """ Splits an operator into """
        return self.process_parts(reg_exp.split(self.operators[OP_no]),OP_no)
        #if there was a match:
        #parts[0] is an empty string 
            #because our strings match the beginning. read documentation of re
        #parts[1] contains our match 
        #parts[2] contains the rest (at least an empty string)
        # if there is no match:
        # parts[0] contains our string

        ##@return returns the return of self.process parts


    def pop_numbers(self):
        """Collects the numbers from operators
        
        Find all numerical prefixes,
        multiply the numbers, eliminate the numbers 
        from the elements and delete those elements that become empty
        return the product
        """
        i=len(self.operators)
        while i>0:
            i-=1
            #going backwards to not mess up indices
            #a for loop may be more intelligent
            
            if self.split(_number_reg_exp,i):
                #splits the ith operator into number and residue
                # if there was a split continue.
                pass
            else:
                #else try the special regexp
                self.split(_special_reg_exp,i)

            if (self.operators[i]==''):
                del self.operators[i]

class _IDXUtil(object):
    """collection of methods to create the idx_sv list"""
    
    def set_idx_sv(self):
        """function to generate a idx_sv list

        """
        counter=it.count(1)
        #generator that returns 
        special_OPs=[]
        special_OP_numbers=[]
        #two lists. 
        #the first contains the operators with ticks ("'") 
        #the second one their respective idx_indices.
        ##\todo make them into one class
        for OP in self.operators:
            if OP.find("'")>=0:
                if OP.replace("'","") in special_OPs:
                    index=special_OPs.index(OP.replace("'",""))
                    self.idx_sv.append(special_OP_numbers[index])
                else:
                    special_OPs.append(OP.replace("'",""))
                    special_OP_numbers.append(counter.next())
                    self.idx_sv.append(special_OP_numbers[-1])
            else:
                self.idx_sv.append(counter.next())



class _BracketBase(_IDXUtil,_NumberCollectUtil):
    def __init__(self):
        self.set_members()

    def set_members(self):
        self.factor=1.0
        self.operators=[]
        self.idx_sv=[]



class _ExpBracket(_BracketBase):
    """Class for expanded Brackets (contain only OP_Products)"""
    def __init__(self,OpProd):
        self.set_members()
        self.operators=OpProd[:]
        self.pop_numbers()
        self.set_idx_sv()

    def __str__(self):
        return str(self.factor)+"*<"+"*".join(self.operators)+">"

    def _strip_operators_from_tics(self):
        """returns the operators without single quotes"""
        return map(lambda x: re.sub("'","",x),self.operators)


    def dic(self,form_dic):
        #builds a dictionary for EXPAND_OP_PRODUCT
        dic={}
        dic['OPERATORS']=self._strip_operators_from_tics()
        dic['FAC']=self.factor
        dic['IDX_SV']=self.idx_sv
        return _combine_dicts(dic,form_dic)

    def set_rule(self,form_dic={},settle=True):
        """Invokes EXPAND_OP_PRODUCT"""
        ges_dic=self.dic(form_dic)
        for argument in ("LABEL", "OP_RES", "OPERATORS", "IDX_SV"):
            if argument not in ges_dic:
                raise MissingArgumentError(argument+"is missing for EXPAND_OP_PRODUCT of "+str(self))
        if settle:
            EXPAND_OP_PRODUCT(ges_dic)
        else:
            return ges_dic

    def show(self,form_dic={}):
        return self.set_rule(form_dic,settle=False)

##########################################################################################
#the surface class Formula for handling formula on the python surface

class Formula(object):
    """Base class for all User interfaces"""


    def set_members(self):
        self.content=[]
        self.part=0
        self.string=""
        self.new=True
        self.OP_res=""
        self.label=""

    def dic(self,special_dic):
        """ creates a dictionary of the member variables
 
        the dictionary is the basis for all EXPORT_OP_PRODUCT dicts of this formula
        """
        dic={}
        dic['NEW']=self.new
        dic['OP_RES']=self.OP_res
        dic['LABEL']=self.label
        return _combine_dicts(special_dic,dic)

    def __init__(self):
        self.set_members()

    def preprocess_string(self,string):
        """ prepares the input for processing
        
        Removes whitespaces and
        saves the input as _InputString in self.string
        """
        string= re.sub(" ", "",string)
        #remove all white spaces
        self.string=_InputString(string)

    def process_string(self):
        while self.string.has_not_ended:
            self.string.goto("<")
            self.content.append(_BracketRep(self.string))
            
        
    def extract(self):
        interm=list(reduce(lambda x,y: x+y.extract() , self.content, [] ))
        self.content=[]
        for item in interm:
            self.content.append(_ExpBracket(item))

#    def show(self,form_dic={}):
#        return map(lambda item: item.show(form_dic), self.content)
        
    def set_rule(self,special_dic={},settle=True):
        new=special_dic.get('NEW',self.new)
        ret=[]
        for item in self.content:
            special_dic['NEW']=new
            if settle:
                item.set_rule(self.dic(special_dic),settle=settle)
            else:
                ret.append( item.set_rule(self.dic(special_dic),settle=settle))
            new=False
        return ret

    def show(self,special_dic={}):
        """Returns the dictionaries that are given to EXPAND_OP_PRODUCT"""
        return self.set_rule(special_dic=special_dic,settle=False)

    def str_short(self):
        """Returns the original input as string"""
        return str(self.string)
        
    def __str__(self):
        return self.strings[0]+"\n+ ".join(self.strings[1:])

    def strings(self):
        """Returns the formula as list of brackets 
        
        
        """
        beginning= [self.label+":"+self.OP_res+"="]
        body=[str(bracket) for bracket in self.content]
        return beginning+body

#User interfaces (classes and functions)
#
class Formula1(Formula):
    """First actual user interfacing class"""
    def extract_label(self):
        """Extracts the formula label from the input string"""
        self.string.set_start()
        if self.string.goto(":")<0:
            raise MissingLabelError("Please precede your formula by 'my_label:' ")
        self.label=self.string.substring
        if len(self.label)==0:
            raise MissingLabelError("This formula label is empty")
        self.string.next()
        
    def extract_OP_res(self):
        """Extracts the result operators label from the input string"""
        self.string.set_start()
        if self.string.goto("=")<1:
            raise MissingOP_ResError("Please precede your formula by 'LABEL:OP_RES=' ")
        self.OP_res=self.string.substring
        if len(self.OP_res)==0:
            raise MissingOP_ResError("This result Operator is empty")
        self.string.next()
        
    def preprocess_string(self,string):
        """Preprocesses the input
        
        Removes any whitespaces (including tabs)
        Extracts Formula Label and OP_res
        Puts input in self.string (as _InputString)
        """
        #remove all (Unicode) white space characters
        string= re.sub("\s", "",string,flags=re.U)
        #set string as _InputString
        self.string=_InputString(string)
        self.extract_label()
        self.extract_OP_res()

    def __init__(self,string):
        """The constructor"""
        self.set_members()     
        self.preprocess_string(string)
        self.process_string()
        self.extract()




def make_dicts(string):
    thisformula=Formula1(string)
    thisformula.show()

def make_formula(string):
    thisformula=Formula1(string)
    thisformula.set_rules()













#---------------------------------------------------------------------------------
#Testing utility
class Test(object):
    """class to access the hidden classes even when importing"""
    def set_member(self):
        self.helpful="undefined"
        #as soon as set this should be a boolean
    @staticmethod
    def make_string(string):
        return _InputString(string)
    @staticmethod
    def make_formula(string):
        return _Formula(string)
    @staticmethod
    def make_bracket(string):
        return _BracketRep(string)
    @staticmethod
    def make_parenthesis(string):
        return _ParRep(string)
    @staticmethod
    def make_mult(string):
        return _MultRep(string)
    @staticmethod
    def make_commutator(string):
        return _ComRep(string)
    @staticmethod
    def make_add(string):
        return _AddRep(string)
    @staticmethod
    def make_OP(string):
        return _OPRep(string)
    @staticmethod
    def make_Formula(string):
        return _Formula(string)
    def ask_helpful(self):
        """asks user to rate the errormessage as (un-)helpful
        
        If self.helpful is not set asks the user above question
        """
        if self.helpful=="undefined":
            while True:
                answer=raw_input("Was this error-message helpful? (y/n)")
                if answer in ("y","Y","yes","Yes","YES"):
                    return True
                elif answer in ("n","N","no","No","NO"):
                    return False
        else: 
            return bool(self.helpful)
            # Hopefully no idiot put something besides True or False into self.helpful
            #especially not "False"
        

    def try_(self,string,exp_res):
        """core function for all testing functions

        can only administer tests for  Rep objects, deriving from _Init
        """

        cls=self.test_class
        try:
            print "input:" +str(string)
            if issubclass(cls ,_UInitOperation):
                handle=cls(self.OP1,string)
            else:
                handle=cls(string)
            print "Sucessfully created "+str(handle)
            #res for result
            res=handle.extract()
            if res==exp_res:
                print "the object was sucessfully extracted:"
                print res
                return 1
            else:
                print "Oh, no extraction didn't work as planned"
                print "expected: " + str(exp_res)
                print "result: " + str(res)
                return 0
        except(Exception):
            print tb.print_exc()

    def test_OPRep(self):
        """ testing _OPRep with valid strings"""
        print self.test_OPRep.__doc__
        
        #self.test_class is a variable that tells try_ which class to test
        self.test_class=_OPRep
        error_counter=0

        #Shortening names (Ts for Test string)
        try_=self.try_
        Ts=_InputStringAlt
        

        #1-try because try gives back 0 or 1 and one (true) 
        # should correspond to a successfull test
        error_counter+=1-try_(Ts("T3)"),[["T3"]])
        error_counter+=1-try_(Ts("3T_evy)"),[["3T_evy"]])
        error_counter+=1-try_(Ts("-GAM)"),[["-GAM"]])
        error_counter+=1-try_(Ts("1d-1T2)"),[["1d-1T2"]])
        error_counter+=1-try_(Ts(".5Gam_g)"),[[".5Gam_g"]])
        error_counter+=1-try_(Ts("C0^+)"),[["C0^+"]])
 



        error_counter+=1-try_(Ts("T3+"),[["T3"]])
        error_counter+=1-try_(Ts("T3^++"),[["T3^+"]])
        error_counter+=1-try_(Ts("T3,"),[["T3"]])
        error_counter+=1-try_(Ts("T3("),[["T3"]])
        error_counter+=1-try_(Ts("T3["),[["T3"]])
        error_counter+=1-try_(Ts("T3]"),[["T3"]])
        print "Testing operator creation, we encountered ",error_counter," errors."        
        return error_counter

    def test_MultRep(self):
        """ testing _MultRep with valid strings"""
        print self.test_MultRep.__doc__
        #self.test_class is a variable that tells try_ which class to test
        self.test_class=_MultRep
        error_counter=0

        #Shortening names (Ts for Test string)
        try_=self.try_
        Ts=_InputStringAlt

        # And another variable used to communicate with self.try_
        self.OP1=_OPRepAlt("T1")

        error_counter+=1-try_(Ts("*T3+"),[["T1","T3"]])
        error_counter+=1-try_(Ts("*3T3+"),[["T1","3T3"]])
        error_counter+=1-try_(Ts("*-T3+"),[["T1","-T3"]])
        error_counter+=1-try_(Ts("*1d-1T3+"),[["T1","1d-1T3"]])
        error_counter+=1-try_(Ts("*.5T_3+"),[["T1",".5T_3"]])
        error_counter+=1-try_(Ts("*T3'^++"),[["T1","T3'^+"]])
        error_counter+=1-try_(Ts("*-T3+"),[["T1","-T3"]])


        error_counter+=1-try_(Ts("T3+"),[["T1","T3"]])
        error_counter+=1-try_(Ts("T3*T2+"),[["T1","T3","T2"]])
        error_counter+=1-try_(Ts("T3)"),[["T1","T3"]])
        error_counter+=1-try_(Ts("T3,"),[["T1","T3"]])
        error_counter+=1-try_(Ts("T3]"),[["T1","T3"]])
        error_counter+=1-try_(Ts("T3-"),[["T1","T3"]])
        
        print "Testing Multiplication, we encountered ",error_counter," errors."
        return error_counter

    def test_AddRep(self):
        """ testing _AddRep with valid strings"""
        print self.test_AddRep.__doc__
        #self.test_class is a variable that tells try_ which class to test
        self.test_class=_AddRep
        error_counter=0

        #Shortening names (Ts for Test string)
        try_=self.try_
        Ts=_InputStringAlt

        # And another variable used to communicate with self.try_
        self.OP1=_OPRepAlt("T1")

        error_counter+=1-try_(Ts("+T3)"),[["T1"],["T3"]])
        error_counter+=1-try_(Ts("+3T3)"),[["T1"],["3T3"]])
        error_counter+=1-try_(Ts("+-T3)"),[["T1"],["-T3"]])
        error_counter+=1-try_(Ts("+1d-1T3)"),[["T1"],["1d-1T3"]])
        error_counter+=1-try_(Ts("+.5T_3)"),[["T1"],[".5T_3"]])
        error_counter+=1-try_(Ts("+T3'^+)"),[["T1"],["T3'^+"]])
        error_counter+=1-try_(Ts("-T3)"),[["T1"],["-T3"]])


        error_counter+=1-try_(Ts("+T3)"),[["T1"],["T3"]])
        error_counter+=1-try_(Ts("+T3]"),[["T1"],["T3"]])
        error_counter+=1-try_(Ts("+T3,"),[["T1"],["T3"]])

        self.OP1=_MultRepAlt(_OPRepAlt("T1"),_OPRepAlt("T2"))
        
        error_counter+=1-try_(Ts("+T3*T4)"),[["T1","T2"],["T3","T4"]])

        print "Testing Addition, we encountered ",error_counter," errors."
        return error_counter



    def test_ParRep(self):
        """ testing _ParRep with valid strings"""
        print self.test_ParRep.__doc__
        #self.test_class is a variable that tells try_ which class to test
        self.test_class=_ParRep
        error_counter=0

        #Shortening names (Ts for Test string)
        try_=self.try_
        Ts=_InputStringAlt

        error_counter+=1-try_(Ts("(T3)"),[["T3"]])
        error_counter+=1-try_(Ts("(3T3)"),[["3T3"]])
        error_counter+=1-try_(Ts("(-0.5T3)"),[["-0.5T3"]])
        error_counter+=1-try_(Ts("(1d-1T3)"),[["1d-1T3"]])
        error_counter+=1-try_(Ts("(.5T_3)"),[[".5T_3"]])
        error_counter+=1-try_(Ts("(T3'^+)"),[["T3'^+"]])
        error_counter+=1-try_(Ts("(-T3)"),[["-T3"]])

        error_counter+=1-try_(Ts("(T1+T3)"),[["T1"],["T3"]])
        error_counter+=1-try_(Ts("(T1*T3)"),[["T1","T3"]])
        error_counter+=1-try_(Ts("(T1+T2*T3)"),[["T1"],
                                                ["T2","T3"]])
        error_counter+=1-try_(Ts("(T1*(T2+T3))"),[["T1","T2"],
                                                  ["T1","T3"]])
        error_counter+=1-try_(Ts("((T1+T2)*(T3+T4)))"),[["T1","T3"],
                                                        ["T1","T4"],
                                                        ["T2","T3"],
                                                        ["T2","T4"]])
        error_counter+=1-try_(Ts("(T1(T2+T3))"),[["T1","T2"],
                                                 ["T1","T3"]])
        error_counter+=1-try_(Ts("((T1+T2)(T3+T4)))"),[["T1","T3"],
                                                       ["T1","T4"],
                                                       ["T2","T3"],
                                                       ["T2","T4"]])
        error_counter+=1-try_(Ts("(T1-T3)"),[["T1"],["-T3"]])

        print "Testing parentheses, we encountered ",error_counter," errors."

        return error_counter

    def test_ComRep(self):
        """ testing _ParRep with valid strings"""
        print self.test_ComRep.__doc__
        #self.test_class is a variable that tells try_ which class to test
        self.test_class=_ComRep
        error_counter=0

        #Shortening names (Ts for Test string)
        try_=self.try_
        Ts=_InputStringAlt

        error_counter+=1-try_(Ts("[T1,T3]"),[["T1","T3"],
                                             ["-1","T3","T1"]])
        error_counter+=1-try_(Ts("[-T1,-T3]"),[["-T1","-T3"],
                                               ["-1","-T3","-T1"]])
        error_counter+=1-try_(Ts("[1.0d-3T1,2.3d-4T3]"),[["1.0d-3T1","2.3d-4T3"],
                                                         ["-1","2.3d-4T3","1.0d-3T1"]])
        error_counter+=1-try_(Ts("[T1+T2,T3]"),[["T1","T3"],
                                                   ["T2","T3"],
                                                   ["-1","T3","T1"],
                                                   ["-1","T3","T2"]
                                               ])

        error_counter+=1-try_(Ts("[T1-T2,T3]"),[["T1","T3"],
                                                   ["-T2","T3"],
                                                   ["-1","T3","T1"],
                                                   ["-1","T3","-T2"]
                                               ])
        error_counter+=1-try_(Ts("[T1,T3*T4]"),[["T1","T3","T4"],
                                                ["-1","T3","T4","T1"]])
        print "Testing commutation, we encountered ",error_counter," errors."
        return error_counter
#        error_counter+=1-try_(Ts("(3T3)"),[["3T3"]])
#        error_counter+=1-try_(Ts("(-0.5T3)"),[["-0.5T3"]])
#        error_counter+=1-try_(Ts("(1d-1T3)"),[["1d-1T3"]])
#        error_counter+=1-try_(Ts("(.5T_3)"),[[".5T_3"]])
#        error_counter+=1-try_(Ts("(T3'^+)"),[["T3'^+"]])
#        error_counter+=1-try_(Ts("(-T3)"),[["-T3"]])
 
    def test_BracketRep(self):
        """ testing _BracketRep with valid strings"""
        print self.test_BracketRep.__doc__
        #self.test_class is a variable that tells try_ which class to test
        self.test_class=_BracketRep
        error_counter=0

        #Shortening names (Ts for Test string)
        try_=self.try_
        Ts=_InputStringAlt

        error_counter+=1-try_(Ts("<T3>"),[["T3"]])
        error_counter+=1-try_(Ts("<(T3)>"),[["T3"]])
        error_counter+=1-try_(Ts("<-T3>"),[["-T3"]])
        error_counter+=1-try_(Ts("<(T1)T3>"),[["T1","T3"]])
        error_counter+=1-try_(Ts("<T1+T3>"),[["T1"],["T3"]])
        error_counter+=1-try_(Ts("<T1-T3>"),[["T1"],["-T3"]])
        error_counter+=1-try_(Ts("<(T1+T2)T3>"),[["T1","T3"],["T2","T3"]])

        print "Testing brackets, we encountered ",error_counter," errors."
        return error_counter

    def run_test(self,helpful="undefined"):
        """ Testsuite for string_to_form V0.8"""
        self.helpful=helpful
        print self.run_test.__doc__
        
        error_counter=0
        error_counter+=self.test_OPRep()
        error_counter+=self.test_MultRep()
        error_counter+=self.test_AddRep()
        error_counter+=self.test_ComRep()
        error_counter+=self.test_BracketRep()
        print "there were ",error_counter," Errors"
        

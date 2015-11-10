from stf_regexp import stf_number_regexp, stf_negative_regexp, stf_div_regexp
# The regular expressions for numbers ... 
from operators import Vertex
#The operator object everything works
from stf_exceptions import *
# custom exceptions



#Instead of declaring some class I just declare all Rep classes and  their base classes
    #better to remember which classes are declared

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

class BracketRep( _InitOperation,_SinglePartOperation):
    pass


#******************************************************************************************
# The InputString class

#why does _InputStringBase not inherit from string?
#Because I didn't know of super then.
class _InputStringBase(object):
    """class for a string with current position and extras """
    def __init__(self,string):
        """The constructor

        @param string The string this object encapsulates
        """
        self.set_member()
        self.string=string

    #functions that change the internal state
    
    def set_start(self,offset=0):
        """Sets a mark to current position+start
        
        Later a substring can be extracted beginning at this postition
        @param offset offset of marker from current position
        """
        self.begin=self.i+offset

    def goto(self,character):
        """Jumps to character
        
        This method jumps to the next occurence of character.
        This method does not jump to a specific index.
        @param character character to which should be jumped
        @return the position to which we jumped
        """
        self.i=self.string.index(character,self.i)
        return self.i

    def reset(self,i=0,begin=0,toggle=False):
        """ Resets the counters

        For testing
        @param i position to which the counter should be set
        @param begin position to which the beginning of current substring should be set
        @param toggle 
        """
        self.begin=begin
        self.i=i
        self.toggle_advance_on_minus=toggle


    def find_first(self,string_list):
        """ Finds which string from input_set occures first in string
 
        @param string_list list of strings that are checked
        @return tuple (which,where):which=index of first occurence of string from delim_list on delim_list
                where=index of first occurence of str from string_list in string
        """
        index=float("inf")
        string_no=-1
        i=0
        while i<len(string_list):
            find=self.string.find(string_list[i] )
            if 0 <= find < index :
                index=find
                string_no=i
            i+=1
        return string_no,index

    
    def extend(self,string):
        """function to append a string to the current content
        
        @param string the string that should be appended"""
        self.string+=string


    
    def test_substring_match(self,re_list):
        """test if current position within the provided regular expressions

        @param re_list list of regular expressions to be tested against
        @return a bool indicating wether any match of a regular expression beginning at the mark set by set_start reaches the current position
        """
        for reg_exp in re_list:
            match=reg_exp.match(self.string[self.begin:])
            if match and  match.end() >= (self.i-self.begin):
                return True
        return False


    #******************************************
    # Overwrite of builtins
    def __len__(self):
        """allows the use of len(string)"""
        return len(self.string)

    def __str__(self):
        """makes object printable"""
        return self.string

    #******************************************
    #display functions (as property getters)
    @property
    def display_error(self): 
        """returns the string and a vizualization of the current character
        
        will get called by the Exception handlers"""
        return "\n"+self.string+"\n"+"-"*self.i+"^"


    #for use in: while string.has_not_ended:
    @property
    def has_not_ended(self):
        """tests if end of string is reached"""
        return self.i<len(self.string)

    @property
    def substring(self):
        """ returns a substring beginning by self.begin"""
        return self.string[self.begin:self.i]

class _HiddenMinusUtil(_InputStringBase):
    """Utility class with methods that hide the minus from objects

    The true_xxx methods are also defined here
    """

    def set_member(self):
        """sets the variables of the class that store information"""
        self.string=''
        #"pointer" to current char
        self.i=0 
        #"pointer" to beginn of current string. is updated at beginn of an object
        self.begin=0  
        # special variable, that controlls the behaviour of next()
        self.toggle_advance_on_minus=False


        # next() is called when an object has processed a character
        # Minus signs must be processed by two objects (_AddRep and _OPRep)
        # therefore only every second time we should advance over a minus the 
        # next() really will advance
    def next(self):
        """advances to next char
        
        Sets pointer to next char
        """

        if self.true_this_char != "-":
            self.i +=1
        else:
            if self.toggle_advance_on_minus:
                self.i+=1
            self.toggle_advance_on_minus=not self.toggle_advance_on_minus

    #***********************************************
    #display only functions as property getters
    @property
    def this_char(self):
        """returns the current character hides minus"""
        if self.string[self.i] != "-":
            return self.string[self.i]
        else:
            #it's basically the same
            return "+"


    @property
    def last_char(self): 
        """returns the previous character """
        #if this is a minus your operator ended on(was a) minus
        return self.string[self.i-1]


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


class InputString(_HiddenMinusUtil):
    """class to handle the string that is converted to the formula"""

    def __init__(self,string):
        """The constructor"""
        self.set_member()
        self.string=string

class _InputStringAlt(InputString):
    """alternative constructor for testing"""
    def __init__(self,string,i=0,begin=0,toggle=False):
        self.set_member()
        self.string=string
        self.reset(i,begin,toggle)

#*****************************************************************************************
#classes for the string expansion

class _StringConvUtil(object):
    """Store all functions helping to convert the string to an object hierachy

    
    """
    #intended order of use of methods:
    #self.check_illegal_char( ... )
    #if self.check_close_char(...)
    #    close
    #elif "special cases"
    #    ...
    #else:
    #    self._create_next_operation(string)


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
        """Tests if the current char closes the object.
        
        @param tuple characters that close the object"""
        if string.this_char in tuple:
            return True
        else:
            return False 

    def _create_next_operation(self,string):
        """Creates and returns a lower tier object"""
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
    """Abstract class for operator containing objects

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
        #part of an old architecture. Does almost nothing now.
        raise NotImplementedError()

    def extract(self):
        #returns a sum of Operator products (a list of lists)
        raise NotImplementedError()




class _SinglePartOperation( _ExpFormBaseClass):
    """ Sets members for all objects, that contain only one other object: (),<>"""

    def set_members(self):
        self.content=[""]
        self.part=0


class _DoublePartOperation( _ExpFormBaseClass):
    """ Sets members for all objects, that contain exactly two other object:  *,[,]"""

    def set_members(self):
        self.content=["",""]
        self.part=0


class _MultiPartOperation( _ExpFormBaseClass):
    """ Sets members for all objects, that contain an arbitrary number of other objects:  +"""

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
        self.check_illegal_char(string, ("<",">","+","","|") ,"first") 
        ##\todo think about which chars are illegal at the beginning of objects 
        while string.has_not_ended:

            self.check_illegal_char(string,("<"))            
            if ( string.true_this_char == "+" and string.last_char == "^" ):
                #Ignore ^+ combinations
                string.next()
            elif ( string.true_this_char=="-" and string.test_substring_match([stf_negative_regexp,
                                                                          stf_number_regexp,
                                                                          stf_div_regexp]) ): 
                string.next()
                #ignore "-" that belong to a number -1.0d-4T1
                #                                  -^----^
                # thanks to toggle minus advance this will be done twice on the second minus
            elif self.check_close_char(string,("[",",","]","(",")","*","+",">","-","|")):
                self.content=[Vertex(string.substring)]
                return 0
            else:
                string.next()
        raise UnexpectedEndOfStringError(string.display_error)

    
    def extract(self):
        return [map(str,self.content)]
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
                self.check_illegal_char(string, ("]",")",",","*","+","|")) 
                self.part +=1
                #string.i can point neither to * nor to + so part may point to out of scope
                self.content.append(self._create_next_operation(string))

            elif self.check_close_char(string,(",","]",")",">","|")):
                return 

            elif string.this_char=="*": 
                self.content[self.part]=_MultRep(self.content[self.part],string)
      
            else:
                #implicit multipication T1+(T2+T3)T4
                # after handling the parethesis the _AddRep sees the T of T4 
                self.content[self.part]=_MultRep(self.content[self.part],string)
        raise UnexpectedEndOfStringError(string.display_error)


    def extract(self):
        #self.content contains a list of objects. upon extraction those return a list of 
        # additively joined lists of mutiplicatively joined operators
        # to extract those lists have simply to be concatenated
        return reduce(lambda x,y: x+y.extract(), self.content, [])


class _ComRep( _InitOperation,_DoublePartOperation):
    """Class to store the commutator parts"""

#    counter=0
#    @classfunction
#    def increase_counter(cls):
#        cls.counter+=1

    def process_string(self,string):                
        string.next() # String pointed to "["
        self.check_illegal_char(string, ("<",">","*","","+") ,"first")
        self.content[self.part]=self._create_next_operation(string)
        while string.i <len(string):
            self.check_illegal_char(string, ("<",")","|") ) 
            
            if string.this_char==",":
                if self.part==0: # If not, it is a [,, combination
                    self.part=1
                    string.next()
                    self.check_illegal_char(string,("<",")","",",","]","|"))
                    self.content[self.part]=self._create_next_operation(string)
                else:
                    raise TooManyCommasError("There may only be one ',' "\
                                             "in any commutator"+
                                             string.display_error)

            elif self.check_close_char(string,("]")):
                if self.part==1: # If not, it is a [] combination
                    self.close()
                    string.next()
                    return 
                else:
                    raise MissingCommaError("All commutators are of the form [...,...]"+\
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
        interm1 = _MultRepAlt(self.content[0],self.content[1])
        interm2 = _MultRepAlt(_OPRepAlt("-1"),_MultRepAlt(self.content[1],self.content[0]))
        return interm1.extract()+interm2.extract()


class _MultRep( _UInitOperation,_DoublePartOperation):
    """Class to store two Factors of a multiplication
    
    _MultRep only handles two factors at a time to allow for easy extraction
    """

    def process_string(self,string,obj=''):
        if string.this_char=="*": 
            #not necessarily true since (F+V)(T1+T2) also triggers a Multiplication
            string.next()

        self.check_illegal_char(string, (")","]","*",",","","+","|") )

        self.content[self.part]=self._create_next_operation(string)

        while string.has_not_ended:
            self.check_illegal_char(string,("<")) 

            if self.check_close_char(string,(",","]",")",">","+","|")):
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
        # (a+b+c)*(d+e+f) = sum over all sets with one of first parenthesis and one out of the second 
        return [x+y for x in self.content[0].extract() for y in self.content[1].extract()]


class _MultRepAlt(_MultRep):
    """ When you already have two objects an only need to multiply them. """
    #look at _ComRep
    def __init__(self,obj1,obj2):
        self.set_members()
        self.content=[obj1,obj2]
        self.part=1

        
class _ParRep( _InitOperation,_SinglePartOperation):
    """Class to store the content of a parenthesis"""

    def process_string(self,string):
        string.next()
        self.check_illegal_char(string, ("<",">","*","]",")","","+"),"first") 
        self.content[self.part]=self._create_next_operation(string)

        while string.has_not_ended:
            self.check_illegal_char(string,("<",",","]","|")) 

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


class BracketRep(_InitOperation,_SinglePartOperation):
 
    def process_string(self,string):
        string.next()
        self.check_illegal_char(string, ("<",",","]",")","") ) #
        self.content[self.part]=self._create_next_operation(string)
        

        
        while string.has_not_ended:
            self.check_illegal_char(string, ("<",",","]",")","") ) # 

            if string.this_char=="+": 
                self.content[self.part]=_AddRep(self.content[self.part],string)
            elif self.check_close_char(string, (">","|") ):
                self.close()   
                return 
            else: 
                self.content[self.part]=_MultRep(self.content[self.part],string)
        raise UnexpectedEndOfStringError(string.display_error)
        
    def close(self):
        return  

    def extract(self):
        return self.content[0].extract()














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
        self.test_class=BracketRep
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
        print( self.run_test.__doc__ )
        
        error_counter=0
        error_counter+=self.test_OPRep()
        error_counter+=self.test_MultRep()
        error_counter+=self.test_AddRep()
        error_counter+=self.test_ComRep()
        error_counter+=self.test_BracketRep()
        print("there were ",error_counter," Errors")
        

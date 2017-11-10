import re #regular expression support
import xml.dom.minidom as minidom
import os # to find the set_keywords file

import unittest

default_keyword_file = os.getenv("GECCO_DIR")+"/data/keyword_registry.xml"

key_root_tag = "key_root"
argument_tag = "argument"
keyword_tag = "keyword"


class ContextError(Exception):
    def __init__(self,*args):
        self.msg = ""
        Exception.__init__(self,*args)




def is_empty(container):
    """ checks if a given container contains no elements

    @param container a container (list, dict, set)
    @return True if the container is empty, False if the container is not empty
    """
    return len(container) == 0
    

     
class Converter(object):
    """class to take care of conversion between xml strings and python objects"""

    def to_bool(self,string):
        """converts an xml string to a boolean


        @param string an xmlstring either 'true' of 'false'
        @return either True or False
        @raises ValueError if the string is not convertible
        """
        if string.strip() == "true":
            return True
        elif string.strip() == "false":
            return False
        else :
            raise ValueError("Unconverible boolean")

    
    def to_float(self, string):
        """converts an string to a float

        Fortran may use d as exponent separator, this function replaces this by e
        @param string a floating point number as a string
        @return the float
        @raises a Value Error, if it was not convertible
        """
        return float( re.sub( "[dD]", "e", string.strip() ) )

    def to_int(self, string):
        """converts an string to an integer

        for completeness
        @param string an integer as a string
        @return the integer
        @raises a Value Error, if it was not convertible
        """
        return int( string.strip() )

    def to_str(self, string):
        """converts an string to an integer
    
        for completeness
        @param string string
        @return the string
        """
        return str(string)


class GeCCoValueExtractor(object):
    """ Holds some GeCCo specific conversions
    """

    
    def __init__(self, converter=Converter()):
        self._conv = converter
        self._convfunc = {8:converter.to_str,
                          2:converter.to_int,
                          1:converter.to_bool,
                          4:converter.to_float,}
        self._typedict = {8:str,
                          2:int,
                          1:bool,
                          4:float,}
        
    def _get_type_convfunc(self, argument):
        """ returns a function that can convert the value of the given argument to the appropriate python type"""
        return self._convfunc[self._getitype(argument)]

    def _get_convfunc(self, argument):
        """ generates a conversion function for elements that are arrays"""
        if self.get_len(argument) >1 and self.get_type != str:
            def conv(arg):
                return [self._get_type_convfunc(argument)(elem) for elem in arg.split(",")]
            return conv
        else:
            return self._get_type_convfunc(argument)
        
    
    def get_type(self, argument):
        """ returns the python class corresponding to the type of an argument

        @param argument an argument node
        @return one of (int, bool, str, float)
        @raises ValueError if the type can not be determined
        """
        return self._typedict[self._getitype(argument)]

    
    def get_len(self, argument):
        """extracts the length of an argument

        Note that by gecco convention the length of a string is its maximal length
        @param argument an argument node
        @return the length of that argument
        @raises Value or TypeError if a length attribute was available but unreadable
        """
        return 1 if ( argument.getAttribute("len") == "" ) else self._conv.to_int( argument.getAttribute("len") )
    
    def _get_itype(self, argument):
        """returns the integer specification of an arguments type"""
        try:
            itype = self._conv.to_int(
                argument.getAttribute("type")
            )
            self._typedict[itype] #trigger key error if not in 
        except ValueError, KeyError:
            raise ValueError("Argument with unknown type"+argument.getAttribute("type"))
        else:
            return itype

    def get_value(self, argument):
        """ returns an arguments value converted to an appropriate python class

        @param argument an argument node
        @return the converted value
        @raises ValueError if the type can not be determined
        """
        convfunc = self._get_convfunc(self, argument)
        try:
            return convfunc(argument.getAttribute("value"))
        except ValueError as err:
            err.msg += "in element "+argument.getAttribute("name")
            raise


class RegistryHandler(object):
    def __init__(self, filename=None):
        filename = filename if ( filename is not None) else default_keyword_file
        self._root = self._parse(filename).getElementsByTagName(key_root_tag).item(0)

    def _parse(self, filename):
        return minidom.parse(filename)
    
    def _find_element_core(self, context_list, parent):
        """recursive core of find element"""
        if is_empty(context_list):
            return parent
        head , context_list = context_list[0], context_list[1:]
        for item in parent.childNodes:
            if (item.nodeType == item.ELEMENT_NODE ) and item.hasAttributes()  \
               and item.getAttribute("name") == head:    # returns "" if name is undefined doesn't throw
                return self._find_element_core(context_list, item)
        else:
            raise ContextError("Could not find {name} in context: ".format(name=head) )

    def _find_element(self, context, name, root):
        """finds an element or argument node in context"""
        context_list = context.strip().split(".")
        context_list = [name] if context_list == [''] else context_list + [name]
        try:
            return self._find_element_core(context_list, root)
        except ContextError as ex:
            ex.msg += context
            raise     # reraise same exception

    def get_value_by_context(self, context, name):
        elem = self._find_element(context, name, self._root)
        return GeCCoValueExtractor().get_value(element)
    
    def get_type_by_context(self, context, name):
        elem = self._find_element(context, name, self._root)
        return GeCCoValueExtractor().get_type(element)

    def does_exist(self, context, name):
        try:
            elem = self._find_element(context, name, self._root)
        except ContextError:
            return False
        else:
            return True

    def is_argument(self, context, name):
        elem = self._find_element(context, name, self._root)
        return elem.tagName == ArgumentTag
    
    def is_keyword(self, context, name):
        elem = self._find_element(context, name, self._root)
        return elem.tagName == KeywordTag
    
    def get_type_convfunc(self, context, name):
        elem = self._find_element(context, name, self._root)
        return GeCCoValueExtractor().get_type_convfunc(element)


        

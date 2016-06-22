import argparse
import sys,os

import xml.dom.minidom as dom

import re




def mybool(string):
    if (bool(string)):
        return "true"
    else:
        return "false"

def myfloat(string):
    string=str(float(string))
    return re.sub("[e,E]","d",string)

_typ_dict={"8":"str",
           "2":"int",
           "1":"log",
           "4":"rl8",}


_typ_dict_r={"str":"8",str:"8","character":"8",
             "int":"2",int:"2","integer":"2",
             "log":"1",bool:"1","bool":"1","logical":"1",
             "rl8":"4",float:"4","float":"4","real":"4"}

class CustomException(Exception):
    pass
class NodeNotFoundError(CustomException):
    pass

class ExistingElementError(CustomException):
    pass

class InconsistentRegistryError(CustomException):
    pass
    

def _format_keyword(key):
    status=key.getAttribute("status")
    name=key.getAttribute("name")
    return "{stat} Keyword {plz}{name}"\
           .format(stat=status,
                   plz="-"*max(0,len(sys.ps1)-9),
                   name=name)

def _format_argument(arg):
    status=arg.getAttribute("status")
    name=arg.getAttribute("name")
    length=arg.getAttribute("length")
    typ=_typ_dict[arg.getAttribute("type")]
    val=arg.getAttribute("value")
    return "{stat} Argument {plz}{name},{typ}({length}):{val}"\
             .format(stat=status,
                     plz="-"*max(0,len(sys.ps1)-10),
                     name=name,
                     typ=typ,
                     length=length,
                     val=val)

def _format_Node(self):
    if self.nodeName == "keyword":
        return _format_keyword(self)
    elif self.nodeName == "argument":
        return _format_argument(self)
    else:
        return self.__repr__()




def _format_comment(self):
    return " "*len(sys.ps1)+self.data

#some monkeypatching (pfui)
dom.Element.__str__=_format_Node
dom.Comment.__str__=_format_comment

class Argument(object):
    cast_dict={"8":str,
               "2":int,
               "1":mybool,
               "4":myfloat,}
    
    
    @staticmethod
    def det_length(value,typ, length =None):
        """ returns the length of an argument with or without given length"""
        if (length is None):
            if (_typ_dict[typ] == "str" or 
                isinstance(value,list) or 
                isinstance(value,tuple) ):
                length = str(len(value))
            else: 
                length="1"
        else:
            length=str(int(length))
        return length

    @classmethod
    def cast_value(cls,value,typ):
        """converts the value into a digestible string

        only accepts string-numericals as typ
        """
        cast=cls.cast_dict[typ]
        if ( isinstance(value,list) or 
             isinstance(value,tuple) ):
            value=",".join(
                [str(cast(val)) for val in value])
        elif value is not None :
            value=str(cast(value))




    
class _Registry(object):
    """ class to represent the global state"""

    #default values
    file_=os.getenv("GECCO_DIR")+"/data/keyword_registry.xml" 
    key_root_tag="key_root"

    def __init__(self):
        self.tree=None
        self.current_keyword=None
        self.key_root_tag=self.key_root_tag #creating object instance
        self.root=None

    def get_tree(self):
        """ returns the registry as DOM-document"""
        return self.tree

    def set_root(self,root_tag=None):
        """ sets the keyword root element"""
        root_tag=root_tag if (root_tag is not None) else self.key_root_tag
        if root_tag is not None:
            self.root=self.tree.getElementsByTagName(root_tag)[0]
        else:
            self.root=self.tree._get_firstChild
        self.key_root_tag=self.root.nodeName #make sure this is consistent!()


    def load_file(self,file=None):
        """loads a (new) xml file, the default file if no file is given"""
        file = file if (file is not None) else self.file_
        self.tree=dom.parse(file)
        self.set_root()

    def saveto_file(self, file):
        file = file  if (file is not None) else self.file_
        with open(file, "w") as f:
            self.tree.writexml(f,indent="",addindent="",newl="")
    
    def get_curContext(self):
        """ returns the context of the current keyword (including itself) as list"""
        current=self.current_keyword
        contextlist=[]
        while  current.nodeName not in  [self.key_root_tag, "#document"]:
            contextlist=[current.getAttribute("name")]+contextlist
            current=current.parentNode
        return contextlist


        
    def get_curChildren(self,tag=None,name=None):
        """ Returns a list of all child elements of the current node"""
        return [ child for child in self.current_keyword.childNodes
                 if ( tag is None or tag == child.nodeName)
                 if (name is None or ( child.nodeType==1  and hasAttribute("name") and
                                       child.getAttribute("name") == name) )]

    def moveto_Root(self):
        self.current_keyword=self.root

    def moveto_Child(self,name):
        """sets moves the currently active node to a child with name name"""
        nodes=self.get_curChildren(tag="keyword", name=name)
        if len(nodes) == 0:
            raise NodeNotFoundError("No keyword found, named:{name}".format(name=name))
        elif len(nodes) ==1 :
            self.current_keyword=nodes[0]
        else:
            raise InconsistentRegistryError("More than one keyword of name {name} found.".format(name=name))


    def moveto_Parent(self):
        self.current_keyword=self.current_keyword.parentNode

    def _append_comment(self,node,comment):
        """appends a comment to the given Node"""
        if comment is not None and len( str(comment).strip() ) != 0 :
            new_comment=self.tree.createComment(str(comment).strip() )
            node.appendChild(new_comment)

    def append_argument(self,name,value,typ,length,description):
        typ=_typ_dict_r[typ]
        name=str(name)
        length=Argument.det_length(value,typ,length)
        value=Argument.cast_value(value,typ)
        if len(self.get_curChildren( name=name)) > 0:
            raise ExistingElementError("This name already exists")
        new_elem=self.tree.createElement("argument")
        new_elem.setAttribute("name",name)
        new_elem.setAttribute("length",length)
        new_elem.setAttribute("status","A")
        new_elem.setAttribute("type",typ)
        new_elem=self.current_keyword.appendChild(new_elem)
        if value is not None:
            new_elem.setAttribute("value", value)
        self._append_comment(new_elem,description)
        return new_elem

    def append_keyword(self,name,comment):
        """appends a possibly commented keyword to current element"""
        name=str(name)
        if len(self.get_curChildren(name=name)) > 0:
            raise ExistingElementError("This name already exists")
        new_elem=self.tree.createElement("keyword")
        new_elem.setAttribute("name",name)
        new_elem.setAttribute("status","A")
        new_elem=self.current_keyword.appendChild(new_elem)
        if comment is not None:
            self._append_comment(new_elem,comment)
        return new_elem

    def delnode(self, name ):
        existing=self.get_curChildren(name=name)
        if len(existing) == 0:
            raise ExistingElementError("Element doesn't exist")
        elif (len(existing) > 1  ):
            raise InconsistentRegistryError("More than one element of name {name} exists.".format(name=name))
        else:
            self.current_keyword.removeChild(existing[0]).unlink()

_registry=_Registry()

_parser = argparse.ArgumentParser(description='Edit the gecco keyword registry')

_parser.add_argument('--file', metavar='file', type=str, default=_registry.file_,
                     dest='file',help='the xml file containing the keyword_registry')

def _get_comment(node):
    for subnode in node.childNodes:
        if subnode.nodeName=="#comment":
            return subnode
    return None




def _update_promt():
    _set_promt(_registry.get_curContext())

def _set_promt(context_list):
    sys.ps1=".".join(context_list)+">" if (len(context_list) >0) else ">" 
    sys.ps2="."*(len(sys.ps1)-0)+">"


def load(file):
    """load(file) load the xml file file"""
    _registry.load_file(file)
    _registry.set_root()
    _registry.moveto_Root()
    _update_promt()
    
def save(file=None):
    """ save(file) save the registry to the given file

    if no file is given: save to the file last loaded from
    """
    _registry.saveto_file(file)
    
def cd(context):
    """cd(context) change current keyword

    context has to be one of 
    ".." -> going one level up
    "~"  -> going back to the root keyword
    a string of type <key>[.<key>[.<key> ...]]
    """
    try:
        if context in ["..","../"] :
            _registry.moveto_Parent()
        elif context in ["~","~/"] :
            _registry.moveto_Root()
        else :
            for key in context.split("."):
                _registry.moveto_Child(key)
        _update_promt()
    except NodeNotFoundError:
        print "That keyword doesn't exist'"
    except InconsistentRegistryError:
        print "something is wrong with this registry"

def ck(context):
    """ck(context) change current keyword

    context has to be one of 
    ".." -> going one level up
    "~"  -> going back to the root keyword
    a string of type <key>[.<key>[.<key> ...]]
    """
    cd(context)


def ls(flags=""):
    """lists the keywords and arguments of the next sublevel

    default format is <status> <kind> <placeholder> <name>[,<type>(<length>):<value>]\n[<comment>] 
    """
    for node in _registry.get_curChildren(tag="keyword"):
        print _format_keyword(node)
        comment=_get_comment(node)
        if comment is not None:
            print _format_comment(comment)
    for node in _registry.get_curChildren(tag="argument"):
        print _format_argument(node)
        comment=_get_comment(node)
        if comment is not None:
            print _format_comment(comment)
            
def get_tree():
    """ returns the currently loaded registry as a dom.minidom object for direct manipulation
    """
    return _registry.get_tree()
                
def main(parser=_parser,args=None):
    "function to set inital values not for public use"
    args=parser.parse_args(args)
    load(args.file)

def mkkey(name,description=None):
    """ mkkey(name,[description])

    append a new keyword as subkeyword to the current keyword
    """
    try:
        elem=_registry.append_keyword(name,description)
    except InconsistentRegistryError as ex:
        print "something went really wrong:"
        print ex
    except ExistingElementError:
        print "Sorry that name already exists."
    else:
        print "New keyword included:"
        print _format_keyword(elem)

def mkarg(name,value,typ,length=None,description=None):
    """ mkarg(name,value,typ,[length],[description])
    
    append a new argument below the current keyword
    name: name of the new argument, may not already exist
    value: default value, may be None for No default lists for array defaults
    typ type of of argumentvalue:
        str or "str" or "character"
        int or "int" or "integer"
        float or "float" or "rl8" or "real"
        bool or "bool" or "logical" or "log"
    length maximum length of the values
        if no value given, current length is assumed maximum
    description: description of the new argument 
    """
    #bugger them a bit, so they don't forget it too often
    if description is None: 
        description=raw_input("Please give a short description for the argument:\n")
    description=None   if (len(description.strip()) == 0) else  description
    try:
        elem=_registry.append_argument(name,value,typ,length,description)
    except InconsistentRegistryError as ex:
        print "something went really wrong:"
        print ex
    except ExistingElementError:
        print "Sorry that name already exists."
    else:
        print "new argument included:"
        print _format_argument(elem)


def delkey(name):
    """delkey(name) deletes a keyword subelement to the current keyword"""
    try:
        _registry.delnode(name)
    except InconsistentRegistryError as ex:
        print ex.msg
        print "Sorry the keyword_editor cannot uniquely identify the element. "
        print "you will have to remove it by hand."
    except ExistingElementError as ex:
        print ex

def delarg(name):
    """delarg(name) deletes an argument subelement to the current keyword"""
    delkey(name)


def describe(func=None):
    """ describes a function. if no function is given, it lists all functions made available by 
    registry_editor"""
    if func is None:
        print "The following commands are currently available"
        for cmd in __all__:
            print cmd
    else :
        string=func.__doc__
        sting=string if ( string != "" ) else "Currently there is no description available for this function"
        print string

__all__=["ls",
         "load",
         "save",
         "cd",
         "ck",
         "mkkey",
         "mkarg",
         "delarg",
         "delkey",
         "get_tree",
         "describe",
]


import argparse
import unittest as ut
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

_typ_cast_dict={"8":str,
                "2":int,
                "1":mybool,
                "4":myfloat,}

_typ_dict_r={"str":"8",str:"8","character":"8",
             "int":"2",int:"2","integer":"2",
             "log":"1",bool:"1","bool":"1","logical":"1",
             "rl8":"4",float:"4","float":"4","real":"4"}

class NodeNotFoundError(Exception):
    pass

class ExistingElementError(Exception):
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
        return self.tree

    def set_root(self,root_tag=None):
        """ sets the keyword root element"""
        root_tag=root_tag if root_tag is not None else self.key_root_tag
        if root_tag is not None:
            self.root=self.tree.getElementsByTagName(root_tag)[0]
        else:
            self.root=self.tree._get_firstChild
        self.key_root_tag=self.root.nodeName #make sure this is consistent!()


    def load_file(self,file=None):
        """loads a (new) xml file, the default file if no file is given"""
        file = self.file_ if file is None else file 
        self.tree=dom.parse(file)
        self.set_root()

    def save_to_file(self, file):
        file = file  if file is not None else self.file_
        with open(file, "w") as f:
            self.tree.write(f,indent="",addindent="",newl="")
    
    def get_currentcontext(self):
        """ returns the context of the current keyword (including itself) as list"""
        current=self.current_keyword
        contextlist=[]
        while  current.tagName not in  [self.key_root_tag, "#document"]:
            contextlist=[current.getAttribute("name")]+contextlist
            current=current.parentNode
        return contextlist

    def get_currentChildren(self):
        """ Returns a list of all child elements of the current node"""
        return self.current_keyword.childNodes

    def moveto_Root(self):
        self.current_keyword=self.root

    def moveto_Child(self,name):
        """sets moves the currently active node to a child with name name"""
        for node in self.current_keyword.childNodes:
            if (node.getAttribute("name") == name):
                if (node.nodeName == "keyword"):
                    self.current_keyword=node
                    return None
        raise NodeNotFoundError("No keyword found, named:{name}".format(name=name))

    def moveto_Parent(self):
        self.current_keyword=self.current_keyword.parentNode

    def _append_comment_to(self,node,comment):
        comment=str(comment)
        if comment is not None:
            new_comment=self.tree.createComment(comment)
            node.appendChild(new_comment)

    def append_argument(self,name,value,typ,length,description):
        typ=_typ_dict_r[typ]
        name=str(name)
        if length is None and \
           (_typ_dict[typ] == "str" or \
            isinstance(value,list) or \
            isinstance(value,tuple) ):
            length = str(len(value))
        elif length is None:
            length="1"
        else:
            length=str(int(length)) # make sure, length can be converted to integer
        cast=_typ_cast_dict[typ]
        if isinstance(value,list) or \
           isinstance(value,tuple) :
            value=",".join(
                [str(cast(val)) for val in value])
        elif value is not None :
            value=str(cast(value))
        existing=[key for key in self.get_currentChildren() if key.getAttribute("name") == name ]
        if len(existing) > 0:
            raise ExistingElementError("This name already exists")
        new_elem=self.tree.createElement("argument")
        new_elem.setAttribute("name",name)
        new_elem.setAttribute("length",length)
        new_elem.setAttribute("status","A")
        new_elem.setAttribute("type",typ)
        new_elem=self.current_keyword.appendChild(new_elem)
        if value is not None:
            new_elem.setAttribute("value", value)
        self._append_comment_to(new_elem,description)
        return new_elem

    def append_keyword(self,name,comment):
        """appends a possibly commented keyword to current element"""
        name=str(name)
        existing=[key for key in self.get_currentChildren() if key.getAttribute("name") == name ]
        if len(existing) > 0:
            raise ExistingElementError("This name already exists")
        new_elem=self.tree.createElement("keyword")
        new_elem.setAttribute("name",name)
        new_elem.setAttribute("status","A")
        new_elem=self.current_keyword.appendChild(new_elem)
        if comment is not None:
            self._append_comment_to(new_elem,comment)
        return new_elem

    def delnode(self, name):
        existing=[key for key in self.get_currentChildren() if key.getAttribute("name") == name ]
        if len(existing) == 0:
            raise ExistingElementError("Element doesn't exist")
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
    _set_promt(_registry.get_currentcontext())

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
    _registry.save_file(file)
    
def cd(context):
    """cd(context) change current keyword

    context has to be one of 
    ".." -> going one level up
    "~"  -> going back to the root keyword
    a string of type <key>[.<key>[.<key> ...]]
    """
    if context in ["..","../"] :
        _registry.moveto_Parent()
    elif context in ["~","~/"] :
        _registry.root()
    else :
        for key in context.split("."):
            _registry.moveto_Child(key)
    _update_promt()

ck=cd


def ls(flags=""):
    """lists the keywords and arguments of the next sublevel

    default format is <status> <kind> <placeholder> <name>[,<type>(<length>):<value>]\n[<comment>] 
    """ 
    for node in _registry.get_currentChildren():
        if node.nodeName == "keyword":
            print _format_keyword(node)
            comment=_get_comment(node)
            if comment is not None:
                print _format_comment(comment)
    for node in _registry.get_currentChildren():
        if node.nodeName == "argument":
            print _format_argument(node)
            comment=_get_comment(node)
            if comment is not None:
                print _format_comment(comment)

def get_tree():
    return _registry.get_tree()
                
_main_run=True
def main(parser=_parser,args=None):
    "function to set inital values not for public use"
    global _main_run
    if _main_run:
        args=parser.parse_args(args)
        load(args.file)
    _main_run=False

def mkkey(name,description=None):
    """ mkkey(name,[description])

    append a new keyword as subkeyword to the current keyword
    """
#    if description is None:
#        description=raw_input("Please give a short description for the use of the keyword")    
    elem=_registry.append_keyword(name,description)
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

    elem=_registry.append_argument(name,value,typ,length,description)
    print _format_argument(elem)

def delarg(name):
    """delarg(name) deletes an argument subelement to the current keyword"""
    _registry.delnode(name)
def delkey(name):
    _registry.delnode(name)

__all__=["ls",
         "main",
         "load",
         "cd",
         "ck",
         "mkkey",
         "mkarg",
         "delarg",
         "delkey",
         "get_tree"]

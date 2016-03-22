import re #regular expression support
import xml.etree.ElementTree as ET #xml and DOM functionality
import os # to find the set_keywords file




keyword_file=os.getenv("GECCO_DIR")+"/input/set_keywords.f"


output_file='./config_xml' # not used here , module originated from a trial to create a config file

regexp_keyword =re.compile(r"^ +call +keyword_add")
regexp_keyword_extr=re.compile(r"^ +call +keyword_add\('(?P<keyword>.+?)'.+?(?:context='(?P<context>.+?)')?")
regexp_argument = re.compile(r"^ +call +argument_add")
regexp_argument_extr = re.compile("^ +call +argument_add *\( *"+
                                  "'(?P<argument>.+?)' *, *"+
                                  "(?:context=)?'(?P<context>.+?)' *, *"+
                                  "type=vtyp_(?P<type>.+?) *, *"+
                                  "(?:len=(?P<len>[0-9]+?)(?: *, *|\))"+")?"+
                                  "(?:"+".def=\(/(?P<value>.+?)/\)"+")?"+
                                  "(?:\)? *! *(?P<comment>.*))?")
#Description of regexp:
# "^" matches the beginning and " +" matches an arbitrary amount (at least on) whitespace.
    #together these make sure that commented lines are not matched.
#Then comes "call argument_add("
# " *" matches zero or more whitespaces 
# Then comes the first argument (?P<something> names the match for easier retrieval
# ".+?" matches one or more arbitrary characters the ? makes sure the match is non greedy 
     #and goes only to the following "'" 
#The first match is the first argument of the function call
     # a string in single quotes ''
#the next terms should be clear.
# the construction (?: )? means zero or one of the matches in the the parenthesis

#matches an empty line, 
# important, does (maybe) not match an empty string
regexp_empty=re.compile("^ *\n")

regexp_comment_line=re.compile("^C|^ +!",flags=re.I)
#matches commented out lines

regexp_continuation=re.compile(r" {5}& *")





#dictionary to transform the types back
#currently not used
_type_dict={"str":("vtyp_str","cdef"),
            "int":("vtyp_int","idef"),
            "log":("vtyp_log","ldef"),
            "rl8":("vtyp_rl8","xdef"),}








class UnknownContextError(Exception):
    def __init__(self,something):
        self.context=something
    



_g_arguments=[]
_g_keywords=[]
#yeah, global variables are bad

root=ET.Element("config", attrib={"name":"root"})
tree=ET.ElementTree(element=root)


def _read_firstline(f,line):
    """ finds the next line that is the beginning of a desired subroutine call"""
    while not (regexp_keyword.match(line) or regexp_argument.match(line) or line=="") :
#dbg
#        if not (regexp_empty.match(line) or regexp_comment_line.match(line)):
#            # print all discarded lines that are not empty or commented out
#            #to make sure nothing important gets dicarded.
#            print line
#dbg
        line=f.readline()
    return line

def _find_continuation(f,line):
    """finds continuation lines and appends them to line"""
    next_line=f.readline()    
    while (regexp_continuation.match(next_line) or regexp_comment_line.match(next_line)) :
        if regexp_continuation.match(next_line):
            line=line+regexp_continuation.sub("",next_line)
        next_line=f.readline()
    return re.sub("\n","",line),next_line

def _extract_parameters(line):
    """extracts keywords and arguments to the global arrays"""
    match_key=regexp_keyword_extr.match(line)
    match_arg=regexp_argument_extr.match(line)
    
    if match_arg:
        _g_arguments.append(match_arg)
#dbg
#            print match_arg.group("argument")
#            print match_arg.group("type") 
#            print match_arg.group("len")
#            print "contest",match_arg.group("context")
#dbg
    elif match_key:
        _g_keywords.append(match_key)
#dbg
#            print match_key.group("keyword")
#            print "context", match_key.group("context")
#dbg
    else:
        raise Exception("line could not be extracted:\n"+line)

def _read_keyword_file(keyword_file):
    with open(keyword_file,'r') as f:
        line="\n"
        next_line="\n"
        i=0
        while line != "":
            i+=1
            # Readline returns "" after EOF
            line=_read_firstline(f,line)
            if line=="":
                break
            line,next_line=_find_continuation(f,line)
            match=_extract_parameters(line)
            line=next_line













def _find_parent_core(context,root):
    """gets the parent Element for a  keyword or argument from context
    """
    for item in root.getchildren():
        #iterates through all children of root.
        regexp_context_extr=re.compile("^"+item.get("name")+"(?:\.|$)"+"(?P<rest>.+)?")
        match=regexp_context_extr.match(context)
        #Looks if context is of the form "item_name.<rest>"
        ##todo It may be faster and more straigthforward to test by string slicing 
        #and extract in the end with only one regexp 
        if match and match.group("rest"):
            #there was a match and 
            return _find_parent_core(match.group("rest"),item)
        elif match:
            #there was a match but all context is resolved
            return item

    raise UnknownContextError(context) 


def _find_parent  (context,root):   
    """wrapper for _find_parent_core
    """
    try:
        return _find_parent_core(context,root)
    except UnknownContextError as ex :
        print "default_keywords.py"
        print "no context found for: " ,match.group("keyword")
        print match.group("context")+"\n at level:"+ex.context
#        raise UnknownContextError(match.group("context")+"\n at level:"+ex.context)



def _find_parent_wrap(match,root):
    """A wrapper for the wrapper so I can reuse the internal wrapper

    """
#dbg
#    try:
#        print "name of keyword",match.group("keyword")
#    except:
#        print "name of argument",match.group("argument")
#dbg
    if match.group("context")==None:
        return root  
    return _find_parent(match.group("context"),root)






def _find_tag(match_key):
    """returns a fitting tag for any keyword

    this implements my naming conventions
    """
    #top level keywords are general contexts
    if match_key.group("context")==None:
        return "context"
    #subkeywords of method are methods
    elif match_key.group("context") =="method" and match_key.group("keyword")!="truncate" :
        return "method"
    else:
        return "keyword"

def _append_keyword(match_key,root):
    """ appends a keyword to the xml tree as an element"""
    new_elem=ET.Element(_find_tag(match_key))
    new_elem.set("name",str(match_key.group("keyword")))
    if (new_elem.tag=="method"):
        new_elem.set("file","none")
        #file that defines the method(relative to GECCO_DIR)
        new_elem.set("lang","none")
        #language of the file (python or fortran)
    assert match_key!=None
    parent=_find_parent_wrap(match_key,root)
    parent.append(new_elem)
    

def _append_argument(match_arg,root):
    """ appends an argument to the xml tree as an element"""
    new_elem=ET.Element("argument")
    new_elem.set("name",str(match_arg.group("argument")))
    new_elem.set("type",str(match_arg.group("type")))    
    new_elem.set("length",str(max(match_arg.group("len"),1)))
    new_elem.text=str(match_arg.group("value"))
#    new_elem.set("name",match_key.group("argument"))
    parent=_find_parent_wrap(match_arg,root)
    parent.append(new_elem)
    if match_arg.group("comment"):
        new_elem.append(ET.Comment(str(match_arg.group("comment"))))
        

def _construct_xml_tree(tree):
    root=tree.getroot()
    for match_item in _g_keywords:
#dbg
#        for item in root.iter():
#            print item.get("name")
#dbg
        assert match_item !=None
        _append_keyword(match_item,root)
    for match_item in _g_arguments:
        _append_argument(match_item,root)
    return tree




def _check_for_close_tag(content,index):
    """ checks wether the next tag is a closing tag"""
    if content.find(">",index)>0 and content.find("/",index)<0:
        raise UnexpectedEndOfFileError("missing '/'")
 #   elif content.find(">",index)<0:
 #       raise NoClosing("missing '>'")
        #not always an error
    elif content.find("/",index)<content.find(">",index):
        return True
    return False

def _check_for_comment(content,index):
    if content.find(">",index)>0 and 0<content.find("--",index)<content.find(">",index):
        return True
    return False

#simplified version of regexp for tag matching
regexp_empty_tag=re.compile("<[^<>]+/>")
regexp_open_tag=re.compile("<[^/][^<>]*[^/]>")
regexp_close_tag=re.compile("</[^<>]*[^/]>")
regext_comment_tag=re.compile("<!-- .* -->")

regexp_parser_instruction=re.compile("<\?[^x][^m][^l][^<>]*>" ,flags=re.I)
#for parser instructions
regexp_xml_def=re.compile("<\?xml[^<>]*>")
regexp_cdata_tag=re.compile("<![CDATA[.*?]]>")
re_list=[regexp_empty_tag,regexp_open_tag,regexp_close_tag,
         regext_comment_tag,regexp_cdata_tag,regexp_parser_instruction,
         regexp_xml_def,]
class EndOfFile(Exception):
    pass
class IllformedXml(Exception):
    pass

def regexp_match(content,index):
    index=find(content,"<",index)
    if index<0:
        #communicating to outer layers via exception: not nice, but useful.
        raise EndOfFile
    i=0
    while i<len(re_list):
        regexp=re_list[i]
        match=regexp.match(content,index)
        if match:
            return match,i
        i+=1
    raise IllformedXml

def restructure_xml(file):
    content=""
    with open(file,"r") as f:
        content=f.read()
    with open(file,"w") as f:
        ident=4
        ident_lvl=0
        index=0
        start_index=content.find("<",index)
        index=start_index
        end_index=start_index 
        match=None #match for current tag,next tag, after the next tag
        match_kind=-1# kind of match corresponds to the index in re_list
        i=0
        while True:
            match,match_kind=regexp_match(content,index)
            if match_kind>=5 :# parser instructions and xml_def
                write_file(content,index,match.start(),ident=ident*ident_lvl) 
                write_file(content,match.start(),match.end(),ident=0)
                start_index=match.end()+1
                index=start_index
                ident=ident+next_ident
                next_ident=0
            elif match_kind>=3: #cdata and comment
                index=match.end()
            elif match_kind==2:
                next_ident-=1

            

        while 0<=index<len(content):
            if _check_for_close_tag(content,index):
                ident_lvl-=1
                #this is an end tag reduce the identation lvl
            elif _check_for_comment(content,index):
                pass
            elif _check_for_close_tag(content,content.find(">",index)+1):
                index= content.find(">",index)+1
            else:
                ident_lvl+=1
            index= content.find(">",index)+1
            f.write("\n"+" "*(ident_lvl*ident)+content[start_index:index])
            start_index=index  



_read_keyword_file(keyword_file)#writes to the global arrays
default_keys=_construct_xml_tree(tree)

def find_by_context(context,name,tree=default_keys):
    full_context=context+"."+name
    root=tree.getroot()

    # finds actually the designated element
    return _find_parent(full_context,root)

def convert_to_bool(string):
    """takes a string of Fortran type bools and converts it to python booleans"""
    if re.match(".true.",string,flags=re.I):
        return True
    elif re.match(".false.",string,flags=re.I):
        return False
    else :
        print "convert_to_bool:"
        print "unconvertible string: "+string
        exit
    
def convert_to_float(string):
    """Fortran uses d as exponent separator, this function replaces this by e"""
    try:
        return float(re.sub("d","e",string,flags=re.I))
    except:
        print "convert_to_float:"
        print "unconvertible string: "+string
        exit

def get_value_by_context(context,name,tree=default_keys):
    elem= find_by_context(context,name,default_keys)
    type_=elem.get("type","str")
    if    type_=="str":
        return str(elem.text)
        #by xml definition elem.text is saved as string
        #but better make this stable
    elif  type_=="int":
        return int(elem.text)
    elif  type_=="log":
        return convert_to_bool(elem.text)
    elif  type_=="rl8":
        return convert_to_float(elem.text)
    else:
        raise Exception("Unknown type: "+type_)






#tree.write(output_file)#writes an xml file
#restructure_xml(output_file)#makes the xml file readable


        

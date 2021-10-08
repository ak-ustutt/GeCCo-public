#
##
###

from os.path import basename
import re
import unittest as ut

__all__ =["export_python_target_graph"]

def export_python_target_graph( filename="python_target_graph.dot", arguments_re=None, target_list=None ):
     if target_list is None: # not PEP8 compliant but totally worth it.
          from python_interface.gecco_interface import _target_list as target_list
     matcher = _TargetMatcher(arguments_re_string=arguments_re)
     nodefactory = _NodeFactory(matcher)
     subgraphs = _create_subgraphs(target_list, nodefactory)
     edges =  _create_edges(target_list, nodefactory)
     external_nodes = nodefactory.implied_nodes
     with open_file(filename,'w') as f:
         _write_graph(f, subgraphs, edges, matcher, external_nodes)

class _TargetMatcher(object):
     """checks for certain pattern in the target 

     can check if either of name, rules, or arguments of rules match a specific pattern (regular expression)
     """
     def __init__(self, name_re_string=None, rules_re_string=None, arguments_re_string=None, default = False):
          """
          @param name_re_string an regular expression to match the name (may be None):
          @param rules_re_string  an regular expression to match any rulename (may be None):
          @param arguments_re_string an regular expression to match any argument of any rule ot the target (may be None):
          @param default return value if all patterns are None (default False):
          """
          self._name_re = None if (name_re_string is None) else re.compile(name_re_string)
          self._rules_re = None if (rules_re_string is None) else re.compile(rules_re_string)
          self._arguments_re = None if (arguments_re_string is None) else re.compile(arguments_re_string)
          self._default = default



     def _is_rules_match(self, target):
          """Checks if any rulename in the target matches """
          for rule in target.rules:
               if self._rules_re.search(str(rule[0] ) ) is not None :
                    return True
          return False

     def _is_name_match(self, target):
          if self._name_re.search(str(target.name) ) is not None :
               return True
          return False

     
     def _is_arguments_match(self, target):
          for rule in target.rules:
               for name in list(rule[1].values()):
                    if isinstance( name,str):
                         if self._arguments_re.search( name  ) is not None :
                              return True
                    else:
                         try:
                              for elem in name:
                                   if self._arguments_re.search(elem ):
                                        return True
                         except TypeError:  #not iterable or not an iterable of  strings/buffers
                              pass
          return False

     def is_match(self, target):
          if self._name_re is None and self._rules_re is None and self._arguments_re is None: 
               return self._default

          if self._name_re is not None and self._is_name_match( target ):
               return True
          elif self._rules_re is not None and self._is_rules_match( target ):
               return True
          elif self._arguments_re is not None and self._is_arguments_match( target ):
               return True
          else:
               return False

     def __str__(self):
          return ", ".join(
               [regex.pattern for regex in
                [ self._name_re, self._rules_re, self._arguments_re ] if regex is not None ]
          )


class _ClusterSubgraph(object):
    num = 0
    def __init__(self, attrs={}):
        self._id = "cluster{num}".format(num=self._get_num())
        self._attrs = attrs
        self._nodes = []

    def _get_num(self):
        self.__class__.num += 1
        return self.__class__.num

    def add_node(self, node):
         self._nodes.append(node)
    def __str__(self):
       id_ = self._id
       attrs = ""+";\n".join(['{attr}="{value}"'.format(attr=attr,value=value) for attr,value in list(self._attrs.items()) ])
       node_list = ""+";\n".join( [str(node) for node in self._nodes] )
       return 'subgraph {id_}{{\n'\
              '{attrs};\n'\
              '{node_list};\n'\
              '}}'.format(id_=self._id,
                          attrs = attrs,
                          node_list = node_list )  



class _TargetNode(object):
     """ Represents a target (Node in the dependency graph)"""
     num = 0
     
     def __init__(self, target, matcher):
          self.id_ = "Target{num}".format(num=self._get_num())
          self.target = target
          self._match = matcher.is_match(target)
          
     def _get_num(self):
          self.__class__.num += 1
          return self.__class__.num
     
     def __str__(self):
          colorattr = "" if ( not self._match ) else 'color="red"'
          return '{id_} [label="{label}" {color}]'.format(id_=self.id_, label=self.target.name, color=colorattr)
    
class _ImpliedNode(object):
     """Represents a Target not set in Python but referenced as a dependency """
     num = 0
     def __init__(self, name):
          self.id_ = "ExternalTarget{num}".format(num=self.get_num())
          self.name = name

     def get_num(self):
          self.__class__.num += 1
          return self.__class__.num

     def __str__(self):
          return '{id_} [label="{label}"]'.format(id_=self.id_, label=self.name)


class _DependencyEdge(object):
     """ Represents an Edge in the graph that is a dependency"""
     def __init__ (self, parent, dependant):
          """ 
          takes two nodes, either _TargetNodes or _ImpliedNodes"""
          self._parent = parent
          self._dependant = dependant

     def __str__(self):
          return "{parent} -> {child}".format(parent=self._parent.id_, child=self._dependant.id_)


class _NodeFactory(object):
     def __init__(self, matcher):
          self._nodes = {}
          self.implied_nodes = {}
          self._matcher = matcher
          self._f_delimiters_re = re.compile("[\w]+") #regular expression for finding all delimiters that end a string in Fortran, only whitespaces atm.
          
     ## A function trying to reproduce how fortran reads the the name
     #
     def _convert_name(self, name):
          return self._f_delimiters_re.findall(name)[0]
      
     def create_node(self, target):
          node = _TargetNode(target, self._matcher)
          name = self._convert_name(target.name)
          self._nodes[name] = node
          return node

     def get_node(self, nodename):
          nodename = self._convert_name(nodename)
          try:
               return self._nodes[nodename]
          except KeyError:
               try:
                    return self.implied_nodes[nodename]
               except KeyError:
                    self.implied_nodes[nodename] = _ImpliedNode(nodename)
                    return self.implied_nodes[nodename]

def _write_graph(f, subgraphs, edges, matcher, external_nodes):
    f.write("digraph dependency_graph{\n")
    if str(matcher) != "" :
        f.write('''label="colored for '{matcher}'"\n'''.format(matcher=str(matcher)))
    f.write("\n".join( [str(subgraph) for subgraph in list(subgraphs.values())]))
    f.write(";\n".join([str(edge) for edge in edges]) )
    f.write(";\n")
    f.write(";\n".join([str(node) for node in list(external_nodes.values()) ]) )
    f.write("\n")
    f.write("}")

def _create_subgraphs(target_list,nodefactory):
     subgraphs = {}
     for target in target_list:
         if target.filename not in subgraphs:
            subgraph = _ClusterSubgraph(attrs={"label":basename(target.filename)})
            subgraphs[target.filename] = subgraph
         else :
            subgraph = subgraphs[target.filename]
         subgraph.add_node( nodefactory.create_node(target) )
     return subgraphs

def _create_edges(target_list, nodefactory):
     edges = []
     for target in target_list:
         edges += [ _DependencyEdge( nodefactory.get_node(target.name), nodefactory.get_node(dependency_name  ) )  for dependency_name in target.dependencies ]
     return edges




def open_file(*args, **kwargs):
     return open( *args, **kwargs)




#---------------------------------------------------------------------------
## some tests.
# Subgraphs and NodeFactory are currently not tested.
# tests are not in the test suite because this is a non essential addon

class _MockRule(list):
     def __init__(self, name, arguments={}):
          list.__init__(self,[name,arguments])

                  
class _MockTarget(object):
     def __init__(self, name,  filename, rules=[], dependencies=[], required=False,  joined=[]):
        self.name = name
        self.required = required
        self.dependencies = dependencies
        self.joined = joined
        self.rules = rules
        self.filename = filename



class Test_TargetMatcher(ut.TestCase):
     def setUp(self):
          rules = []
          rules.append(_MockRule("EVALUATE",{"FORM":"FORM_ME_1234"}))
          rules.append(_MockRule("REPLACE", {"OPERATORS":["OP_A","OP_B"],"REPLACEMENT":["OP_C","OP_D"]} ) )
          rules.append(_MockRule("ABORT", {}) )
          self.target = _MockTarget("DEF_FORM_1234", "xyztest.py", rules=rules, dependencies=["FORM_1243"] )
          self.empty_target = _MockTarget("CollectOperators", "zzztest.py", dependencies=["FORM_1234","FORM_1423"])
                       
     def test_nomatch(self):
          matcher = _TargetMatcher( ) #testing default settings
          for target in [self.target, self.empty_target]:
               self.assertFalse(matcher.is_match(target) )
               
     def test_nomatch2(self):
          matcher = _TargetMatcher( name_re_string="DEF_FORM_4321",
                                    rules_re_string="RELAPSE",
                                    arguments_re_string="form_me_1234", #no match, case mismatch
                                    default=False ) #testing default settings
          for target in [self.target, self.empty_target]:
               self.assertFalse(matcher.is_match(target) )
               

     def test_noargumentmatch(self):
          matcher = _TargetMatcher( name_re_string=None,
                                    rules_re_string=None,
                                    arguments_re_string="OPERATORS",
                                    default=False )
          self.assertFalse(matcher.is_match(self.target))

     def test_argumentmatch(self):
          matcher = _TargetMatcher( name_re_string=None,
                                    rules_re_string=None,
                                    arguments_re_string="FORM_ME",
                                    default=False )
          self.assertTrue(matcher.is_match(self.target) )

     def test_rulematch(self):
          matcher = _TargetMatcher( name_re_string=None,
                                    rules_re_string="ABORT",
                                    arguments_re_string=None,
                                    default=False )
          self.assertTrue(matcher.is_match(self.target) )

     def test_name(self):
          matcher = _TargetMatcher( name_re_string="FORM_1234",
                                    rules_re_string=None,
                                    arguments_re_string=None,
                                    default=False )
          self.assertTrue(matcher.is_match(self.target) )

     def test_printable(self):
          matcher = _TargetMatcher( name_re_string="DEF_FORM_4321",
                                    rules_re_string="RELAPSE",
                                    arguments_re_string="form_me_1234", #no match, case mismatch
                                    default=False )
          self.assertEqual(str(matcher), "DEF_FORM_4321, RELAPSE, form_me_1234")

class Test_TargetNode(ut.TestCase):
     def setUp(self):
          rules = []
          rules.append(_MockRule("EVALUATE",{"FORM":"FORM_ME_1234"}))
          rules.append(_MockRule("REPLACE", {"OPERATORS":["OP_A","OP_B"],"REPLACEMENT":["OP_C","OP_D"]} ) )
          rules.append(_MockRule("ABORT", {}) )
          self.target = _MockTarget("DEF_FORM_1234", "xyztest.py", rules=rules, dependencies=["FORM_1243"] )
          self.empty_target = _MockTarget("CollectOperators", "zzztest.py", dependencies=["FORM_1234","FORM_1423"])
          _TargetNode.num = 0
     def tearDown(self):
          _TargetNode.num = 0 # resetting id generation

     def test_unique(self):
          node1  = _TargetNode(self.target,_TargetMatcher() )
          node2  = _TargetNode(self.target,_TargetMatcher() )
          self.assertNotEqual(node1.id_ ,node2.id_)
     

     def test_print_nomatch(self):
          node = _TargetNode(self.target,_TargetMatcher() )
          self.assertEqual(str(node), 'Target1 [label="DEF_FORM_1234" ]' )
     
     def test_print_match(self):
          node = _TargetNode(self.target,_TargetMatcher(default = True ) ) 
          self.assertEqual(str(node), 'Target1 [label="DEF_FORM_1234" color="red"]' )

class Test_ImpliedNode(ut.TestCase):
     def setUp(self):
          _ImpliedNode.num = 0
     def tearDown(self):
          _ImpliedNode.num = 0 # resetting id generation

     def test_unique(self):
          node1  = _ImpliedNode("aqbc" )
          node2  = _ImpliedNode("aqbc")
          self.assertNotEqual(node1.id_ ,node2.id_)

     def test_print(self):
          node1  = _ImpliedNode("aqbc" )
          self.assertEqual(str(node1), 'ExternalTarget1 [label="aqbc"]' )


class Test_DependencyEdge(ut.TestCase):
     def setUp(self):
          _ImpliedNode.num = 0
          _TargetNode.num = 0
     def tearDown(self):
          _ImpliedNode.num = 0
          _TargetNode.num = 0
     def test_print(self):
          node1  = _TargetNode( _MockTarget("CollectOperators",
                                            "zzztest.py",
                                            dependencies=["FORM_1234","FORM_1423"]),
                                _TargetMatcher()
          ) 
          node2  = _ImpliedNode("aqbc" )
          edge = _DependencyEdge(node1, node2) 
          self.assertEqual(str(edge) , "Target1 -> ExternalTarget1")
     
if __name__ == "__main__":
     ut.main()

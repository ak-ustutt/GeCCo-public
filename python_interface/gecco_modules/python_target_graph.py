#
##
###

from python_interface.gecco_interface import _target_list
from os.path import basename
import re

class _TargetMatcher(object):
     def __init__(self, name_re_string=None, rules_re_string=None, arguments_re_string=None, default = False):
          self._name_re = None if (name_re_string is None) else re.compile(name_re_string)
          self._rules_re = None if (rules_re_string is None) else re.compile(rules_re_string)
          self._arguments_re = None if (arguments_re_string is None) else re.compile(arguments_re_string)
          self._default = default


     def is_arguments_match(self, target):
          for rule in target.rules:
               for name in rule[1].values():
                    if self._arguments_re.search(str(name ) ) is not None :
                         return True
          return False
     def is_rules_match(self, target):
          for rule in target.rules:
               if self._rules_re.search(str(rule[0] ) ) is not None :
                    return True
          return False

     def is_name_match(self, target):
          if self._name_re.search(str(target.name) ) is not None :
               return True
          return False

     
     def is_arguments_match(self, target):
          for rule in target.rules:
               for name in rule[1].values():
                    if self._arguments_re.search(str(name ) ) is not None :
                         return True
          return False

     def is_match(self, target):
          if self._name_re is None and self._rules_re is None and self._arguments_re is None: 
               return self._default

          if self._name_re is not None and self.is_name_match( target ):
               return True
          elif self._rules_re is not None and self.is_rules_match( target ):
               return True
          elif self._arguments_re is not None and self.is_arguments_match( target ):
               return True
          else:
               return False

     def __str__(self):
          return " ".join([regex.pattern for regex in [self._name_re, self._rules_re, self._arguments_re ]  if regex is not None ])

class _ClusterSubgraph(object):
    num = 0
    def __init__(self, attrs={}):
        self._id = "cluster{num}".format(num=self.get_num())
        self._attrs = attrs
        self._nodes = []

    def get_num(self):
        self.__class__.num += 1
        return self.__class__.num

    def add_node(self, node):
         self._nodes.append(node)
    def __str__(self):
       id_ = self._id
       attrs = ""+";\n".join(['{attr}="{value}"'.format(attr=attr,value=value) for attr,value in self._attrs.items() ])
       node_list = ""+";\n".join( [str(node) for node in self._nodes] )
       return 'subgraph {id_}{{\n'\
              '{attrs};\n'\
              '{node_list};\n'\
              '}}'.format(id_=self._id,
                          attrs = attrs,
                          node_list = node_list )  



class _TargetNode(object):
      num = 0
      
      def __init__(self, target, matcher):
         self.id_ = "Target{num}".format(num=self.get_num())
         self.target = target
         self._match = matcher.is_match(target)
      def get_num(self):
        self.__class__.num += 1
        return self.__class__.num
      def __str__(self):
         colorattr = "" if ( not self._match ) else 'color="red"'
         return '{id_} [label="{label}" {color}]'.format(id_=self.id_, label=self.target.name, color=colorattr)

class _ImpliedNode(object):
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
     def __init__ (self, parent, dependant):
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
    f.write("\n".join( [str(subgraph) for subgraph in subgraphs.values()]))
    f.write(";\n".join([str(edge) for edge in edges]) )
    f.write(";\n")
    f.write(";\n".join([str(node) for node in external_nodes.values() ]) )
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
     for target in _target_list:
         edges += [ _DependencyEdge( nodefactory.get_node(target.name), nodefactory.get_node(dependency_name  ) )  for dependency_name in target.dependencies ]
     return edges

def export_target_graph( filename="python_target_graph.dot", arguments_re=None ):
     target_list = _target_list # Explicit local from global
     matcher = _TargetMatcher(arguments_re_string=arguments_re)
     nodefactory = _NodeFactory(matcher)
     subgraphs = _create_subgraphs(target_list, nodefactory)
     edges =  _create_edges(target_list, nodefactory)
     external_nodes = nodefactory.implied_nodes
     with open(filename,'w') as f:
         _write_graph(f, subgraphs, edges, matcher, external_nodes)

 


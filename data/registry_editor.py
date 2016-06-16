import argparse
import unittest as ut
import sys

import xml.etree.ElementTree as ET

_g_tree=None
_g_root=None
_g_current_keyword=None

_parser = argparse.ArgumentParser(description='Edit the gecco keyword registry')
_parser.add_argument('--file', metavar='N', type=str, default="./keyword_registry",
                     dest='file',help='the xml file containing the keyword_registry')
_parser.add_argument('--mode', dest='mode', default="user",
                   help='level of privileges granted in editing session')



def _load(file):
    global _g_tree
    _g_tree=ET.parse(file)

def _change_prompt(level,new_key):
    pass

def _set_promt():
    sys.ps1=

def 
    
    
def pwd():
    print

def ls():

def main(parser=_parser)

class _mock_parser()

class _TestMain(ut.TestCase):
    def test_load(self):
        

import unittest as ut
from stf_top import *
from stf_top import _OPProduct,_Bracket
from python_interface.gecco_modules.stf_pack.operators import Vertex


#NOTE: For manual testing you have to comment out the import of  gecco_interface in string_to_form2.py as gecco_interface exits if the correct environment is not found

#Running non-essential tests? (string representations ...)
non_essential=True


    
class Test_OPProduct(ut.TestCase):
    def setUp(self):
        logger.clear()
        self.logger=logger
        
    def test_creation(self):
        """ test a simple creation and set_rule functionality"""
        op_list=["C0","H0","T1"]
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
        _OPProduct([Vertex(x)for x in op_list]).set_rule(inp)
        self.assertEqual(self.logger.events[0],
                         ("EXPAND_OP_PRODUCT",exp_res)
        )
        self.assertEqual(len(self.logger.events), 1)

    def test_idx_sv_creation(self):
        """tests the creation of a vertex index list"""
        op_list=["C0","H0","T1"]
        inp={LABEL:"LAG", 
             OP_RES:"L", 
             FAC:1, 
             FAC_INV:1 }            
        exp_res= {"OP_RES": 'L',
                  'OPERATORS': ['C0', 'H0', 'T1'], 
                  'LABEL': 'LAG', 
                  'FAC_INV': 1, 
                  'FAC': 1, 
                  'IDX_SV': [1, 2, 3]}
        _OPProduct([Vertex(x)for x in op_list]).set_rule(inp)
        self.assertEqual(self.logger.events[0],
                         ("EXPAND_OP_PRODUCT",exp_res)
        )
        
    def test_factors_int(self):
        """test the result with an integer prefactor
        on two operators"""
        op_list=["2C0", "2H0"]
        inp={LABEL:"LAG", 
             OP_RES:"L", 
             FAC:1, 
             FAC_INV:1 }
        exp_res={'OP_RES': 'L', 
                 'OPERATORS': ['C0',"H0"], 
                 'LABEL': 'LAG', 
                 'FAC_INV': 1, 
                 'FAC': 4.0, 
                 'IDX_SV': [1,2]}
        
        _OPProduct([Vertex(x)for x in op_list]).set_rule(inp)
        self.assertEqual(self.logger.events[0],
                         ("EXPAND_OP_PRODUCT",exp_res)
        )

    def test_factors_float(self):
        op_list=["2.4C0", "2H0"]
        inp={LABEL:"LAG", 
             OP_RES:"L", 
             FAC:1, 
             FAC_INV:1 }
        exp_res={'OP_RES': 'L', 
                 'OPERATORS': ['C0',"H0"], 
                 'LABEL': 'LAG', 
                 'FAC_INV': 1, 
                 'FAC': 4.8, 
                 'IDX_SV': [1,2]}
        _OPProduct([Vertex(x)for x in op_list]).set_rule(inp)
        self.assertEqual(self.logger.events[0],
                         ("EXPAND_OP_PRODUCT",exp_res)
        )

    def test_factors_neg_float(self):
        """ test recogincion of an operator with an negative float prefactor"""
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
        _OPProduct([Vertex(x)for x in op_list]).set_rule(inp)
        self.assertEqual(self.logger.events[0],
                         ("EXPAND_OP_PRODUCT",exp_res)
        )
    def test_factors_neg_float2(self):
        """ test recogincion of an operator with an negative float prefactor format 2"""
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
        _OPProduct([Vertex(x)for x in op_list]).set_rule(inp)
        self.assertEqual(self.logger.events[0],
                         ("EXPAND_OP_PRODUCT",exp_res)
        )
    def test_factors_neg_float2(self):
        """tests the correct evaluation with two negative prefactors"""
        op_list=["-2C0","-2H0"]
        inp={LABEL:"LAG", 
             OP_RES:"L", 
             FAC:1, 
             FAC_INV:1 }
        exp_res={'OP_RES': 'L', 
                 'OPERATORS': ['C0',"H0"], 
                 'LABEL': 'LAG', 
                 'FAC_INV': 1, 
                 'FAC': 4.0, 
                 'IDX_SV': [1,2]}
        _OPProduct([Vertex(x)for x in op_list]).set_rule(inp)
        self.assertEqual(self.logger.events[0],
                         ("EXPAND_OP_PRODUCT",exp_res)
        )

    def test_factors_float_exp(self):
        """tests the recognicion of an floating point prefactor with an exponent"""
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
        _OPProduct([Vertex(x)for x in op_list]).set_rule(inp)
        self.assertEqual(self.logger.events[0],
                         ("EXPAND_OP_PRODUCT",exp_res)
        )
    def test_factors_fractional(self):
        """tests fractional prefactors"""
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
        _OPProduct([Vertex(x)for x in op_list]).set_rule(inp)
        self.assertEqual(self.logger.events[0],
                         ("EXPAND_OP_PRODUCT",exp_res)
        )
    def test_factors_neg_sign(self):
        """tests the case where the prefactor is simply a negative presign"""
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
        _OPProduct([Vertex(x)for x in op_list]).set_rule(inp)
        self.assertEqual(self.logger.events[0],
                         ("EXPAND_OP_PRODUCT",exp_res)
        )

    def test_factors_combined(self):
        """tests a combination of prefactors"""
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
        _OPProduct([Vertex(x)for x in op_list]).set_rule(inp)
        self.assertEqual(self.logger.events[0],
                         ("EXPAND_OP_PRODUCT",exp_res)
        )
    def test_primed_sv_list_generation(self):
        """ttests sv_list generation for two primed vertices"""
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
        _OPProduct([Vertex(x)for x in op_list]).set_rule(inp)
        self.assertEqual(self.logger.events[0],
                         ("EXPAND_OP_PRODUCT",exp_res)
        )


    def test_primed_sv_list_generation2(self):
        """ttests sv_list generation for two primed vertices of different operators"""
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
        _OPProduct([Vertex(x)for x in op_list]).set_rule(inp)
        self.assertEqual(self.logger.events[0],
                         ("EXPAND_OP_PRODUCT",exp_res)
        )
    def test_primed_sv_list_generation3(self):
        """ttests sv_list generation for two primed vertices of same operator with different prime positions"""
        
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
        _OPProduct([Vertex(x)for x in op_list]).set_rule(inp)
        self.assertEqual(self.logger.events[0],
                         ("EXPAND_OP_PRODUCT",exp_res)
        )
    def test_string_representation(self):
        """tests the string representation (including primed vertices)"""
        op_list=["GAM'","H","GA'M"]       
        i=_OPProduct([Vertex(x)for x in op_list])
        self.assertEqual(str(i),'GAM*H*GAM')

    def test_string_representation2(self):
        """tests the string representation (including primed vertices and prefactors)"""
        op_list=["1/2GAM'","H","3GA'M"]       
        i=_OPProduct([Vertex(x)for x in op_list])
        self.assertEqual(str(i),'3.0/2.0*GAM*H*GAM')



class Test_Bracket(ut.TestCase):
    def setUp(self):
        logger.clear()
        self.logger=logger
    def test_creation(self):
        """Tests the creation of a _Bracket object"""
        string="<H>"
        inp={"LABEL":"L", 
             "OP_RES":"L", 
             "NEW":True}

        exp_res={'OP_RES': 'L', 'OPERATORS': ['H'], 'LABEL': 'L', 'FAC_INV': 1, 'NEW': True, 'FAC': 1, 'IDX_SV': [1]}
        _Bracket(InputString(string)).set_rule(inp)
        self.assertEqual(
            self.logger.events[0],
            ("EXPAND_OP_PRODUCT",exp_res) )
    def test_factors_float(self):
        """testst the bracket with a floatinf point prefactor"""
        string="2.4*<H>"
        inp={"LABEL":"L", 
             "OP_RES":"L", 
             "NEW":True}
        exp_res={'OP_RES': 'L', 
                 'OPERATORS': ['H'], 
                 'LABEL': 'L', 
                 'FAC_INV': 1, 
                 'NEW': True, 
                 'FAC': 2.4, 
                 'IDX_SV': [1]
        }
        _Bracket(InputString(string)).set_rule(inp)
        self.assertEqual(
            self.logger.events[0],
            ("EXPAND_OP_PRODUCT",exp_res) )

    def test_factors_fractional(self):
        """tests fractional prefactors"""
        string="2/7<H>"
        inp={"LABEL":"L", 
             "OP_RES":"L", 
             "NEW":True}
        exp_res={'OP_RES': 'L', 
                 'OPERATORS': ['H'], 
                 'LABEL': 'L', 
                 'FAC_INV': 7.0, 
                 'NEW': True, 
                 'FAC': 2.0, 
                 'IDX_SV': [1]
        }
    def test_ignore(self):
        """tests if objects between < and | and  | and > respectively are ignored"""
        string="<C0*D|H|C0+-*[[[>"
        inp={"LABEL":"L", 
             "OP_RES":"L", 
             "NEW":True}

        exp_res={'OP_RES': 'L', 
                 'OPERATORS': ['H'], 
                 'LABEL': 'L', 
                 'FAC_INV': 1, 
                 'NEW': True, 
                 'FAC': 1, 
                 'IDX_SV': [1]
        }
    def test_string_representation(self):
        string="<C0^+*(H+1/2V)C0>"
        self.assertEqual(
            str(
                _Bracket(InputString(string))
            ),
            '<\nC0^+*H*C0\n+1.0/2.0*C0^+*V*C0\n>'
        )

    def test_Bracket_with_restriction(self):
        string="<C0^+*(H*T'*T''*T'''*T'''')*C0>"
        inp={"LABEL":"L",
             "OP_RES":"L",
             "NEW":True}
        exp_res={'OP_RES': 'L', 'OPERATORS': ['C0^+','H','T','T','T','T','C0'], 'LABEL': 'L', 'FAC_INV': 1, 'NEW': True, 'FAC': 1, 'IDX_SV': [1,2,3,4,5,6,7], 'AVOID':[3,5,4,6]}
        b_with_restr = _Bracket(InputString(string))
        b_with_restr.avoid("T'","T'''")
        b_with_restr.avoid("T''","T''''")
        b_with_restr.set_rule(inp)
        self.assertEqual(
            self.logger.events[0],
            ("EXPAND_OP_PRODUCT",exp_res) )


        
if __name__=="__main__":
    suite = ut.\
            TestLoader().\
            loadTestsFromTestCase(Test_OPProduct)
    suite.addTests(ut.\
                   TestLoader().\
                   loadTestsFromTestCase(
                       Test_Bracket))
    ut.TextTestRunner(verbosity=2).run(suite)

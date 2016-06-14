"""Utility classes and functions"""
import itertools as it
import re
from stf_regexp import stf_div_extract_regexp,stf_number_extract_regexp,stf_negative_extract_regexp



class _IDXUtil(object):
    """Collection of methods to create the idx_sv list"""
    
    def generate_idx_sv(self, OPs):
        """function to generate a idx_sv list
        """
        #generator that returns a row of increasing numbers
        OPs_e=[] #list of already encountered OPs
        idx_sv=[]
        for OP in OPs:
            for ii in xrange(len(OPs_e)):
                if OP.same(OPs_e[ii]):                    
                    idx_sv.append(ii+1)       #xrange starts with 0 idx_sv starts with 1
                    break
            else:
                idx_sv.append(len(OPs_e)+1)
                OPs_e.append(OP)
        return idx_sv


class _NumberCollectUtil(object):
    """Utility class to derive a factor from operator prefactors"""

    ##
    #@param string str to be matched
    #@param nominator
    #@param denominator
    #@return tuple (x,y,z): x new nominator, y new denominator, 
    #z nonnumber part of string.
    def _eval_numbers(self, string, nominator=1, denominator=1):
        """converts the number matched by the regexps to floats

 
        """
        residue=""
        match_div=stf_div_extract_regexp.match(string)
        match_num=stf_number_extract_regexp.match(string)
        match_neg=stf_negative_extract_regexp.match(string)
        if match_div is not None :
            nominator=nominator*float(match_div.group("nominator") )
            denominator=denominator*float(match_div.group("denominator") )
            residue=match_div.group("operator")
        elif match_num is not None :
            nominator=nominator*float(match_num.group("number"))
            residue=match_num.group("operator")
        elif match_neg is not None :
            nominator=-nominator
            residue=match_neg.group("operator")
        else:
            residue=string
        nominator,denominator=self._normalize(nominator,denominator)
        return nominator,denominator,residue

    def pop_numbers(self, OP_List, nominator, denominator):
        """Collects the numbers from operators
        
        Find all numerical prefixes,
        multiply the numbers, eliminate the numbers 
        from the elements and delete those elements that become empty
        
        """
        #indexing the list from the end to avoid messing with indices during deletion
        i=-len(OP_List) # 0
        while i<0:
            i+=1
            a,b,c=self._eval_numbers(OP_List[i])
            nominator*=a
            denominator*=b
            OP_List[i]=c
            if OP_List[i] == '' or OP_List[i] is None :
                del OP_List[i]
        return nominator,denominator
    

    def _normalize(self,nominator,denominator):
        """function called to normalize nominator and denominator

        Does nothing right now.
        """
        return nominator,denominator


def combine_dicts(first,second):
    """concatenates two dictionaries

    if both have an identical key. the first dictionary takes precedence
    """
    combined = {}
    combined.update(second)
    combined.update(first)
    return combined


def remove_whites(string):    
    """removes all (Unicode) white space characters"""
    return re.sub("\s", "",string,flags=re.U)



def ggT(a,b):
    """ calculates the largest common denominator

    uses Euclids (?) Algorithm
    """
    if b == 0 :
        return a
    else:
        return ggT(b,a%b)


    ##@param num1,num2 a tuple of the form (fac,fac_inv)
    #@return a tuple of the form (fac,fac_inv)
def fraction_multiplication(num1,num2):
    """Multiplies two fractions"""
    fac1,fac_inv1=num1
    fac2,fac_inv2=num2
    return fac1*fac2,fac_inv1*fac_inv2


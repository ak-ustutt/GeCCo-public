

##\param t_shape a string of the form: ",|,|,| ... ,|,"
#\return a string of the form   ",;,|,;,|...|,;,"
def omg_generator(t_shape):
    """ Generates the shape of OMG from the corresponding T-shape

    Generates the LABEL_DESCR of OMG from the LABEL_DESCR 
    of the corresponding T for multireference cases. 
    e.g: P,H|PV,HV=>,;P,H|,;,V;PV,H
    This function is for cases without alternatives of the form [HV].
    """
    #manages one bracket [VH] in the anihilation part ()
    o_shape=''
    t_shape=t_shape.replace('v','V')
    for block in t_shape.split('|'):
        ca_string=block.split(',')
        o_shape=o_shape+''+','+'V'*ca_string[1].count('V')+';'+\
                        ca_string[0]+','+ca_string[1].replace('V','')+'|'
    return o_shape[:-1]

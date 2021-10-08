from functools import reduce
#Exceptions for the string_to_form package



class ExpFormError(Exception):
    #base class for all exceptions special to this module
    def __init__(self,*args):
        if len(args)==0:
            self.msg=""
        elif len(args)==1:
            self.msg=str(args)
        else:
            self.msg=reduce(lambda a,y:a+y ,args,"")

class DevError(ExpFormError):
    #class for all Errors in Develoment
    pass

class NotImplementedError(DevError):
    pass

class InputError(ExpFormError):
    # For all bugs on the other side of the monitor
    pass

class IllegalCharError(InputError):
    pass

class IllegalFirstCharError(IllegalCharError):
    pass

class MissingCommaError(IllegalCharError):
    pass

class TooManyCommaError(IllegalCharError):
    pass

class UnexpectedEndOfStringError(InputError):
    pass

class BadSplitError(InputError):
    pass

class MissingArgumentError(InputError):
    pass

class MissingLabelError(MissingArgumentError):
    pass

class MissingOP_ResError(MissingArgumentError):
    pass

class MissingBracket(InputError):
    pass

class UnknownEntity(InputError):
    pass

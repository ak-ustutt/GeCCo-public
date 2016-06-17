*----------------------------------------------------------------------*
*     return number of appearences of keyword in currently
*     active keyword_history
*     the keyword is given as "context" i.e. as string including all
*     higher level keywords, e.g. key.key.key
*----------------------------------------------------------------------*
      integer function is_keyword_set(context)
*----------------------------------------------------------------------*

      use parse_input2, only : inp_key_from_context,Node
      use FoX_dom, only : getNextSibling
      implicit none
      include 'stdunit.h'

      integer,parameter::
     &     ntest=00
      character(len=*),parameter ::
     &     i_am="is_keyword_set"

      character, intent(in) ::
     &     context*(*)


      type(Node), pointer ::
     &     current


      integer ::
     &     icount

      if (ntest.ge.100) then
         call write_title(lulog,wst_dbg_subr,i_am)
      end if
      icount=-1
      current => inp_key_from_context(context,latest=.false.,
     &     keycount=icount)
      

      
      is_keyword_set = -icount-1
      if (ntest.ge.100)
     &     write (lulog,*) i_am," has found",is_keyword_set,
     &     "occurences of ",trim(context)
      return
      end

      include 'par_vtypes.h'

      interface

      subroutine keyword_add(key,context,required,status)
      implicit none
      character, intent(in) ::
     &     key*(*)
      character, intent(in), optional ::
     &     context*(*)
      logical, intent(in), optional ::
     &     required
      integer, intent(in), optional ::
     &     status
      end subroutine

      subroutine argument_add(argkey,context,type,len,
     &     idef,xdef,ldef,cdef)
      implicit none
      character, intent(in) ::
     &     argkey*(*), context*(*)
      integer, intent(in), optional ::
     &     type,len
      integer, intent(in), optional ::
     &     idef(*)
      real(8), intent(in), optional ::
     &     xdef(*)
      logical, intent(in), optional ::
     &     ldef(*)
      character, intent(in), optional ::
     &     cdef(*)
      end subroutine

      integer function is_keyword_set(context)
      implicit none
      character, intent(in) ::
     &     context*(*)
      end function

      integer function is_argument_set(context,argkey,keycount)
      implicit none
      character, intent(in) ::
     &     context*(*), argkey*(*)
      integer, intent(in), optional ::
     &     keycount
      end function

      subroutine get_argument_dimension(dim,
     &     context,argkey,keycount,argcount,type)
      implicit none
      integer, intent(out) ::
     &     dim
      character, intent(in) ::
     &     context*(*), argkey*(*)
      integer, intent(in), optional ::
     &     keycount,argcount
      integer, intent(out), optional ::
     &     type
      end subroutine

      subroutine get_argument_value(
     &     context,argkey,keycount,argcount,
     &     ival,iarr,lval,larr,xval,xarr,str)
      implicit none
      character, intent(in) ::
     &     context*(*), argkey*(*)
      integer, intent(in), optional ::
     &     keycount,argcount
      logical, intent(out), optional ::
     &     lval, larr(*)
      integer, intent(out), optional ::
     &     ival, iarr(*)
      real(8), intent(out), optional ::
     &     xval, xarr(*)
      character, intent(out), optional ::
     &     str*(*)

      end subroutine

      end interface

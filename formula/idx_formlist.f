*----------------------------------------------------------------------*
      integer function idx_formlist(formname,form_info)
*----------------------------------------------------------------------*
*     given an formula name, search form_info and return index of
*     of corresponding operator
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'
      include 'mdef_formula_info.h'

      integer, parameter ::
     &     ntest = 00

      character(*), intent(in) ::
     &     formname
      type(formula_info), intent(in) ::
     &     form_info

      integer ::
     &     iform

      if (ntest.ge.100) then
        write(luout,*) '--------------------'
        write(luout,*) 'this is idx_formlist'
        write(luout,*) '--------------------'
        write(luout,*) ' looking for: "',trim(formname),'"'
      end if

      idx_formlist = -1
      do iform = 1, form_info%nform
        if (trim(formname).eq.
     &      trim(form_info%form_arr(iform)%form%label)) then
          idx_formlist = iform
          exit
        end if
      end do

      if (ntest.ge.100) then
        write(luout,*) 'result: ',idx_formlist
      end if

      return
      end

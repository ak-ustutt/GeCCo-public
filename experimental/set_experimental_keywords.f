*----------------------------------------------------------------------*
      subroutine set_experimental_keywords()
*----------------------------------------------------------------------*
*     set here additional sub-keyword of "calculate.experimental"
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'ifc_input.h'

      ! maximum perturbation order
      call argument_add('order','calculate.experimental',type=vtyp_int,
     &                  idef=(/0/))
      ! perturbation operator (XDIPLEN,YDIPLEN,ZDIPLEN)
      call argument_add('pert','calculate.experimental',type=vtyp_str,
     &                  len=7,cdef=(/'Z','D','I','P','L','E','N'/))
      ! irrep of perturbation operator
      call argument_add('pert_sym','calculate.experimental',
     &                  type=vtyp_int,idef=(/1/))
      ! frequency
      call argument_add('freq','calculate.experimental',
     &                  type=vtyp_rl8,xdef=(/0d0/))


      ! call keyword_add('new_kwd',context='calculate.experimental')
      ! see set_keywors for how to set up things ...

      return
      end

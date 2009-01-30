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
      ! perturbation operator (X,Y,Z)
      call argument_add('pert','calculate.experimental',type=vtyp_str,
     &                  len=1,cdef=(/'Z'/))
      ! frequency
      call argument_add('freq','calculate.experimental',
     &                  type=vtyp_rl8,xdef=(/0d0/))

      ! maximum excitation
      call argument_add('maxexc','calculate.experimental',type=vtyp_int,
     &                  idef=(/2/))

      ! treat V(1) intermediates as fock operator in approximation C
      call argument_add('treat_BV','calculate.experimental',
     &     type=vtyp_log,ldef=(/.true./))

      ! call keyword_add('new_kwd',context='calculate.experimental')
      ! see set_keywors for how to set up things ...

      return
      end

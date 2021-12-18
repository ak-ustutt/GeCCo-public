*----------------------------------------------------------------------*
      subroutine init_operator_0(op)
*----------------------------------------------------------------------*
*     allocate operator sub-arrays to 0/null()
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_operator.h'
      include 'def_orbinf.h'
      include 'ifc_memman.h'
      include 'par_globalmarks.h'

      type(operator), intent(inout) ::
     &     op

      integer ::
     &     nblk, nblkt, iblk, nsym, ndis, nexc, ncount, ifree

      op%ngas = 0
      op%nspin = 0

      op%name(1:len_opname) = ' '
      op%assoc_list(1:2*len_opname+2) = ' '

      op%order = -1
      op%species = -1
      op%dagger = .false.

      op%hermitian = 0

      op%ifreq => null()

      op%ica_occ => null()
      op%igasca_restr => null()
      op%formal_blk => null()
      op%blk_version => null()

      return
      end

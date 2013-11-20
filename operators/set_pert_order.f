*----------------------------------------------------------------------*
      subroutine set_pert_order(op,iorder,spec,freq_idx)
*----------------------------------------------------------------------*
*     set perturbation order, species of operator, and frequency indices
*     or operator block version
*     matthias, 2008
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'def_operator.h'

      type(operator), intent(inout) ::
     &     op

      integer, intent(in) ::
     &     iorder, spec, freq_idx(*)

      if (spec.gt.0) then
        op%order = iorder            ! perturbation order
c        print *,'order = ',iorder
        op%species = spec            ! operator species
c        print *,'species = ',spec
        if (associated(op%ifreq)) deallocate(op%ifreq)
        allocate(op%ifreq(iorder))
        op%ifreq=freq_idx(1:iorder)  ! frequency indices
c        print *,'ifreq = ',op%ifreq(1:iorder)
      else
        if (op%n_occ_cls.ne.iorder) then
          write(luout,*) 'number of blocks: ',op%n_occ_cls
          write(luout,*) 'your input:       ',iorder
           call quit(1,'set_pert_order',
     &                          'invalid argument value')
        end if
        if (associated(op%blk_version)) deallocate(op%blk_version)
        allocate(op%blk_version(iorder))
        op%blk_version=freq_idx(1:iorder) ! op. block versions
      end if

      return
      end

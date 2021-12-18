      subroutine auto_freeze(shell_def,nfreeze,orb_info)

      implicit none

      include 'stdunit.h'
      include 'def_orbinf.h'

      integer, intent(out) ::
     &     shell_def(*)
      integer, intent(inout) ::
     &     nfreeze
      type(orbinf), intent(in), target ::
     &     orb_info

      integer, pointer ::
     &     isym_bound_orbs(:)
      integer ::
     &     nfrz, idx, nsym

      nfrz = nfreeze
      if (nfreeze.lt.0) nfrz = orb_info%n_freeze_rcmd
      nfreeze = nfrz

      if (nfrz.gt.orb_info%n_bound_orbs)
     &     call quit(0,'auto_freeze',
     &     'freezing more orbitals than possible!?!')

      nsym = orb_info%nsym
      isym_bound_orbs => orb_info%isym_bound_orbs
      shell_def(1:nsym) = 0
      do idx = 1, nfrz
        shell_def(isym_bound_orbs(idx)) =
     &       shell_def(isym_bound_orbs(idx))+1
      end do

      return
      end

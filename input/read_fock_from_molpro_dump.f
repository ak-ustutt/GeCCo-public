      subroutine read_fock_from_molpro_dump(luinp,ecore,fock,nfock)

      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'
      include 'par_dalton.h'

      integer, intent(in) ::
     &     luinp
      integer, intent(in) ::
     &     nfock
      real(8), intent(inout) ::
     &     ecore,fock(nfock)

      real(8) ::
     &     val     
      integer ::
     &     iii, jjj, imo, jmo, ijmo


      real(8), external ::
     &     dnrm2

      ! file is assumed open
      ! fock is preinitialized with additional two-electron contributions
      ! the same is true for ecore

      do
        read(luinp,*,err=1234,end=1234) val, imo, jmo, iii, jjj

        if (imo.eq.0.and.jmo.eq.0) then
          ecore = ecore + val
          exit
        end if

        ijmo = (imo-1)*imo/2+jmo
        fock(ijmo) = fock(ijmo)+val

      end do

      if (dnrm2(nfock,fock,1).lt.1d-12)
     &   call quit(0,'read_fock_from_molpro_dump',
     &               'No sensible fock matrix found!')

      return

 1234 call quit(0,'read_fock_from_molpro_dump','unexpected EOF')
      end

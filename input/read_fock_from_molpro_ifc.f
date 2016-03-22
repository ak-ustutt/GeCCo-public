      subroutine read_fock_from_molpro_ifc(luinp,ifpos,
     &                                     ecore,fock,nfock,icoff)

      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'
      include 'par_dalton.h'

      integer, intent(in) ::
     &     luinp, ifpos, icoff
      integer, intent(in) ::
     &     nfock
      real(8), intent(inout) ::
     &     ecore,fock(nfock)

      real(8) ::
     &     val    
      integer(2) ::
     &     idxrd(4) 
      integer ::
     &     imo, jmo, ijmo, idxfl


      real(8), external ::
     &     dnrm2

      ! file is assumed open
      ! fock is preinitialized with additional two-electron contributions
      ! the same is true for ecore

      idxfl = ifpos
      do
        read(luinp,pos=idxfl,err=1234,end=1234) val, idxrd(1:4)
        idxfl = idxfl+16

        imo = idxrd(1);  jmo = idxrd(2) 

        if (imo.eq.0.and.jmo.eq.0) then
          ecore = ecore + val
          exit
        end if

        imo = imo + icoff; jmo = jmo + icoff

        ijmo = (imo-1)*imo/2+jmo
        fock(ijmo) = fock(ijmo)+val

      end do

      if (dnrm2(nfock,fock,1).lt.1d-12)
     &   call quit(0,'read_fock_from_molpro_ifc',
     &               'No sensible fock matrix found!')

      return

 1234 call quit(0,'read_fock_from_molpro_ifc','unexpected EOF')
      end

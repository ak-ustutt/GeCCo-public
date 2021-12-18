      subroutine read_fock_from_file_gms(eref,fock,nfock)

      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'
      include 'par_gamess.h'

      integer, intent(in) ::
     &     nfock
      real(8), intent(out) ::
     &     eref,fock(nfock)

      type(filinf) ::
     &     ffinp
      integer ::
     &     luinp, ii

      ! 64 bit version of GAMESS writes integer*8
      integer(8) ::
     &     irecst,ioda(len_da),ifilen(len_da),is,ipk
      real(8) ::
     &     potnuc,emcscf,sz,s2,etot2,eerd,e1,e2,ven,vee,epot,
     &     eelct,ekin,estate(mxrt),statn,ecore

      real(8), external ::
     &     dnrm2

      ! first we get the reference energy from DICTNRY
      call file_init(ffinp,dictnry,ftyp_da_unf,irecln)
      call file_open(ffinp)

      luinp = ffinp%unit

      ! read first record: information about where to find which record
      read (luinp,rec=1) irecst,ioda,ifilen,is,ipk

      ! read energy quantities
      read (luinp,rec=ioda(6)) potnuc,eelct,etot2,sz,s2,ecore,emcscf,
     &                          eerd,e1,e2,ven,vee,epot,ekin,
     &                          estate,statn

      ! reference energy is sum of core energy and nuclear potential energy:
      eref = ecore + potnuc
      call file_close_keep(ffinp)

      ! now we read fock matrix from MOINTS
      call file_init(ffinp,moints,ftyp_sq_unf,0)
      call file_open(ffinp)

      fock(1:nfock) = 0d0

      luinp = ffinp%unit
      rewind luinp
      
      read(luinp) fock(1:nfock)

      call file_close_keep(ffinp)

      if (dnrm2(nfock,fock,1).lt.1d-12)
     &   call quit(0,'read_fock_from_file_gms',
     &               'No sensible fock matrix found!')

      return
      end

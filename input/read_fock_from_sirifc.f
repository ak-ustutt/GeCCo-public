      subroutine read_fock_from_sirifc(eref,fock,nfock)

      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'
      include 'par_dalton.h'

      integer, intent(in) ::
     &     nfock
      real(8), intent(out) ::
     &     eref,fock(nfock)

      type(filinf) ::
     &     ffsir
      logical ::
     &     closeit
      real(8) ::
     &     potnuc,emy,eactiv,emcscf
      ! DALTON writes integer*4, so we must take care of that
      integer(4) ::
     &     istate,ispin,nactel,lsym
      
      integer ::
     &     lusir, luerr, irec

      real(8), external ::
     &     dnrm2

      ! open files
      call file_init(ffsir,sirifc,ftyp_sq_unf,0)
      call file_open(ffsir)

      lusir = ffsir%unit
      rewind lusir
      
      luerr = luout
      call mollab('SIR IPH ',lusir,luerr)

      read (lusir) potnuc,emy,eactiv,emcscf,istate,ispin,nactel,lsym
      ! overread a few records (depends on DALTON version)
      do irec = 1, nskip_in_sirifc-1
        read(lusir) 
      end do
      ! and finally: the inactive fock matrix in symmetry-blocked
      ! upper triangular form
      read (lusir,err=16) fock(1:nfock)
      ! first element must be negative
      if (fock(1).ge.0d0) goto 16
      goto 1
      ! on error try next record
 16   write(luout,*) 'Trying new DALTON format ...'
      read (lusir) fock(1:nfock) 

  1   call file_close_keep(ffsir)

c dbg
c      print *,'fock',fock(1:nfock)
c dbg

      if (dnrm2(nfock,fock,1).lt.1d-12)
     &   call quit(0,'read_fock_from_sirifc',
     &               'No sensible fock matrix found!')

      if (nactel.eq.0) then
        eref = emcscf
      else
        eref = emy + potnuc    ! core energy
      end if

      return
      end

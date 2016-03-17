      subroutine read_fock_from_file_cfour(ecore,fock,nfock,fname)

      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'
      include 'par_dalton.h'

      integer, intent(in) ::
     &     nfock
      real(8), intent(inout) ::
     &     ecore,fock(nfock)
      character(len=*), intent(in) ::
     &     fname

      real(8) ::
     &     val    
      integer(2) ::
     &     idxrd(4) 
      integer ::
     &     lu, imo, jmo, ijmo, idxfl

      type(filinf) :: fffock

      real(8), external ::
     &     dnrm2


      call file_init(fffock,trim(fname),ftyp_sq_frm,0)
      call file_open(fffock)
      lu = fffock%unit

      read(lu,*) ecore

      do
        read(lu,*,err=1234,end=4321) imo, jmo, val

        !imo = imo + icoff; jmo = jmo + icoff

        ijmo = (imo-1)*imo/2+jmo
        fock(ijmo) = val

      end do

 4321 continue


      if (dnrm2(nfock,fock,1).lt.1d-12)
     &   call quit(0,'read_fock_from_file_cfour',
     &               'No sensible fock matrix found!')

      call file_close_keep(fffock)

      return

 1234 call quit(0,'read_fock_from_file_cfour','Error during reading')
      end

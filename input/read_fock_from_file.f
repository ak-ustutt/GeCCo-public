      subroutine read_fock_from_file(eref,fock,nfock,name)

      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'
      include 'par_dalton.h'

      character*(*), intent(in) ::
     &     name
      integer, intent(in) ::
     &     nfock
      real(8), intent(out) ::
     &     eref,fock(nfock)

      type(filinf) ::
     &     ffinp
      integer(8) ::
     &     blk_len
      real(8) ::
     &     escf, enuc      
      integer ::
     &     luinp, iii, imo, jmo, ijmo

      integer(2), pointer ::
     &     idx(:,:)
      real(8), pointer ::
     &     val(:)

      real(8), external ::
     &     dnrm2

      ! open files
      call file_init(ffinp,trim(name),ftyp_sq_unf,0)
      call file_open(ffinp)

      fock(1:nfock) = 0d0

      luinp = ffinp%unit
      rewind luinp
      
      read(luinp) blk_len, escf, enuc

      eref = escf

      allocate(idx(2,blk_len),val(blk_len))

      do 

        read(luinp) blk_len,idx(1:2,1:blk_len),val(1:blk_len)
        if (blk_len.le.0) exit
        
        do iii = 1, blk_len
          imo = idx(1,iii)
          jmo = idx(2,iii)
          ijmo = (imo-1)*imo/2+jmo
          fock(ijmo) = val(iii)
        end do

      end do

      deallocate(idx,val)

      call file_close_keep(ffinp)

      if (dnrm2(nfock,fock,1).lt.1d-12)
     &   call quit(0,'read_fock_from_file',
     &               'No sensible fock matrix found!')

      return
      end

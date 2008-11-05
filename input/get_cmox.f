      subroutine get_cmox(cmo,unit)

      implicit none

      real(8), intent(out) ::
     &     cmo(*)
      integer, intent(in) ::
     &     unit

      integer(8) ::
     &     nblk, nrow(8), ncol(8)

      integer ::
     &     idxst, idxnd, iblk

      rewind(unit)
      read(unit) nblk,nrow(1:nblk),ncol(1:nblk)

      idxst = 1
      do iblk = 1, nblk
        idxnd = idxst-1 + nrow(iblk)*ncol(iblk)
        read(unit) cmo(idxst:idxnd)
        idxst = idxnd+1
      end do

      end

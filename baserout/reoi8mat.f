*----------------------------------------------------------------------*
      subroutine reoi8mat(imat,ireo,nrow,ncol,mode)
*----------------------------------------------------------------------*

      ! mode = 1: reorder rows (= within each column) dim(ireo) = ncol
      ! mode = 2: reorder colums (= within each row)  dim(ireo) = nrow
      ! mode = 3: both rows and columns (nrow==ncol ! )

      implicit none

      integer, intent(in) ::
     &     mode, nrow, ncol, ireo(*)
      integer(8), intent(inout) ::
     &     imat(nrow,ncol)

      integer(8) ::
     &     iscr(max(nrow,ncol))
      integer ::
     &     idx, jdx

      if (mode.lt.1.or.mode.gt.3) then
        call quit(1,'reoi8mat','unknown mode parameter')
      end if
        
      if (mode.eq.1.or.mode.eq.3) then
        do idx = 1, ncol
          do jdx = 1, nrow
            iscr(jdx) = imat(ireo(jdx),idx)
          end do
          imat(1:nrow,idx) = iscr(1:nrow)
        end do
      end if

      if (mode.eq.2.or.mode.eq.3) then
        do jdx = 1, nrow
          do idx = 1, ncol
            iscr(idx) = imat(jdx,ireo(idx))
          end do
          imat(jdx,1:ncol) = iscr(1:ncol)
        end do
      end if

      return
      end

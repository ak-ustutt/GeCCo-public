
      subroutine optc_trmat(imet,isym,xmat,ndim,xmat0,ndim0,ndel,ld0)

* obtain reduced space matrix in new metric
* if matrix was symmetry packed: unpack it on the way      

      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     imet, isym, ndim, ndim0, ndel, ld0
      real(8), intent(in) ::
     &     xmat0(*)
      real(8), intent(out) ::
     &     xmat(*)

      integer ::
     &     ii, jj, ii0, jj0, iioff, iioff0, iip1f0, iindf0, jjnd,
     &     idxijp1, idxin

      if (ntest.gt.0) then
        write(lulog,*) '------------'
        write(lulog,*) ' optc_trmat'
        write(lulog,*) '------------'
        write(lulog,*) ' ndim, ndim0, ld0 : ',ndim, ndim0, ld0
        write(lulog,*) ' ndel : ',ndel 
      end if

      if (ntest.ge.10) then
        write(lulog,*) ' imet, isym: ',imet,isym
        write(lulog,*) ' xmat0: '
        if (isym.eq.0) then
          call wrtmat2(xmat0,ndim0,ndim0,ld0,ld0)
        else
          call prtrlt(xmat0,ndim0)
        end if
      end if

      if (imet.eq.0) then
      
        do ii = 1, ndim

          ii0 = ndel+ii

          if (isym.eq.0) then
            iioff  = (ii-1)*ndim
            iioff0 = (ii0-1)*ld0
            iip1f0 = (ii0  )*ld0
            jjnd = ndim
          else
            iioff  = (ii-1)*ndim
            iioff0 = (ii0-1)*ii0/2
            iip1f0 = ii0*(ii0+1)/2
            jjnd = ii
          end if
          do jj = 1, jjnd
            jj0 = jj+ndel

            if (isym.eq.0) then
              idxijp1 = iioff0+jj0+1
            else if (ii0.ge.jj0+1) then
              idxijp1 = iioff0+jj0+1
            else
              idxijp1 = jj0*(jj0+1)/2 + ii0
            end if
              
            xmat(iioff + jj) = xmat0(iioff0+jj0) 
     &                       - xmat0(idxijp1)
     &                       - xmat0(iip1f0+jj0)
     &                       + xmat0(iip1f0+jj0+1)

          end do
        end do

      else

        do ii = 1, ndim

          ii0 = ndel+ii

          if (isym.eq.0) then
            iioff  = (ii-1)*ndim
            iioff0 = (ii0-1)*ld0
            iindf0 = (ndim0-1)*ld0
            jjnd = ndim
          else
            iioff  = (ii-1)*ndim
            iioff0 = (ii0-1)*ii0/2
            iindf0 = (ndim0-1)*ndim0/2
            jjnd = ii
          end if
          do jj = 1, jjnd
            jj0 = jj+ndel

            if (isym.eq.0) then
              idxin = iioff0 + ndim0
            else
              idxin = iindf0 + ii0
            end if
            xmat(iioff + jj) = xmat0(iioff0+jj0) 
     &                       - xmat0(idxin)
     &                       - xmat0(iindf0+jj0)
     &                       + xmat0(iindf0+ndim0)

          end do
        end do
        
      end if

      if (isym.ne.0) then
        do ii = 1, ndim
          iioff = (ii-1)*ndim
          do jj = 1, ii-1
            xmat((jj-1)*ndim+ii) = xmat(iioff+jj) 
          end do
        end do
        
      end if

      if (ntest.ge.10) then
        write(lulog,*) ' final xmat: '
        call wrtmat2(xmat,ndim,ndim,ndim,ndim)
      end if

      return
      end

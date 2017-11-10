*----------------------------------------------------------------------*
      subroutine svd_get_left(nrow,ncol,mat)
*----------------------------------------------------------------------*
*     performs a singular value decomposition of a rectangular matrix
*     and returns the nonredundant part of the left unitary matrix
*
*     matthias, jun 2013
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'routes.h'
      include 'ifc_memman.h'

      integer, parameter ::
     &     ntest = 00
      real(8), parameter ::
     &     warn_sv = 1d-12 ! give a warning for small singular values
      integer, intent(in) ::
     &     nrow, ncol
      real(8), intent(inout) ::
     &     mat(nrow,ncol)

      real(8) ::
     &     singval(min(nrow,ncol)),
     &     dum1,dum2
      real(8),pointer ::
     &      wrk(:)
      integer ::
     &     idx, lwrk, info,
     &     ifree

      if (nrow*ncol.eq.0) return
      if (nrow.lt.ncol) call quit(1,'svd_get_left',
     &       'not yet tested for nrow < ncol')

      ifree = mem_setmark('svd_get_left')
      lwrk=max(1024,max(nrow,ncol)**2)
      ifree = mem_alloc_real(wrk, lwrk,'wrk')
      info = 0

      if (ntest.ge.100) then
        write(lulog,'(x,a)') '--------------------'
        write(lulog,'(x,a)') 'svd_get_left at work'
        write(lulog,'(x,a)') '--------------------'
        write(lulog,*) 'input matrix:'
        call wrtmat2(mat,nrow,ncol,nrow,ncol)
      end if

      ! calculate left-hand U and singular values:
      call dgesvd('O','N',nrow,ncol,
     &     mat,nrow,singval,
     &     dum1,1,dum2,1,
     &     wrk,lwrk,info)

      if (info.ne.0) then
        write(lulog,*) 'WARNING in svd_get_left: SVD in trouble'
      end if

      if (ntest.ge.100) then
        write(lulog,*) 'left-hand eigenvector matrix U:'
        call wrtmat2(mat,nrow,ncol,nrow,ncol)
        write(lulog,*) 'singular values:'
        call wrtmat2(singval,min(nrow,ncol),1,min(nrow,ncol),1)
      end if

      ! zero the columns with negligible singular values
      do idx = 1, min(nrow,ncol)
        if (singval(idx).lt.sv_thresh) then
          mat(1:nrow,idx) = 0d0
          if (singval(idx).gt.warn_sv) then
            write(lulog,*) 'small singular value: ',singval(idx)
            call warn('svd_get_left','small singular values')
          end if
        end if
      end do
      ifree = mem_flushmark()


      return
      end

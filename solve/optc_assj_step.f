
*----------------------------------------------------------------------*
      subroutine optc_assj_step(imet,xvgmat,xvvmat,xvec,ndim,mxdim,ndel,
     &     idamp,trrad,xamat,xsmat)
*----------------------------------------------------------------------*
*
*     idamp = 0 :  no damping (no S-matrix needed)
*     idamp = 1 :  w. damping (S-matrix needed)
*
*     imet = 0 :  use |t^i>-|t^{i+1}> as basis
*     imet = 1 :  use |t^i>-|t^{n}> as basis
*
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 100

      integer, intent(in) ::
     &     ndim, mxdim, imet, idamp
      integer, intent(out) ::
     &     ndel
      real(8), intent(in) ::
     &     trrad, xvvmat(mxdim,*), xvgmat(mxdim,*)
      real(8), intent(out) ::
     &     xvec(*)
      real(8), intent(inout) ::
     &     xamat(*), xsmat(*)

      logical ::
     &     converged, again
      integer ::
     &     ixd_iter, ierr, irg, job, ncdim,
     &     ii, jj
      real(8) ::
     &     xdlt, xdlt_min, cond, xnrm, eig_min, xval, xdum

* local O(N)-scratch
      integer ::
     &     kpiv(ndim)
      real(8) ::
     &     xscrv(ndim), xscrv2(ndim), xscrv3(ndim)

      if (ntest.ge.0) then
        write(luout,*) '===================='
        write(luout,*) ' ASSJ modul entered'
        write(luout,*) '===================='
        write(luout,*) ' ndim : ', ndim
      end if

      converged = .false.
      ixd_iter = 0
      xdlt = 0d0
      xdlt_min = 0d0
      ncdim = ndim-1
      ndel = 0

      do while(.not.converged.and.ncdim.gt.0)

        if (idamp.gt.0) then
          ! set up reduced S-matrix
          call optc_trmat(imet,1,xsmat,ncdim,xvvmat,ndim,ndel,mxdim)

          ! test for positive definite metric
          if (ixd_iter.eq.0) then
            irg = 0
            call rg(ncdim,ncdim,xsmat,xscrv,xscrv2,irg,
     &           xdum,kpiv,xscrv3,ierr)
            if (ierr.ne.0) then
              write(luout,*) 'error code from rg: ',ierr
              stop 'optc_assj_step'
            end if
            if (ntest.ge.20) then
              write(luout,*) 'eigenvalues of metric:'
              do ii = 1, ncdim
                write(luout,*) ii,xscrv(ii),
     &               xscrv2(ii)
              end do
            end if

            eig_min = +100d0
            do ii = 1, ncdim
              eig_min = min(eig_min,xscrv(ii))
            end do

            if (eig_min.le.0d0) then
              ndel = ndel+1
              ncdim = ncdim-1
              cycle
            end if

            ! restore reduced S-matrix
            call optc_trmat(imet,1,xsmat,ncdim,xvvmat,ndim,ndel,mxdim)
          end if
        end if

        ! set up reduced A-matrix
        call optc_trmat(imet,0,xamat,ncdim,xvgmat,ndim,ndel,mxdim)
        

        ! set up reduced gradient
        call optc_trvec(imet,xvec,ncdim,xvgmat(ndel+1,ndim))

        ! check for negative eigenvalues:
        if (idamp.gt.0.and.ixd_iter.eq.0) then
          irg = 0 
          call rgg(ncdim,ncdim,xamat,xsmat,
     &         xscrv,xscrv2,xscrv3,irg,xdum,ierr)
          if (ierr.ne.0) then
            write(luout,*) 'error code from rgg: ',ierr
            stop 'optc_assj_step'
          end if
          if (ntest.ge.20) then
            write(luout,*) 'eigenvalues of reduced approx. jacobian:'
            do ii = 1, ncdim
              write(luout,*) ii,xscrv(ii)/xscrv3(ii),
     &                     xscrv2(ii)/xscrv3(ii)
            end do
          end if

          eig_min = 0d0
          do ii = 1, ncdim
            eig_min = min(eig_min,xscrv(ii)/xscrv3(ii))
          end do

          if (eig_min.lt.0d0) then
            xdlt = -eig_min+0.1d0
            xdlt_min = -eig_min
          end if

          ! rebuild A and S:
          call optc_trmat(imet,0,xamat,ncdim,xvgmat,ndim,ndel,mxdim)        
          call optc_trmat(imet,1,xsmat,ncdim,xvvmat,ndim,ndel,mxdim)
        end if

        if (idamp.gt.0.and.xdlt.ne.0d0) then
          ! make [A + dlt S]
          xamat(1:ncdim*ncdim) = xamat(1:ncdim*ncdim)
     &                     +xdlt*xsmat(1:ncdim*ncdim)
        end if

        ! solve subspace problem
        call dgeco(xamat,ncdim,ncdim,kpiv,cond,xscrv)
        if (ntest.ge.10) write(luout,*)'dimension ,condition number: ',
     &         ncdim,cond

        if (ntest.ge.50) then
          write(luout,*) 'Factorized matrix from dgeco:'
          call wrtmat2(xamat,ncdim,ncdim,ncdim,ncdim)
          write(luout,*) 'Pivot array:'
          call iwrtma(kpiv,1,ncdim,1,ncdim)
        end if

        again = ixd_iter.eq.0.and.cond.lt.1d-14

        if (again) then
          ndel = ndel+1
          ncdim = ncdim-1
          cycle
        end if

        if (ntest.ge.20) then
          write(luout,*) 'RHS:'
          call wrtmat2(xvec,ncdim,1,ncdim,1)
        end if

        job = 0
        call dgesl(xamat,ncdim,ncdim,kpiv,xvec,job)

        if (ntest.ge.20) then
          write(luout,*) 'Resulting vector:'
          call wrtmat2(xvec,ncdim,1,ncdim,1)
        end if
        
        if (idamp.gt.0) then
          xnrm = 0d0
          do ii = 1, ncdim
            do jj = 1, ncdim
              xnrm = xnrm + xvec(ii)*xsmat((ii-1)*ncdim+jj)*xvec(jj)
            end do
          end do
          xnrm = sqrt(xnrm)
          
          ! call control kernel
          xval = xnrm - trrad
          call optc_xdamp_ctl(xdlt,xdlt_min,xval,converged,ixd_iter)

        else
          converged = .true.

        end if

      end do

      if (ncdim.eq.0) then
        xvec(1:mxdim) = 0d0
      else
        xscrv(1:ncdim) = xvec(1:ncdim)
        xvec(1:ndel) = 0d0
        if (imet.eq.0) then
          xvec(ndel+1) = - xscrv(1)
          do ii = 2, ncdim
            xvec(ii+ndel) = xscrv(ii-1)-xscrv(ii)
          end do
          xvec(ndim) = xscrv(ncdim)
        else
          xvec(ndel+1:ndel+ncdim) = -xscrv(1:ncdim)
          xvec(ndim) = 0d0
          do ii = 1, ncdim
            xvec(ndim) = xvec(ndim) + xscrv(ii)
          end do
        end if

      end if
      
      if (ntest.ge.20) then
        write(luout,*) 'Re-transformed vector (ndel = ',ndel,'):'
        call wrtmat2(xvec,ndim,1,ndim,1)
      end if

      return
      end

*----------------------------------------------------------------------*
      subroutine optc_prc_special2_inner
     &     (grad, beyond_A, scale, w_shift,
     &      scr1, scr2,  bmat, xmat, xdia,
     &      ndim_bx, ndim_c, ndim_a,
     &      nidx_cstr,ms_cstr,gam_cstr,restr_cstr,mostnd_cstr,ngas_cstr,
     &      nidx_astr,ms_astr,gam_astr,restr_astr,mostnd_astr,ngas_astr,
     &      igamorb,ngam,ngas)
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      logical ::
     &     beyond_A, scale
      integer, intent(in) ::
     &     ndim_bx, ndim_c, ndim_a,
     &     nidx_cstr,ms_cstr,gam_cstr,ngas_cstr,
     &     nidx_astr,ms_astr,gam_astr,ngas_astr,ngam,ngas
      integer, intent(in) ::
     &     igamorb(*),
     &     restr_cstr(2,ngas,2), mostnd_cstr(2,ngam,ngas_cstr),
     &     restr_astr(2,ngas,2), mostnd_astr(2,ngam,ngas_astr)
      real(8), intent(in) ::
     &     w_shift
      real(8), intent(in) ::
     &     bmat(ndim_bx,ndim_bx), xmat(ndim_bx,ndim_bx), xdia(*)
      real(8), intent(inout) ::
     &     scr1(ndim_bx,ndim_a), scr2(ndim_bx,ndim_bx),
     &     grad(ndim_c,ndim_bx,ndim_a)

      logical ::
     &     first
      integer ::
     &     idx_c, idx_a, idx_b1, idx_b2, istr, idx, jdx
      real(8) ::
     &     xsum_c(ndim_c), xsum_a(ndim_a), orbsum
      integer ::
c     &     idxorb(max(ndim_c,ndim_a)),
c     &     idxspn(max(ndim_c,ndim_a)),idxdss(max(ndim_c,ndim_a))
     &     idxorb(20),
     &     idxspn(20),idxdss(20)

      logical, external ::
     &     next_string

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'optc_prc_special2_inner')
        write(lulog,*) 'beyond_A: ',beyond_A
        write(lulog,*) 'ndim_bx, ndim_c, ndim_a: ',
     &                  ndim_bx, ndim_c, ndim_a
      end if

      if (.not.beyond_A) then
        do idx_c = 1, ndim_c
          
          ! MM-multiply
          do idx_a = 1, ndim_a
            do idx_b1 = 1, ndim_bx
              scr1(idx_b1,idx_a) = 0d0
              do idx_b2 = 1, ndim_bx
                scr1(idx_b1,idx_a) =
     &               scr1(idx_b1,idx_a) +
     &               bmat(idx_b1,idx_b2)*grad(idx_c,idx_b2,idx_a)
              end do
            end do
          end do

          ! scatter scr1 to grad
          do idx_a = 1, ndim_a
            do idx_b1 = 1, ndim_bx
              grad(idx_c,idx_b1,idx_a) = scr1(idx_b1,idx_a)
            end do
          end do
          
        end do

      else

        if (nidx_cstr.gt.0) then
          first = .true.
          istr = 0
          ! set up orbital sum for C-strings
          cstr_loop: do
            if(.not.next_string(idxorb,idxspn,idxdss,
     &                 nidx_cstr,ms_cstr,gam_cstr,first,
     &                 restr_cstr,mostnd_cstr,
     &                 igamorb,ngam,ngas_cstr)) exit cstr_loop
            first = .false.
            istr = istr+1
            xsum_c(istr) = 0d0
            do idx = 1, nidx_cstr
              xsum_c(istr) = xsum_c(istr) + xdia(idxorb(idx))
            end do
          end do cstr_loop
        else
          xsum_c(1) = 0d0
        end if

        first = .true.
        istr = 0
        ! set up orbital sum for A-strings
        astr_loop: do
          if(.not.next_string(idxorb,idxspn,idxdss,
     &                 nidx_astr,ms_astr,gam_astr,first,
     &                 restr_astr,mostnd_astr,
     &                 igamorb,ngam,ngas_astr)) exit astr_loop
          first = .false.
          istr = istr+1
          xsum_a(istr) = 0d0
          do idx = 1, nidx_astr
            xsum_a(istr) = xsum_a(istr) + xdia(idxorb(idx))
          end do
        end do astr_loop

        do idx_c = 1, ndim_c
          do idx_a = 1, ndim_a
            orbsum = xsum_c(idx_c)-xsum_a(idx_a) + w_shift
            ! effective Beff = B + (eps_C - eps_A)X
            scr2 = bmat + orbsum*xmat 

c test
            if (scale) then
              scr2 = 4d0*scr2
c              do idx = 1, ndim_bx
c                do jdx = 1, idx-1
c                  scr2(idx,jdx) = 0d0
c                  scr2(jdx,idx) = 0d0
c                end do
c                scr2(idx,idx) = scr2(idx,idx)
c              end do
            end if
c
            ! invert
            call inv_svd(scr1(1,idx_a),scr2,grad(idx_c,1,idx_a),
     &                   ndim_bx,ndim_c)

c            call gaussj(scr2,ndim_bx,ndim_bx)

c            ! multiply with gradient
c            do idx_b1 = 1, ndim_bx
c              scr1(idx_b1,idx_a) = 0d0
c              do idx_b2 = 1, ndim_bx
c                scr1(idx_b1,idx_a) =
c     &               scr1(idx_b1,idx_a) +
c     &               scr2(idx_b1,idx_b2)*grad(idx_c,idx_b2,idx_a)
c              end do
c            end do
            
          end do

          ! scatter scr1 to grad
          do idx_a = 1, ndim_a
            do idx_b1 = 1, ndim_bx
              grad(idx_c,idx_b1,idx_a) = scr1(idx_b1,idx_a)
            end do
          end do

        end do

      end if

      return
      end

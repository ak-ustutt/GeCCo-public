*----------------------------------------------------------------------*
      subroutine expand_occ(occ,
     &                      occ_tmpl,
     &                      nc, na,
     &                      occ_csub, occ_asub,
     &                      hpvx_csub,hpvx_asub,
     &                      njoined)
*----------------------------------------------------------------------*
*     expand occupation-like information in condensed arrays 
*     occ_csub/occ_asub to full occupation form;
*     a template occupation is needed, however, in order to reconstruct
*     the actual ordering in case of njoined > 1
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     nc, na,
     &     occ_csub(nc), occ_asub(na), 
     &     hpvx_csub(nc), hpvx_asub(na),
     &     njoined, occ_tmpl(ngastp,2,njoined)
      integer, intent(out) ::
     &     occ(ngastp,2,njoined)

      integer ::
     &     ijoin, hpvx, idx_hpvx, idx


      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'expand_occ')
        write(luout,'(3x,"C: occ  ",10i5)') occ_csub(1:nc) 
        write(luout,'(3x,"   hpvx ",10i5)') hpvx_csub(1:nc) 
        write(luout,'(3x,"A: occ  ",10i5)') occ_asub(1:na) 
        write(luout,'(3x,"   hpvx ",10i5)') hpvx_asub(1:na) 
      end if

      occ(1:ngastp,1:2,1:njoined) = 0
      do idx = 1, nc
        hpvx = hpvx_csub(idx)
        do ijoin = 1, njoined
          if (occ_tmpl(hpvx,1,ijoin).gt.0.and.
     &        occ     (hpvx,1,ijoin).eq.0     ) then
            occ(hpvx,1,ijoin) = occ_csub(idx)
            exit
          end if
        end do
      end do
      do idx = 1, na
        hpvx = hpvx_asub(idx)
        do ijoin = 1, njoined
          if (occ_tmpl(hpvx,2,ijoin).gt.0.and.
     &        occ     (hpvx,2,ijoin).eq.0     ) then
            occ(hpvx,2,ijoin) = occ_asub(idx)
            exit
          end if
        end do
      end do
      
      if (ntest.ge.100) then
        call wrt_occ_n(luout,occ,njoined)
      end if

      return
      end
*----------------------------------------------------------------------*

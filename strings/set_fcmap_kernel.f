*----------------------------------------------------------------------*
      subroutine set_fcmap_kernel(strmap,
     &     iocc,irestr1,irestr2,
     &     idxms,igam,
     &     g2_y4sg,g2_yinf,
     &     g2_yssg,g2_wssg,
     &     g2_off_dgm,g2_ndis,
     &     mostnd_cur,nsym,ngas_cur,igamorb)
*----------------------------------------------------------------------*
*     core routine for setting up the (str)->(str_spin_flipped) mapping
*
*     andreas, july 2008
*      
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 000

      integer, intent(out) ::
     &     strmap(*)
      integer, intent(in) ::
     &     iocc, irestr1(2,ngas_cur,2), irestr2(2,ngas_cur,2),
     &     idxms, igam,
     &     nsym, ngas_cur,
     &     mostnd_cur(2,nsym),igamorb(*),
     &     g2_y4sg(*),g2_yinf(*),
     &     g2_yssg(*),g2_wssg(*),
     &     g2_off_dgm(*),g2_ndis

      logical ::
     &     first
      integer ::
     &     idorb(iocc),idspn(iocc),idgam(iocc),idss(iocc),
     &     istr, nstr, idx, ms, idxmap, isgn

      integer, external ::
     &     idx4sg
      logical, external ::
     &     next_string, allow_sbsp_dis

      ms = iocc-(idxms-1)*2

      idxmap = 0

      first = .true.
      ! loop over strings of graph1
      do while(next_string(idorb,idspn,idss,
     &     iocc,ms,igam,first,
     &     irestr1,
     &     mostnd_cur,igamorb,
     &     nsym,ngas_cur))

        first = .false.

        idxmap = idxmap+1

        ! allowed for new restriction??
        if(.not.allow_sbsp_dis(idss,iocc,ngas_cur,irestr2)) then
          strmap(idxmap) = 0
          cycle
        end if

        ! set symmetry info
        do idx = 1, iocc
          idgam(idx) = igamorb(idorb(idx))
        end do

        ! get index on string with different restriction
        strmap(idxmap) = 
     &        (idx4sg(iocc,idss,idorb,idspn,idgam,
     &                g2_y4sg,g2_yinf,
     &                g2_yssg,g2_wssg,
     &                g2_off_dgm,g2_ndis,
     &                mostnd_cur,
     &                iocc,nsym,ngas_cur)+1)

      end do

      if (ntest.ge.1000) then
        write(lulog,*) 'set_fcmap_kernel: the result of my work is'
        nstr = idxmap
        do istr = 1, nstr
          write(lulog,'(x,i5,"->",i5)') istr,strmap(istr)
        end do
      end if

      return
      end


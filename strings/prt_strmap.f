*----------------------------------------------------------------------*
      subroutine prt_strmap(strmap,iocc1,iocc2,nstr1,nstr2)
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'hpvxseq.h'

      integer, intent(in) ::
     &     strmap(*),
     &     iocc1(ngastp), iocc2(ngastp), nstr1(ngastp), nstr2(ngastp)

      integer ::
     &     imap, idx, igastp, hpvx

      write(lulog,'(x,a,4i4,a,4i4,a)')
     &               'maps for: (',iocc1(1:ngastp),' )(',
     &                             iocc2(1:ngastp),' )'
      imap = 0
      idx = 1
      do igastp = 1, ngastp
        hpvx = hpvxseq(igastp)
        if (iocc1(hpvx)+iocc2(hpvx).eq.0) cycle
        imap = imap+1
        write(lulog,*) 'map #',imap
        call wrtimat2(strmap(idx),nstr1(hpvx),nstr2(hpvx),
     &                            nstr1(hpvx),nstr2(hpvx))
        idx = idx + nstr1(hpvx)*nstr2(hpvx)
      end do

      return
      end

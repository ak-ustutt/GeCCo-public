*----------------------------------------------------------------------*
      subroutine set_ientry(ientry,iwexit,maxms,ngam,nel)
*----------------------------------------------------------------------*
*     another little wrapper routine
*----------------------------------------------------------------------*

      integer, intent(in) ::
     &     maxms,ngam,nel,
     &     iwexit(-maxms:maxms,ngam,0:nel)

      integer, intent(out) ::
     &     ientry(-maxms:maxms,ngam)

c dbg
c      print *,'iwexit passed to subr:'
c      call iwrtma(iwexit,maxms*2+1,nel+1,maxms*2+1,nel+1)
c dbg

      ! get appropriate line from iwexit
      ientry(-maxms:maxms,1:ngam) = iwexit(-maxms:maxms,1:ngam,nel)

      return
      end

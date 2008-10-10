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

      ! get appropriate line from iwexit
      ientry(-maxms:maxms,1:ngam) = iwexit(-maxms:maxms,1:ngam,nel)

c dbg
      print *,'iwexit',iwexit(-maxms:maxms,1:ngam,nel)
c dbg
      return
      end

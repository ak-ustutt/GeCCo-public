

*----------------------------------------------------------------------*
*     a few auxiliary routines follow
*----------------------------------------------------------------------*
      subroutine icp_save(isrc,idst,n)
      
      implicit none

      integer, intent(in) ::
     &     n
      integer, intent(out) ::
     &     idst(n)
      integer, intent(in) ::
     &     isrc(n)

      integer ::
     &     ibuff(n)

      ibuff(1:n) = isrc(1:n)
      idst(1:n) = ibuff(1:n)

      return
      end

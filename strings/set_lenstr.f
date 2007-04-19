*----------------------------------------------------------------------*
      subroutine set_lenstr(lenstr,idis,iw4sg,ndis,nel,norb,ngam,maxms)
*----------------------------------------------------------------------*
*     one more wrapper for resorting information
*     get length of strings per distribution, IRREP, Ms from
*     last row in weight array
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     ndis, idis, nel, ngam, maxms, norb,
     &     iw4sg(0:nel,-maxms:maxms,ngam,norb)

      integer, intent(out) ::
     &     lenstr(ndis,ngam,1:maxms+1)

      integer ::
     &     ims_cnt, ims, igam

      do ims_cnt = 1, maxms+1
        ims = -maxms+2*(ims_cnt-1)
        do igam = 1, ngam
          lenstr(idis,igam,ims_cnt) = iw4sg(nel,ims,igam,norb)
        end do
      end do

      return
      end

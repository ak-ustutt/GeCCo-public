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
c dbg
c          print *,'for lenstr',idis,igam,ims_cnt
c          print *,'  looking at iw4sg at ',nel,ims,igam,norb
c          print *,'  value = ',iw4sg(nel,ims,igam,norb)
c dbg
          lenstr(idis,igam,ims_cnt) = iw4sg(nel,ims,igam,norb)
        end do
      end do

c dbg
c      print *,'for idis = ',idis 
c      print *,'lenstr = ',lenstr(idis,1:ngam,1:maxms+1)
c dbg

      return
      end

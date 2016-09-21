*----------------------------------------------------------------------*
      subroutine optc_fix_signs2(ff_vec,irecvec,
     &     opti_info,iopt,
     &     nwfpar,xbuf1)
*----------------------------------------------------------------------*
*
*  Matthias' sign fix for internal contraction stuff
*  version that just applies sign-fix to a single vector
*     
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'
      include 'def_optimize_info.h'

      integer, parameter ::
     &     ntest = 1000

      integer, intent(in) ::
     &     irecvec
      type(filinf) ::
     &     ff_vec
      integer, intent(in) ::
     &     nwfpar, iopt
      real(8), intent(inout) ::
     &     xbuf1(*)
      type(optimize_info) ::
     &     opti_info


      integer ::
     &     ioff, nsec, idx, irec, isec
      integer, pointer ::
     &     nwfpsec(:), idstsec(:)
      real(8), pointer ::
     &     signsec(:)

      real(8), external ::
     &     ddot
      


      nsec = opti_info%nsec(iopt)
      if(ntest.gt.10)then
         write(lulog,*) "fixing signs on", ff_vec%name
      end if 
c dbg
c      print *,'nsec = ',nsec
c dbg
      if (nsec.le.1) return
      if (nsec.eq.1.and.opti_info%signsec(1).eq.1d0) return

      ioff = sum(opti_info%nsec(1:iopt))-nsec
      nwfpsec => opti_info%nwfpsec(ioff+1:ioff+nsec)
      idstsec => opti_info%idstsec(ioff+1:ioff+nsec)
      signsec => opti_info%signsec(ioff+1:ioff+nsec)
      call vec_from_da(ff_vec,irecvec,xbuf1,nwfpar)
      do isec = 1, nsec
        xbuf1(idstsec(isec):idstsec(isec)+nwfpsec(isec)-1) = 
     &       signsec(isec)*
     &       xbuf1(idstsec(isec):idstsec(isec)+nwfpsec(isec)-1)
      end do
      call vec_to_da(ff_vec,irecvec,xbuf1,nwfpar)

      return
      end

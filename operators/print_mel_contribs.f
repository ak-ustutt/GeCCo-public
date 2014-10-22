      subroutine print_mel_contribs(luout,mel,fmt,xvals,nblk_in)

      implicit none

      include "mdef_me_list.h"

      integer, intent(in) ::
     &    luout, nblk_in
      type(me_list), intent(in) ::
     &    mel
      real(8), intent(in) ::
     &    xvals(nblk_in)
      character(len=*), intent(in) ::
     &    fmt

      integer, parameter ::
     &    len_descr = 20
     
      character(len=len_descr) ::
     &    descr
      integer, pointer ::
     &    occ(:,:,:), length(:)
      integer ::
     &    nblk, iblk, idxblk, nj

      nblk = mel%op%n_occ_cls
      nj   = mel%op%njoined

      occ => mel%op%ihpvca_occ(:,:,:)

      do iblk = 1, nblk
        idxblk = (iblk-1)*nj+1
        call occ2descr(descr,len_descr,occ(1,1,idxblk),nj)
        if (.not.mel%op%formal_blk(iblk)) then
          write(luout,'(5x,i3,2x,a,'//trim(fmt)//')') 
     &        iblk,descr,xvals(iblk)
        end if
      end do

      end 


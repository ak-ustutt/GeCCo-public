      subroutine print_mel_info(luout,mel)

      implicit none

      include "mdef_me_list.h"

      integer, intent(in) ::
     &    luout
      type(me_list), intent(in) ::
     &    mel

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
      length => mel%len_op_occ(:)

      write(luout,'(2x,47("-"))')
      write(luout,'(3x,"Operator:  ",a)') trim(mel%op%name)
      write(luout,'(3x,"List:      ",a)') trim(mel%label)
      if(associated(mel%fhand))then
        write(luout,'(3x,"Active_record: ",I3)')mel%fhand%current_record
      else
        write(luout,'(3x,"No file associated")')
      end if
      write(luout,'(2x,47("-"))')
      do iblk = 1, nblk
        idxblk = (iblk-1)*nj+1
        call occ2descr(descr,len_descr,occ(1,1,idxblk),nj)
        if (mel%op%formal_blk(iblk)) then
          write(luout,'(5x,i3,2x,a," - formal -")') iblk,descr
        else
          write(luout,'(5x,i3,2x,a,i14)') iblk,descr,length(iblk)
        end if
      end do
      write(luout,'(2x,24("- "))')
      write(luout,'(5x,"total number of elements ",i14)') mel%len_op
      write(luout,'(2x,47("-"))')

      end 

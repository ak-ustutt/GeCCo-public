      subroutine print_mel_info(luout,mel)

      implicit none

      include "mdef_me_list.h"
      include "molpro_out.h"

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
      character(len=2) ::
     &    ci_name

      nblk = mel%op%n_occ_cls
      nj   = mel%op%njoined

      occ => mel%op%ihpvca_occ(:,:,:)
      length => mel%len_op_occ(:)

      if (lmol .and. trim(mel%op%name)=='T2g' .or.
     &    trim(mel%op%name)=='T') then
      ! Special molpro output for the cluster operators
      write(luout,*)
      write(luout,*) "Summery of cluster operator dimensions: "
      write(luout,*)
      write(luout,*) "   Name      Excitation        Size"
      do iblk = 1, nblk
        idxblk = (iblk-1)*nj+1
        call occ2descr(descr,len_descr,occ(1,1,idxblk),nj)

        select case (trim(descr))
          case('V,H')
             ci_name = 'I1'
          case('VV,HV')
             ci_name = 'I1'
          case('VV,HH')
             ci_name = 'I2'
          case('P,V')
             ci_name = 'S0'
          case('PV,VV')
             ci_name = 'S0'
          case('P,H')
             ci_name = 'S1'
          case('PV,HV')
             ci_name = 'S1'
          case('PV,HH')
             ci_name = 'S2'
          case('PP,VV')
             ci_name = 'P0'
          case('PP,HV')
             ci_name = 'P1'
          case('PP,HH')
             ci_name = 'P2'
          case default
             ci_name = 'XX'
        end select

        write(luout,'(4x,a2,11x,a5,i14)')
     &         ci_name,descr,length(iblk)
      end do
      write(luout,*)
      write(luout,'(" Total number of elements: ",i9)') mel%len_op
      write(luout,*)

      else if (lmol .and. trim(mel%op%name)=='T1') then
      ! Do nothing for the T1 operator
      else
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
      end if

      end

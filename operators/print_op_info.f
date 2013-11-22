*----------------------------------------------------------------------*
      subroutine print_op_info(lulog,modestr,op_info)
*----------------------------------------------------------------------*

      implicit none

      include 'mdef_operator_info.h'

      integer, intent(in) ::
     &     lulog
      character(*), intent(in) ::
     &     modestr
      type(operator_info), intent(in) ::
     &     op_info

      integer ::
     &     idx
      type(operator_array), pointer ::
     &     op_arr(:)
      type(me_list_array), pointer ::
     &     mel_arr(:)

      select case(trim(modestr))
      case('op','ops','operator','operators')
        op_arr => op_info%op_arr
        write(lulog,*) 'Number of operators defined: ',op_info%nops
        write(lulog,'(x,40("-"))')
        write(lulog,*)
     &  ' idx     op  vtx blk   current list'
        write(lulog,'(x,40("-"))')
        do idx = 1, op_info%nops
          write(lulog,'(2x,i3,a8,3x,i1,2x,i2,2x,a16)')
     &         idx,trim(op_arr(idx)%op%name),
     &         op_arr(idx)%op%njoined, op_arr(idx)%op%n_occ_cls,
     &         trim(op_arr(idx)%op%assoc_list)          
        end do
        write(lulog,'(x,40("-"))')
      case('mel','me_list','lists','ME-lists')
        mel_arr => op_info%mel_arr
        write(lulog,*) 'Number of lists defined: ',op_info%nmels
        write(lulog,'(x,78("-"))')
        write(lulog,*)
     &  ' idx           list  sym spn      length   op.      '//
     &       '       file'     
        write(lulog,'(x,78("-"))')
        do idx = 1, op_info%nmels
          if (associated(mel_arr(idx)%mel%fhand)) then
            write(lulog,'(2x,i3,a16,2x,i1,x,i2,x,i12,x,a8,x,a29)')
     &           idx,trim(mel_arr(idx)%mel%label),
     &           mel_arr(idx)%mel%gamt,
     &           mel_arr(idx)%mel%mst,
     &           mel_arr(idx)%mel%len_op,
     &           trim(mel_arr(idx)%mel%op%name),
     &           trim(mel_arr(idx)%mel%fhand%name)
          else
            write(lulog,'(2x,i3,a16,2x,i1,x,i2,x,i12,a8)')
     &           idx,trim(mel_arr(idx)%mel%label),
     &           mel_arr(idx)%mel%gamt,
     &           mel_arr(idx)%mel%mst,
     &           mel_arr(idx)%mel%len_op,
     &           trim(mel_arr(idx)%mel%op%name)
          end if
        end do
        write(lulog,'(x,78("-"))')
      end select

      return
      end

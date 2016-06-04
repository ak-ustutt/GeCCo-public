*----------------------------------------------------------------------*
      integer function idx_mel_list(mel_label,op_info)
*----------------------------------------------------------------------*
*     given an me_list name, search op_info and return index of
*     of corresponding me_list
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'mdef_operator_info.h'

      integer, parameter ::
     &     ntest = 00

      character, intent(in) ::
     &     mel_label*(*)
      type(operator_info), intent(in) ::
     &     op_info

      integer ::
     &     imel, llabel, llabel2

      if (ntest.ge.100) then
        write(lulog,*) '--------------------'
        write(lulog,*) 'this is idx_mel_list'
        write(lulog,*) '--------------------'
        write(lulog,*) ' looking for: "',trim(mel_label),'"'
      end if

      ! note: currently this routine is called quite often
      ! via update_op_arr ... try to be more efficient and avoid trim
      idx_mel_list = -1
      llabel = len_trim(mel_label)
      if (llabel.eq.0) return
      do imel = 1, op_info%nmels
        if (mel_label(1:llabel).eq.
     &      op_info%mel_arr(imel)%mel%label(1:llabel)) then
          if (op_info%mel_arr(imel)%mel%label(llabel+1:llabel+1).eq.' ')
     &    then
            idx_mel_list = imel
            exit
          end if
        end if
      end do

      if (ntest.ge.100) then
        write(lulog,*) 'result: ',idx_mel_list
      end if

      return
      end

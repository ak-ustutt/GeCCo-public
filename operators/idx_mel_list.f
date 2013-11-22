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
     &     imel

      if (ntest.ge.100) then
        write(lulog,*) '--------------------'
        write(lulog,*) 'this is idx_mel_list'
        write(lulog,*) '--------------------'
        write(lulog,*) ' looking for: "',trim(mel_label),'"'
      end if

      idx_mel_list = -1
      do imel = 1, op_info%nmels
        if (trim(mel_label).eq.
     &      trim(op_info%mel_arr(imel)%mel%label)) then
          idx_mel_list = imel
          exit
        end if
      end do

      if (ntest.ge.100) then
        write(lulog,*) 'result: ',idx_mel_list
      end if

      return
      end

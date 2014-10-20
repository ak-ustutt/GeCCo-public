*----------------------------------------------------------------------*
      function get_mel_record(mel)
*----------------------------------------------------------------------*
*     return the current record of the mel
*----------------------------------------------------------------------*
      implicit none

      include 'mdef_operator_info.h'

      type(me_list), intent(inout) ::
     &     mel
      
      integer :: 
     &     get_mel_record

      get_mel_record = mel%fhand%current_record

      return
      end

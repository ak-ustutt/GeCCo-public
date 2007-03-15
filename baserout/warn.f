      subroutine warn(rout,str)

      implicit none

      include 'stdunit.h'

      character, intent(in) ::
     &     rout*(*),str*(*)

      write(luout,'(/x,"WARNING IN <",a,">: ",a/)')
     &       rout,str

      return
      end

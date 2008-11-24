      subroutine warn(rout,str)

      implicit none

      include 'stdunit.h'
      include 'warnings.h'

      character, intent(in) ::
     &     rout*(*),str*(*)
      
      write(luwarn,'(/x,"WARNING IN <",a,">: ",a/)')
     &       rout,str
      write(luout,'(/x,"WARNING IN <",a,">: ",a/)')
     &       rout,str

      return
      end

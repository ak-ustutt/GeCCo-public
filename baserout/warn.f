      subroutine warn(rout,str)

      implicit none

      include 'stdunit.h'
      include 'warnings.h'

      character, intent(in) ::
     &     rout*(*),str*(*)
      
      nwarn = nwarn + 1
      write(luwarn,'(/x,"WARNING IN <",a,">: ",a/)')
     &       rout,str
      write(lulog,'(/x,"WARNING IN <",a,">: ",a/)')
     &       rout,str

      return
      end

      subroutine quit(level,rout,str)

      implicit none

      include 'stdunit.h'

      integer, intent(in) ::
     &     level
      character, intent(in) ::
     &     rout*(*),str*(*)

      character(24) ::
     &     date
      character(256) ::
     &     host

      call hostname(host)
      call datum(date)

      if (level.eq.0) then
        write(luout,'(/x,"ERROR IN <",a,">: ",a/)') rout,str
        write(luout,'(x,"run ends at ",a,"   host: ",a)')
     &     trim(date),trim(host)
        stop 'error exit'
      else
        write(luout,'(/x,"INTERNAL ERROR IN <",a,">: ",a/)') rout,str
        write(luout,'(x,"run ends at ",a,"   host: ",a)')
     &     trim(date),trim(host)
c dbg
c        call tracebackqq()
c dbg
        stop 'internal error exit'
      end if

      return
      end

      subroutine quit_mem(need,free,rout,str)

      implicit none

      include 'stdunit.h'

      integer, intent(in) ::
     &     need,free
      character, intent(in) ::
     &     rout*(*),str*(*)

      character(24) ::
     &     date
      character(256) ::
     &     host

      call hostname(host)
      call datum(date)

      write(luout,'(/x,"needed: ",i12," R*8 words  ",e8.1," Mb")')
     &       need,dble(need)/128d0/1024d0
      write(luout,'(/x,"free:   ",i12," R*8 words  ",e8.1," Mb")')
     &       free,dble(free)/128d0/1024d0
      write(luout,'(/x,"INSUFFICIENT MEMORY IN <",a,">: ",a/)')
     &       rout,str

      write(luout,'(x,"run ends at ",a,"   host: ",a)')
     &     trim(date),trim(host)
      stop 'insufficient memory'

      return
      end

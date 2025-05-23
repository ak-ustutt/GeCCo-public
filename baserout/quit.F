      subroutine quit(level,rout,str)

#if defined(__INTEL_COMPILER)
      use ifcore, only : tracebackqq
#endif // if __INTEL_COMPILER     
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
        write(lulog,'(/x,"ERROR IN <",a,">: ",a/)') rout,str
        write(lulog,'(x,"run ends at ",a,"   host: ",a)')
     &     trim(date),trim(host)
        if (lulog.ne.luout) then
          write(luout,'(/x,"ERROR IN <",a,">: ",a/)') rout,str
          write(luout,'(x,"run ends at ",a,"   host: ",a)')
     &     trim(date),trim(host)
       end if
#if defined(__INTEL_COMPILER)
       call tracebackqq(user_exit_code = -1)
#endif // if __INTEL_COMPILER
        stop 'error exit'
      else
        write(lulog,'(/x,"INTERNAL ERROR IN <",a,">: ",a/)') rout,str
        write(lulog,'(x,"run ends at ",a,"   host: ",a)')
     &     trim(date),trim(host)
        if (luout.ne.lulog) then
          write(luout,'(/x,"INTERNAL ERROR IN <",a,">: ",a/)') rout,str
          write(luout,'(x,"run ends at ",a,"   host: ",a)')
     &     trim(date),trim(host)
        end if
c dbg
#if defined(__INTEL_COMPILER)
        call tracebackqq(user_exit_code = -1)
#endif // if __INTEL_COMPILER
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

      write(lulog,'(/x,"needed: ",i12," R*8 words  ",e8.1," Mb")')
     &       need,dble(need)/128d0/1024d0
      write(lulog,'(/x,"free:   ",i12," R*8 words  ",e8.1," Mb")')
     &       free,dble(free)/128d0/1024d0
      write(lulog,'(/x,"INSUFFICIENT MEMORY IN <",a,">: ",a/)')
     &       rout,str
      write(lulog,'(x,"run ends at ",a,"   host: ",a)')
     &     trim(date),trim(host)

      if (lulog.ne.luout) then
        write(luout,'(/x,"needed: ",i12," R*8 words  ",e8.1," Mb")')
     &       need,dble(need)/128d0/1024d0
        write(luout,'(/x,"free:   ",i12," R*8 words  ",e8.1," Mb")')
     &       free,dble(free)/128d0/1024d0
        write(luout,'(/x,"INSUFFICIENT MEMORY IN <",a,">: ",a/)')
     &       rout,str
        write(luout,'(x,"run ends at ",a,"   host: ",a)')
     &     trim(date),trim(host)

      end if

      stop 'insufficient memory'

      return
      end

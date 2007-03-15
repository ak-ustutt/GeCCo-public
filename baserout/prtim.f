      subroutine prtim(luout,str,cpu,sys,wall)

      implicit none

      integer, intent(in) ::
     &     luout
      character, intent(in) ::
     &     str*(*)
      real*8, intent(in) ::
     &     cpu, sys, wall
      integer, parameter ::
     &     mxlen = 34
      character ::
     &     str_scr*(mxlen)
      integer ::
     &     len

      len = min(mxlen,len_trim(str))
      str_scr(1:len) = str(1:len)
      if (len.lt.mxlen)
     &     str_scr(len+1:len+mxlen-len) = ' '
      write(luout,
     &     '(x,"@ ",a,"  cpu/sys/wall: ",f11.2,"/",f11.2,"/",'//
     &     'f11.2" s")') 
     &     str_scr,cpu,sys,wall

      return
      end

*----------------------------------------------------------------------*
      subroutine wrt_occ_rstr(luout,idx,iocc,irstr,ngas)
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'

      integer, intent(in) ::
     &     luout, iocc(ngastp,2), ngas, irstr(2,ngas,2,2), idx
      character ::
     &     fmtr*8, fmto*3
      integer, external ::
     &     ielsqsum

      write(fmto,'(i1,"i3")') ngastp
      write(fmtr,'(i1,"(x,2i3)")') ngas 
      write(luout,'(x,i4,2x,"/",'//fmto//',"'//char(92)//'",'//
     &             fmtr//',3x,'//fmtr//')')
     &     idx,iocc(1:ngastp,1),
     &     irstr(1:2,1:ngas,1,1),
     &     irstr(1:2,1:ngas,1,2)
      write(luout,'(7x,"'//
     &             char(92)//'",'//fmto//',"/",'//
     &             fmtr//',3x,'//fmtr//')')
     &     iocc(1:ngastp,2),
     &     irstr(1:2,1:ngas,2,1),
     &     irstr(1:2,1:ngas,2,2)

      return
      end

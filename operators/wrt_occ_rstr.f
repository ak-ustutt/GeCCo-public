*----------------------------------------------------------------------*
      subroutine wrt_occ_rstr(luout,idx,iocc,irstr,ngas,nspin)
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'

      integer, intent(in) ::
     &     luout, iocc(ngastp,2), ngas, nspin,
     &     irstr(2,ngas,2,2,nspin), idx
      character ::
     &     fmtr*8, fmto*3, fmtx*5
      integer, external ::
     &     ielsqsum

      write(fmto,'(i1,"i3")') ngastp
      write(fmtx,'(i1,"(3x)")')  ngastp
      write(fmtr,'(i1,"(x,2i3)")') ngas 
      write(luout,'(x,i4,2x,"/",'//fmto//',"'//char(92)//'",'//
     &             fmtr//',3x,'//fmtr//')')
     &     idx,iocc(1:ngastp,1),
     &     irstr(1:2,1:ngas,1,1,1),
     &     irstr(1:2,1:ngas,1,2,1)
      if (nspin.eq.2) then
        write(luout,'(x,i4,2x,"|",'//fmtx//',"|",'//
     &             fmtr//',3x,'//fmtr//')')
     &     irstr(1:2,1:ngas,1,1,2),
     &     irstr(1:2,1:ngas,1,2,2)
      end if
      write(luout,'(7x,"'//
     &             char(92)//'",'//fmto//',"/",'//
     &             fmtr//',3x,'//fmtr//')')
     &     iocc(1:ngastp,2),
     &     irstr(1:2,1:ngas,2,1,1),
     &     irstr(1:2,1:ngas,2,2,1)
      if (nspin.eq.2) then
        write(luout,'(x,i4,2x," ",'//fmtx//'," ",'//
     &             fmtr//',3x,'//fmtr//')')
     &     irstr(1:2,1:ngas,2,1,2),
     &     irstr(1:2,1:ngas,2,2,2)
      end if

      return
      end

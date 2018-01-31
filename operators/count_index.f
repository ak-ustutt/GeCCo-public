*----------------------------------------------------------------------*
      subroutine count_index(lulog,idx,iocc,irstr,ngas,nspin,nops)
*----------------------------------------------------------------------*
!     Count index of operator
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'

      integer, intent(in) ::
     &     lulog, iocc(ngastp,2), ngas, nspin,
     &     irstr(2,ngas,2,2,nspin), idx
      integer, intent(inout) ::
     &     nops(4,2)                     ! Matrix of index info
      character ::
     &     fmtr*8, fmto*3, fmtx*5
      integer, external ::
     &     ielsqsum
      integer ::
     &     pc=0, hc=0, vc=0, xc=0,       ! Number of creation and annhilation operators
     &     pa=0, ha=0, va=0, xa=0

!      write(fmto,'(i1,"i3")') ngastp
!      write(fmtx,'(i1,"(3x)")')  ngastp
!      write(fmtr,'(i1,"(x,2i3)")') ngas 
!      write(lulog,'(x,i4,2x,"/",'//fmto//',"'//char(92)//'",'//
!     &             fmtr//',3x,'//fmtr//')')
!     &     idx,iocc(1:ngastp,1),
!     &     irstr(1:2,1:ngas,1,1,1),
!     &     irstr(1:2,1:ngas,1,2,1)
!      if (nspin.eq.2) then
!        write(lulog,'(x,i4,2x,"|",'//fmtx//',"|",'//
!     &             fmtr//',3x,'//fmtr//')')
!     &     irstr(1:2,1:ngas,1,1,2),
!     &     irstr(1:2,1:ngas,1,2,2)
!      end if
!      write(lulog,'(7x,"'//
!     &             char(92)//'",'//fmto//',"/",'//
!     &             fmtr//',3x,'//fmtr//')')
!     &     iocc(1:ngastp,2),
!     &     irstr(1:2,1:ngas,2,1,1),
!     &     irstr(1:2,1:ngas,2,2,1)
!      if (nspin.eq.2) then
!        write(lulog,'(x,i4,2x," ",'//fmtx//'," ",'//
!     &             fmtr//',3x,'//fmtr//')')
!     &     irstr(1:2,1:ngas,2,1,2),
!     &     irstr(1:2,1:ngas,2,2,2)
!      end if

      write(lulog, *) 'Count index', iocc(1:ngastp,1)
      write(lulog, *) 'Count index', iocc(1:ngastp,2)

      nops(1,1)=iocc(1,1)
      nops(2,1)=iocc(2,1)
      nops(3,1)=iocc(3,1)
      nops(4,1)=iocc(4,1)
      nops(1,2)=iocc(1,2)
      nops(2,2)=iocc(2,2)
      nops(3,2)=iocc(3,2)
      nops(4,2)=iocc(4,2)

      write(lulog, *) 'hc:', nops(1,1)
      write(lulog, *) 'pc:', nops(2,1)
      write(lulog, *) 'vc:', nops(3,1)
      write(lulog, *) 'xc:', nops(4,1)
      write(lulog, *) 'ha:', nops(1,2)
      write(lulog, *) 'pa:', nops(2,2)
      write(lulog, *) 'va:', nops(3,2)
      write(lulog, *) 'xa:', nops(4,2)

! return an array with this info

      return
      end



      subroutine index_array(lulog,tensor,i_array)

      include 'mdef_operator_info.h'

      integer, intent(in) ::
     &     lulog
      type(operator), intent(in) ::
     &     tensor
      character, intent(inout)  ::
     &     i_array(2,2)
      integer ::
     &     idx, iblk,
     &     nops(4,2)     ! Matrix of index info
      character, dimension(4) ::
     &     hol=(/ 'i','j','k','l' /),
     &     par=(/ 'a','b','c','d' /),
     &     val=(/ 'u','v','w','x' /)


       do idx = 1, tensor%n_occ_cls*tensor%njoined
         iblk = (idx-1)/tensor%njoined + 1
         call count_index(lulog,iblk,
     &       tensor%ihpvca_occ(1,1,idx),
     &       tensor%igasca_restr(1,1,1,1,1,idx),
     &       tensor%ngas,tensor%nspin,nops)
       end do

       ! Annhilators first
       do i=1, nops(1,2)
           i_array(1,i)=hol(i)
       end do
       do i=1, nops(2,2)
         i_array(1,i+nops(1,2))=par(i)
       end do
       do i=1, nops(3,2)
         i_array(1,i+nops(1,2)+nops(2,2))=val(i)
       end do
       ! Creations second
       do i=1, nops(1,1)
         i_array(2,i)=hol(i+nops(1,2))
       end do
       do i=1, nops(2,1)
         i_array(2,i+nops(1,1))=par(i+nops(2,2))
       end do
       do i=1, nops(3,1)
         i_array(2,i+nops(1,1)+nops(2,1))=val(i+nops(3,2))
       end do

       return
       end


*----------------------------------------------------------------------*
      subroutine count_index2(lulog,idx,iocc,irstr,ngas,nspin,nops)
*----------------------------------------------------------------------*
!     Count index of operator
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'

      integer, intent(in) ::
     &     lulog, iocc(ngastp,2), ngas, nspin,
     &     irstr(2,ngas,2,2,nspin), idx
      integer, intent(inout) ::
     &     nops(4,2)                     ! Matrix of index info
      character ::
     &     fmtr*8, fmto*3, fmtx*5
      integer, external ::
     &     ielsqsum
      integer ::
     &     pc=0, hc=0, vc=0, xc=0,       ! Number of creation and annhilation operators
     &     pa=0, ha=0, va=0, xa=0

      write(lulog, *) 'Count index', iocc(1:ngastp,1)
      write(lulog, *) 'Count index', iocc(1:ngastp,2)

      nops(1,1)=iocc(1,1)
      nops(2,1)=iocc(2,1)
      nops(3,1)=iocc(3,1)
      nops(4,1)=iocc(4,1)
      nops(1,2)=iocc(1,2)
      nops(2,2)=iocc(2,2)
      nops(3,2)=iocc(3,2)
      nops(4,2)=iocc(4,2)

      write(lulog, *) 'hc:', nops(1,1)
      write(lulog, *) 'pc:', nops(2,1)
      write(lulog, *) 'vc:', nops(3,1)
      write(lulog, *) 'xc:', nops(4,1)
      write(lulog, *) 'ha:', nops(1,2)
      write(lulog, *) 'pa:', nops(2,2)
      write(lulog, *) 'va:', nops(3,2)
      write(lulog, *) 'xa:', nops(4,2)

      return
      end



      subroutine index_array2(lulog,tensor,i_array,j_array)

      include 'opdim.h'
      include 'def_contraction.h'

      integer, intent(in) ::
     &     lulog
      type(binary_contr), intent(in) ::
     &     tensor
      character, intent(inout)  ::
     &     i_array(2,2),
     &     j_array(2,2)
      integer ::
     &     idx, iblk,
     &     nops(4,2)     ! Matrix of index info
      character, dimension(4) ::
     &     hol=(/ 'i','j','k','l' /),
     &     par=(/ 'a','b','c','d' /),
     &     val=(/ 'u','v','w','x' /)


       do i = 1, tensor%nj_op1
        call count_index2(lulog,i,
     &       tensor%occ_op1(1:,1:,i),
     &       tensor%rst_op1(1:,1:,1:,1:,1:,i),tensor%ngas,
     &       tensor%nspin,nops)
      end do
       

       ! Annhilators first
       do i=1, nops(1,2)
           i_array(1,i)=hol(i)
       end do
       do i=1, nops(2,2)
         i_array(1,i+nops(1,2))=par(i)
       end do
       do i=1, nops(3,2)
         i_array(1,i+nops(1,2)+nops(2,2))=val(i)
       end do
       ! Creations second
       do i=1, nops(1,1)
         i_array(2,i)=hol(i+nops(1,2))
       end do
       do i=1, nops(2,1)
         i_array(2,i+nops(1,1))=par(i+nops(2,2))
       end do
       do i=1, nops(3,1)
         i_array(2,i+nops(1,1)+nops(2,1))=val(i+nops(3,2))
       end do
        
        do i = 1, tensor%nj_op1
        call count_index2(lulog,i,
     &       tensor%occ_op2(1:,1:,i),
     &       tensor%rst_op2(1:,1:,1:,1:,1:,i),tensor%ngas,
     &       tensor%nspin,nops)
      end do
       
       ! Annhilators first
       do i=1, nops(1,2)
           j_array(1,i)=hol(i)
       end do
       do i=1, nops(2,2)
         j_array(1,i+nops(1,2))=par(i)
       end do
       do i=1, nops(3,2)
         j_array(1,i+nops(1,2)+nops(2,2))=val(i)
       end do
       ! Creations second
       do i=1, nops(1,1)
         j_array(2,i)=hol(i+nops(1,2))
       end do
       do i=1, nops(2,1)
         j_array(2,i+nops(1,1))=par(i+nops(2,2))
       end do
       do i=1, nops(3,1)
         j_array(2,i+nops(1,1)+nops(2,1))=val(i+nops(3,2))
       end do

       return
       end

      subroutine import_r12_dalton(hop,ffr12,ffinp,
     &     mode,str_info,orb_info)
*-----------------------------------------------------------------------
*     Routine to read in and reorder 2-electron integrals required for 
*     R12 calculations. Integrals are reordered to be in type order.
*     GWR July 2007.
*-----------------------------------------------------------------------
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_operator.h'
      include 'hpvxseq.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'def_filinf.h'
      include 'par_dalton.h'
      include 'ifc_memman.h'

      integer, parameter::
     &     ntest=00

      type(orbinf),intent(in),target ::
     &     orb_info
      type(strinf),intent(in) ::
     &     str_info
      type(operator),intent(in) ::
     &     hop
      type(filinf),intent(inout) ::
     &     ffr12
      character*256, intent(in) ::
     &     ffinp
      integer, intent(in) ::
     &     mode

      integer ::
     &     idx,jdx,kdx,igas,isym,i,j,nsym,ngas,ntoob,caborb,lu2in,
     &     ip,iq,ir,is,nstr,istr,ifree,lbuf
      integer ::
     &     index(4),igam(4),idss(4),igtp(4),idxstr(8)
      real(8) ::
     &     int
      logical ::
     &     fexist,closeit,ierr
      type(filinf) ::
     &     ffmor

      integer, pointer ::
     &     ihpvgas(:),igamorb(:),igasorb(:),idx_gas(:),iad_gas(:),
     &     reost(:)
      integer, allocatable ::
     &     tosym(:),totyp(:),koffs(:),reord(:)
      real*8,pointer ::
     &     r12scr(:)

      ifree=mem_setmark('import_r12')

      lbuf=hop%len_op
      if(lbuf.gt.ifree)
     &     call quit(1,'import_r12','insufficient memory')
      ifree=mem_alloc_real(r12scr,lbuf,'r12buf')

      r12scr(1:lbuf)=0d0

      nsym=orb_info%nsym
      ngas=orb_info%ngas
      ntoob=orb_info%ntoob
      caborb=orb_info%caborb

      allocate(tosym(ntoob+caborb),totyp(ntoob+caborb),koffs(nsym),
     &     reord(ntoob+caborb))

      ! If an R12 calculation is requested then must do an initial 
      ! reordering, as the external orbitals are extracted separately
      ! from the others. All other orbitals are written together in 
      ! symmetry ordering (i.e. h/p/v in one list). It is necessary to
      ! combine this list with that of the x-space (also symmetry 
      ! ordered) before subsequent arrays can act on the total set.
      reost=>orb_info%ireost
      do j=1,ntoob
        reord(j)=reost(j)
      enddo
      do j=ntoob+1,ntoob+caborb
        reord(j)=j
      enddo  

      ! Open the required MO-integral file.
      inquire(file=ffinp,exist=fexist)
      if(.not.fexist)
     &     call quit(1,'import_r12_dalton','No MO integral file.')

c      call file_init(ffmor,mo_r,ftyp_sq_frm,0)
c      lu2in=ffmor%unit
c      call file_open(ffmor)
c      write(luout,*)lu2in
      lu2in=99
      open(unit=lu2in,file=ffinp,status='old',form='formatted',
     &     access='sequential')
      rewind(lu2in)
      read(lu2in,*)
      read(lu2in,*)

      ! Point to specific parts of orb_info.
      ihpvgas=>orb_info%ihpvgas
      igamorb=>orb_info%igamorb
      igasorb=>orb_info%igasorb
      idx_gas=>orb_info%idx_gas
      iad_gas=>orb_info%iad_gas

      ! Open output file.
      if(ffr12%unit.le.0)then
        call file_open(ffr12)
        closeit = .true.
      else
        closeit = .false.
      endif  

      ! Loop over the written integrals.
      do 
        ! Read in the integral and its indices. Replace the indices with 
        ! those due to the type ordering array.
        read(unit=lu2in,fmt=1000,end=999,err=998)ip,iq,ir,is,int
        if((ip.lt.iq.and.iq.le.ntoob).or.ir.lt.is)then
          cycle
        endif  

        ! Do two loops as we must "manually" permute p <=> q.
        do i=1,2
          if(i.eq.1)then
            index(1)=reord(ip)
            index(2)=reord(ir)
            index(3)=reord(iq)
            index(4)=reord(is)
          elseif(i.eq.2)then
            if(ir.gt.ntoob)cycle
            if(ip.eq.iq.or.ir.eq.is)cycle
            index(1)=reord(ip)
            index(2)=reord(is)
            index(3)=reord(iq)
            index(4)=reord(ir)
c          elseif(i.eq.3)then
c            if(ip.gt.ntoob.or.iq.gt.ntoob)cycle
c            if(ip.eq.iq.or.ir.eq.is)cycle
c            index(1)=reord(iq)
c            index(2)=reord(ir)
c            index(3)=reord(ip)
c            index(4)=reord(is)
c          else
c            if(ip.gt.ntoob.or.iq.gt.ntoob)cycle
c            if(ip.eq.iq.or.ir.eq.is)cycle
c            index(1)=reord(iq)
c            index(2)=reord(is)
c            index(3)=reord(ip)
c            index(4)=reord(ir)
          endif  
c dbg
c          if (ip.eq.9.and.iq.eq.4.and.ir.eq.8.and.is.eq.2)
c     $         print *,'1 da isset ',int
c          if ((ip.eq.9.and.iq.eq.7.and.ir.eq.6.and.is.eq.5).or.
c     &         (ip.eq.7.and.iq.eq.6.and.ir.eq.9.and.is.eq.5).or.
c     &         (ip.eq.9.and.iq.eq.6.and.ir.eq.7.and.is.eq.5))then
c            print *,'2 da isset ',int
c            write(luout,'(a,4i5)')'raw indices = ',ip,iq,ir,is
c            write(luout,'(a,4i5)')'reo indices = ',index(1:4)
c          endif
c dbg            

          ! Place into arrays the irrep. and active space indices of the 
          ! orbitals in the integral under consideration.

          do j=1,4
            igam(j)=igamorb(index(j))
            idss(j)=igasorb(index(j))
            igtp(j)=ihpvgas(idss(j))            
          enddo

          ! Loop if not all orbitals are active.
          if(iad_gas(idss(1)).ne.2.or.iad_gas(idss(2)).ne.2.or.
     &         iad_gas(idss(3)).ne.2.or.iad_gas(idss(4)).ne.2) cycle

          do j=1,4
            idss(j)=idss(j)-idx_gas(igtp(j))+1
          enddo
  
          ! Generate the string addresses of the current integral by 
          ! reference to the spin-orbital basis.
          if(mode.eq.0)then
            call idx42str(nstr,idxstr,
     &           index,igam,idss,igtp,
     &           orb_info,str_info,hop,hpvxseq,ierr)
          else
            call idx42str2(nstr,idxstr,mode,
     &           index,igam,idss,igtp,
     &           orb_info,str_info,hop,hpvxseq,ierr)
          endif
  
          ! Store integral in r12scr.
          do istr=1,nstr
c dbg
c            if (abs(idxstr(istr)).lt.1.or.abs(idxstr(istr)).gt.lbuf)then
c              write(luout,*) idxstr(istr),lbuf,nstr
c              stop 'range'
c            end if
c dbg

            ! Add the integrals to the string array.
            if(.not.ierr)then
              if(idxstr(istr).gt.0)
     &             r12scr(idxstr(istr))=
     &             r12scr(idxstr(istr))+int
              if(idxstr(istr).lt.0)
     &             r12scr(-idxstr(istr))=
     &             r12scr(-idxstr(istr))-int
            endif  
          enddo  
        enddo
  
      enddo  

 998  call quit(1,'import_r12_dalton','Error reading integrals.')
 999  continue

      call put_vec(ffr12,r12scr,1,lbuf)

      if(closeit)
     &     call file_close_keep(ffr12)
c      call file_close_keep(ffmor)
      close(unit=lu2in,status='keep')
      deallocate(tosym,totyp,koffs,reord)

      if(ntest.ge.1000)then
        write(luout,*)'R12-integral import results:'
        call wrt_op_file(luout,5,ffr12,hop,
     &       1,hop%n_occ_cls,
     &       str_info,orb_info)
      endif  

      ifree=mem_flushmark()

      return
 1000 format(4i5,e25.15)
      end

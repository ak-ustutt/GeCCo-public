      subroutine import_r12_dalton(oplist,fname_inp,
     &     mode,str_info,orb_info)
*-----------------------------------------------------------------------
*     Routine to read in and reorder 2-electron integrals required for 
*     R12 calculations. Integrals are reordered to be in type order.
*     GWR July 2007.
*-----------------------------------------------------------------------
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'hpvxseq.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'par_dalton.h'
      include 'ifc_memman.h'

      integer, parameter::
     &     ntest= 100

      type(orbinf),intent(in),target ::
     &     orb_info
      type(strinf),intent(in) ::
     &     str_info
      type(me_list),intent(inout) ::
     &     oplist
      character*(*), intent(in) ::
     &     fname_inp
      integer, intent(in) ::
     &     mode

      integer ::
     &     idx,jdx,kdx,igas,isym,i,j,nsym,ngas,ntoob,caborb,lu2in,
     &     ip,iq,ir,is,nstr,istr,ifree,lbuf,asym_el
      integer ::
     &     index(4),igam(4),idss(4),igtp(4),idxstr(8)
      real(8) ::
     &     int
      logical ::
     &     fexist,closeit,error,antisym
      type(filinf) ::
     &     ffinp

      integer, pointer ::
     &     ihpvgas(:,:),igamorb(:),igasorb(:),idx_gas(:),iad_gas(:),
     &     reost(:)
      integer, allocatable ::
     &     tosym(:),totyp(:),koffs(:),reord(:)
      real*8,pointer ::
     &     r12scr(:)

      type(operator), pointer ::
     &     op
      type(filinf), pointer ::
     &     ffop
      real(8) ::
     &     cpu0, sys0, wall0, cpu, sys, wall


      call quit(1,'import_r12_dalton','obsolete!')

      ifree=mem_setmark('import_r12')

      op =>   oplist%op
      ffop => oplist%fhand

      if(ntest.ge.100)then
        write(lulog,*)'=================='
        write(lulog,*)'Import-r12-dalton '
        write(lulog,*)'=================='
        write(lulog,*)'Operator-list: ',trim(oplist%label)
      endif

      call atim_csw(cpu0,sys0,wall0)

      lbuf=oplist%len_op
      if(lbuf.gt.ifree)
     &     call quit(1,'import_r12','insufficient memory')
      ifree=mem_alloc_real(r12scr,lbuf,'r12buf')

      r12scr(1:lbuf)=0d0

      nsym=orb_info%nsym
      ngas=orb_info%ngas
      ntoob=orb_info%ntoob
      caborb=orb_info%caborb

c      ! Define some flags for different integral types.
c      ! Mode 3 and 4 deal with (pq|[Ti,r12]|rs) type integrals whilst
c      ! the other modes deal with 'normal' 2-electron integrals. 
c      if(mode.eq.3)then
c        antisym = .true.
c        asym_el = 2
c      elseif(mode.eq.4)then
c        antisym=.true.
c        asym_el = 1
c      else
c        antisym = .false.
c        asym_el = 0
c      endif  

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
      inquire(file=fname_inp,exist=fexist)
      if(.not.fexist)
     &     call quit(1,'import_r12_dalton','No MO integral file: '//
     &       trim(fname_inp))

      call file_init(ffinp,fname_inp,ftyp_sq_frm,0)
      call file_open(ffinp)

      lu2in=ffinp%unit
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
      if(ffop%unit.le.0)then
        call file_open(ffop)
        closeit = .true.
      else
        closeit = .false.
      endif  

      ! Loop over the written integrals.
      int_loop: do 
        ! Read in the integral and its indices. Replace the indices with 
        ! those due to the type ordering array.
c        read(unit=lu2in,fmt=1000,end=999,err=998)ip,iq,ir,is,int
        read(unit=lu2in,fmt=*,end=999,err=998)ip,iq,ir,is,int
        if(mode.ne.2)then
          if((ip.lt.iq.and.iq.le.ntoob).or.ir.lt.is)cycle int_loop    
        endif  

c dbg
c        if(mode.eq.1)then
c          if(ip.eq.3.and.iq.eq.1.and.ir.eq.2.and.is.eq.1)then
c            write(lulog,*)'testing'
c            write(lulog,*)int
c          endif
c        endif
c dbg

        ! Do two loops as we must "manually" permute p <=> q.
        do i=1,2
          if(i.eq.1)then
            index(1)=reord(ip)
            index(2)=reord(ir)
            index(3)=reord(iq)
            index(4)=reord(is)
          elseif(i.eq.2)then
            if(mode.eq.2)cycle int_loop
            if(ir.gt.ntoob)cycle int_loop
            if(ip.eq.iq.or.ir.eq.is)cycle int_loop
            index(1)=reord(ip)
            index(2)=reord(is)
            index(3)=reord(iq)
            index(4)=reord(ir)
          endif  

          ! Place into arrays the irrep. and active space indices of the 
          ! orbitals in the integral under consideration.

          do j=1,4
            igam(j)=igamorb(index(j))
            idss(j)=igasorb(index(j))
            igtp(j)=ihpvgas(idss(j),1) ! OPEN SHELL: adapt            
          enddo

          ! we need some of the inactive orbitals as well ...
c          ! Loop if not all orbitals are active.
c          if(iad_gas(idss(1)).ne.2.or.iad_gas(idss(2)).ne.2.or.
c     &         iad_gas(idss(3)).ne.2.or.iad_gas(idss(4)).ne.2) cycle

          do j=1,4
            idss(j)=idss(j)-idx_gas(igtp(j))+1
          enddo

c dbg
c          if(mode.eq.1)then
c            if(ip.eq.3.and.iq.eq.1.and.ir.eq.2.and.is.eq.1)then
c              write(lulog,*)'testing 2'
c              write(lulog,*)index(1:4)
c              write(lulog,*)int
c            endif
c          endif
c dbg

          ! Generate the string addresses of the current integral by 
          ! reference to the spin-orbital basis.
          call idx42str2(nstr,idxstr,mode,
     &         index,igam,idss,igtp,
     &         orb_info,str_info,oplist,hpvxseq,error)

c dbg
c          if(mode.eq.1)then
c            if(index(1).eq.8.and.index(2).eq.2.and.index(3).eq.2
c     &           .and.index(4).eq.2)then
c              write(lulog,*)'nstr = ',nstr,error
c            endif
c          endif
c dbg

          ! Store integral in r12scr.
          do istr=1,nstr

            ! Add the integrals to the string array.
            if(.not.error)then
              if(idxstr(istr).gt.0)
     &             r12scr(idxstr(istr))=
     &             r12scr(idxstr(istr))+int
              if(idxstr(istr).lt.0)
     &             r12scr(-idxstr(istr))=
     &             r12scr(-idxstr(istr))-int
            endif  

          enddo  
        enddo

      enddo int_loop 

 998  call quit(1,'import_r12_dalton','Error reading integrals.')
 999  continue

      call put_vec(ffop,r12scr,1,lbuf)

      if(closeit)
     &     call file_close_keep(ffop)
      call file_close_keep(ffinp)

      deallocate(tosym,totyp,koffs,reord)

      ifree=mem_flushmark()

      call atim_csw(cpu,sys,wall)

      if (iprlvl.ge.5) 
     &     call prtim(lulog,'time in 2int import',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      return
 1000 format(4i5,e25.15)
      end

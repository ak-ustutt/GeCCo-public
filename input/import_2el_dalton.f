      subroutine import_2el_dalton(oplist,fname_inp,
     &     mode,scaling,anti,str_info,orb_info)
*-----------------------------------------------------------------------
*     Routine to read in and reorder 2-electron integrals required for 
*     R12 calculations. Integrals are reordered to be in type order.
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
     &     ntest= 00

      type(orbinf),intent(in),target ::
     &     orb_info
      type(strinf),intent(in) ::
     &     str_info
      type(me_list),intent(inout) ::
     &     oplist
      character*(*), intent(in) ::
     &     fname_inp
      integer, intent(in) ::
     &     mode, scaling
      logical, intent(in) ::
     &     anti

      real(8), parameter ::
     &     scal_aa = 0.250d0,
     &     scal_ab = 0.375d0

      integer ::
     &     rt, ct, idx_int, idx_typ, ntypes, ifree, nbuffer, igas,
     &     ngam, ngas, caborb, ntoob, naux_max, iaux_max, idx,
     &     nh_min, nh_max
      integer ::
     &     nints(8)
      real(8) ::
     &     fac_s, fac_t
      logical ::
     &     fexist,closeit
      type(filinf) ::
     &     ffinp

      integer, pointer ::
     &     ihpvgas(:,:),igamorb(:),igasorb(:),idx_gas(:),
     &     mnmxspc(:,:)
      real*8,pointer ::
     &     buffer(:)

      type(operator), pointer ::
     &     op
      type(filinf), pointer ::
     &     ffop
      logical ::
     &     hhgg_list
      real(8) ::
     &     cpu0, sys0, wall0, cpu, sys, wall, cpux, sysx 
      integer, external ::
     &     idx_int_graph
      integer, pointer ::
     &     iy_int(:,:,:), typetab(:)

      integer, external ::
     &     max_rank_op

      ifree=mem_setmark('import_r12')

      op =>   oplist%op
      ffop => oplist%fhand

      if(ntest.ge.100)then
        write(luout,*)'=================='
        write(luout,*)'Import-r12-dalton '
        write(luout,*)'=================='
        write(luout,*)'Operator-list: ',trim(oplist%label)
      endif

      ngam  = orb_info%nsym
      ngas  = orb_info%ngas
      ntoob = orb_info%ntoob
      caborb = orb_info%caborb
      igamorb => orb_info%igamorb
      ihpvgas => orb_info%ihpvgas

      call atim_csw(cpu0,sys0,wall0)

      ! Open the required MO-integral file.
      inquire(file=fname_inp,exist=fexist)
      if(.not.fexist)
     &     call quit(1,'import_r12_dalton','No MO integral file: '//
     &       trim(fname_inp))

      if (mode.gt.0) then
        call file_init(ffinp,fname_inp,ftyp_sq_unf,0)
      else
        call file_init(ffinp,fname_inp,ftyp_sq_frm,0)
      end if
      call file_open(ffinp)

      ! array for weights
      allocate(iy_int(4,orb_info%nsym,orb_info%ntoob+orb_info%caborb),
     &     typetab(24),mnmxspc(2,ngas))

      naux_max = max_rank_op('X',oplist%op,.true.)
      hhgg_list = (max_rank_op('AP',oplist%op,.true.)+
     &             max_rank_op('AX',oplist%op,.true.) ) .eq. 0 .or.
     &            (max_rank_op('CP',oplist%op,.true.)+
     &             max_rank_op('CX',oplist%op,.true.) ) .eq. 0

      nh_min = 0
      nh_max = 4
      if (hhgg_list) nh_min = 2

c dbg
c      print *,'import naux_mx',naux_max
c dbg

      do iaux_max = 0, naux_max
        
c          mnmxspc(1,1:ngas) = 0
c          mnmxspc(2,1:ngas) = 4
        select case(iaux_max)
        case(0)
          mnmxspc(1,1:ngas) = 0
          mnmxspc(2,1:ngas) = 4
          do igas = 1, ngas-1
            if (hhgg_list.and.ihpvgas(igas,1).eq.IHOLE) then
              if (ihpvgas(igas+1,1).ne.IHOLE)
     &             mnmxspc(1:2,igas) = (/2,4/) 
            end if
            if (ihpvgas(igas+1,1).eq.IEXTR) then
              mnmxspc(1,igas) = 4
            end if
          end do
          if (ihpvgas(ngas,1).eq.IEXTR) mnmxspc(1,ngas) = 4
        case(1)
          mnmxspc(1,1:ngas) = 0
          mnmxspc(2,1:ngas) = 3
          do igas = 1, ngas-1
            if (hhgg_list.and.ihpvgas(igas,1).eq.IHOLE) then
              if (ihpvgas(igas+1,1).ne.IHOLE)
     &             mnmxspc(1:2,igas) = (/2,3/) 
            end if
            if (ihpvgas(igas+1,1).eq.IEXTR) then
              mnmxspc(1,igas) = 3
            end if
            if (ihpvgas(igas,1).eq.IEXTR) mnmxspc(2,igas) = 4              
          end do
          if (ihpvgas(ngas,1).eq.IEXTR) mnmxspc(1:2,ngas) = 4
        case(2)
          mnmxspc(1,1:ngas) = 0
          mnmxspc(2,1:ngas) = 2
          do igas = 1, ngas-1
            if (hhgg_list.and.ihpvgas(igas,1).eq.IHOLE) then
              if (ihpvgas(igas+1,1).ne.IHOLE)
     &             mnmxspc(1:2,igas) = (/2,2/) 
            end if
            if (ihpvgas(igas+1,1).eq.IEXTR) then
              mnmxspc(1,igas) = 2
            end if
            if (ihpvgas(igas,1).eq.IEXTR) mnmxspc(2,igas) = 4              
          end do
          if (ihpvgas(ngas,1).eq.IEXTR) mnmxspc(1:2,ngas) = 4
        end select
      
        ! set up DA-buffer for 2el integrals
        ntypes = 3
        if (abs(mode).eq.2) ntypes = 6
        if (abs(mode).eq.3) ntypes = 12
        call set_integral_graph(iy_int,nints,4,mnmxspc,orb_info)
        call set_integral_typetab(ntypes,4,typetab)

        nbuffer = nints(1)*ntypes
c dbg
        write(luout,*) 'size of buffer: ',nints(1),ntypes,nbuffer
c dbg
        ifree = mem_setmark('DAimport')
c FUSK FIX - we currently allow up to 3 times of ifree:
        if (nbuffer.gt.3*ifree) then
          write(luout,*) nbuffer,ifree
          call quit(0,'DAimport buffer','definitely too large')
        end if
        !ifree = mem_alloc_real(buffer,nbuffer,'DAbuffer')
        allocate(buffer(nbuffer))        
c
        buffer(1:nbuffer) = 0d0

        call atim_cs(cpux,sysx)

        if (mode.gt.0) then
          call read_mo_file_dalton_spc2(buffer,ffinp,
     &       iaux_max,iaux_max,nh_min,nh_max,
     &       iy_int,typetab,ntypes,orb_info)
        else
          call read_mo_file_dalton_spc(buffer,ffinp,
     &       iaux_max,iaux_max,nh_min,nh_max,
     &       iy_int,typetab,ntypes,orb_info)
        end if

        call atim_cs(cpu,sys)

        if (iprlvl.ge.5) 
     &       call prtim(luout,'time for reading integrals',
     &       cpu-cpux,sys-sysx,-1d0)
     
c dbg
c        idx = 0
c        do idx_int = 1, nints(1)
cc          write(luout,*) 'idx_int = ',idx_int
c          do idx_typ = 1, ntypes
c            idx = idx+1
cc            if (abs(buffer(idx)).gt.10d0)
c                write(luout,*) idx, idx_int, idx_typ, buffer(idx)
c          end do
c        end do
c dbg
        call atim_cs(cpux,sysx)

        ! scaling options
        select case(scaling)
        case(0)
          fac_s = 1d0
          fac_t = 1d0
        case(1)
          fac_s = scal_ab !0.50d0
          fac_t = scal_aa !0.25d0
        case(2)
          fac_s = scal_ab*scal_ab    !0.50d0*0.50d0
          fac_t = scal_aa*scal_aa    !0.25d0*0.25d0
        case default
          call quit(1,'import_2el_dalton','unknown scaling parameter')
        end select

        ! sort into actual list
        call import_list_from_buffer(oplist,buffer,
     &       iaux_max,iaux_max,
     &       fac_s,fac_t,scaling,
     &       iy_int,typetab,ntypes,
     &       anti,
     &       str_info,orb_info)

        call atim_cs(cpu,sys)

        if (iprlvl.ge.5) 
     &       call prtim(luout,'time for sorting integrals',
     &       cpu-cpux,sys-sysx,-1d0)

c FUSK FIX
        deallocate(buffer)
c
        ifree = mem_flushmark('DAimport')
      
      end do

      call file_close_keep(ffinp)

      ifree=mem_flushmark()

      call atim_csw(cpu,sys,wall)

      if (iprlvl.ge.5) 
     &     call prtim(luout,'time in 2int import',
     &     cpu-cpu0,sys-sys0,wall-wall0)


      return

      end

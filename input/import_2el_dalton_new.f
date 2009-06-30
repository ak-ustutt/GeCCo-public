*-----------------------------------------------------------------------
      subroutine import_2el_dalton_new(oplist,fname_inp,
     &     mode,scaling,anti,str_info,orb_info)
*-----------------------------------------------------------------------
*     Routine to read in and reorder 2-electron integrals required for 
*     R12 calculations. Integrals are reordered to be in type order.
*-----------------------------------------------------------------------
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'ioparam.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'hpvxseq.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'par_dalton.h'
      include 'ifc_memman.h'
      include 'def_sort_buffer.h'
      include 'routes.h'

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

      logical ::
     &     fexist, set_ccaa, set_ca, set_abx, set_pq_rs,
     &     incore_sort,
     &     close_again
      integer ::
     &     ifree, lenlist, lenbin, lenbuffer, nbin,
     &     nl1cache,nl2cache, len_rec, n_per_rec,
     &     len_bin, len_bin_min, len_bin_max,
     &     max_mem, maxchain, npass, nbin_per_pass
      real(8) ::
     &     fac_ccaa, fac_aa, fac_ab, fac_abx,
     &     cpu0, sys0, wall0, cpu, sys, wall, cpux, sysx, wallx 

      real(8), pointer ::
     &     buffer(:)
      type(sort_buffer) ::
     &     sbuffer

      type(filinf) ::
     &     ffinp, ffpre, ffchain

      type(filinf), pointer ::
     &     fflist
      type(operator), pointer ::
     &     op


      op =>   oplist%op
      fflist => oplist%fhand

      if(ntest.ge.100)then
        write(luout,*)'=================='
        write(luout,*)'Import-r12-dalton '
        write(luout,*)'=================='
        write(luout,*)'Operator-list: ',trim(oplist%label)
      endif

      ! translate mode and scaling into parameters
      fac_ccaa = 1d0
      fac_aa = 1d0
      fac_ab = 1d0
      fac_abx = 0d0
      set_abx = .false.
      set_ccaa = .true.
      set_ca = .true.
      if (mode.eq.2) then
        set_ccaa = .false.
c        set_ccaa = .true.
        set_ca = .false.
c        fac_ccaa = -1d0
      end if
      if (mode.eq.3) then
        set_ccaa = .false.
        set_ca = .false.
      end if

      if (scaling.eq.1) then
        fac_aa = scal_aa
        fac_ab = scal_ab
        fac_abx = -scal_ab/3d0
        set_abx = .true.
      end if
      if (scaling.eq.2) then
        fac_aa = scal_aa*scal_aa
        fac_ab = (10d0/9d0)*scal_ab*scal_ab
        fac_abx = -(2d0/3d0)*scal_ab*scal_ab
        set_abx = .true.
      end if

      set_pq_rs = .false.
      if (.not.anti) then
        set_pq_rs = .true.
      end if

      call atim_csw(cpu0,sys0,wall0)

      ! set up batching:

      ! get available memory
      ifree=mem_setmark('import_r12')
      ! get length of operator
      lenlist = oplist%len_op

      ! assumed size of L2-cache:
      nl2cache = 2*1024*1024/8
      ! assumed size of L1-cache:
      nl1cache = 16*1024/8

      ! we use 90 per cent of core for batches of final list
      max_mem = ifree*9/10
      if (force_ooc_sort.gt.1) max_mem = lenlist/force_ooc_sort

      incore_sort = max_mem.gt.lenlist

      if (incore_sort) then
        lenbuffer = lenlist
        npass = 1
        len_rec = 32*1024 ! dummy
        if (iprlvl.gt.0) then
          write(luout,'(x,a,g8.3,a)')
     &         'in-core sort active -- allocating ',
     &         dble(lenbuffer)*8d0/(1024d0*1024d0),' Mb'
        end if
      else
        len_rec = 32*1024
        if (force_ooc_sort.gt.1) len_rec = (max_mem/100)*4
        len_bin_min = 10*len_rec
        len_bin_max = lenlist/10
        call import_2el_set_batch(len_bin,n_per_rec,nbin,maxchain,
     &       npass,nbin_per_pass,
     &       lenlist,max_mem,len_rec,
     &       len_bin_min,len_bin_max)
        if (iprlvl.gt.0) then
          write(luout,*) 'out-of-core sort active'
          write(luout,*) ' number of bins:            ',nbin
          write(luout,*) ' size of bin:               ',len_bin
          write(luout,*) ' passes over integral file: ',npass
          write(luout,*) ' number of bins per pass:   ',nbin_per_pass
        end if
      end if

      if (ntest.ge.100) then
        write(luout,*) 'current setting:'
        write(luout,*) 'nl1cache (words) = ',nl1cache
        write(luout,*) 'nl2cache (words) = ',nl2cache
      end if

      if (incore_sort) then
        ! get buffer
        ifree = mem_alloc_real(buffer,lenbuffer,'buffer')
      else
        buffer => null()
        call init_sort_buffer(sbuffer,nbin_per_pass,n_per_rec,maxchain)
      end if

      ! Open the required MO-integral file.
      inquire(file=fname_inp,exist=fexist)
      if(.not.fexist)
     &     call quit(1,'import_2_dalton','No MO integral file: '//
     &       trim(fname_inp))

      call file_init(ffinp,fname_inp,ftyp_sq_unf,0)
      call file_open(ffinp)

      if (.not.incore_sort) then
        call file_init(ffpre,  'LIST_PRESORT',ftyp_da_unf,len_rec)
        call file_init(ffchain,'LIST_CHAIN',  ftyp_sq_unf,0)
        call file_open(ffpre)
        call file_open(ffchain)
      end if

      ! presort integrals from file
      ! if incore, the presorted integrals are on buffer
      call import_2el_presort(buffer,ffpre,ffchain,
     &     incore_sort,sbuffer,npass,len_bin,len_rec,
     &     ffinp,oplist,
     &     fac_aa,fac_ab,fac_abx,fac_ccaa,
     &     set_abx,set_ccaa,set_ca,set_pq_rs,
     &     str_info,orb_info)

      call file_close_keep(ffinp)

      if (.not.incore_sort) then
        call clean_sort_buffer(sbuffer)
      end if

      ! open ME-list file
      close_again = .false.
      if (fflist%unit.le.0) then
        close_again = .true.
        call file_open(fflist)
      end if

      if (incore_sort) then
        call put_vec(fflist,buffer,1,lenlist)
      else
        call atim_csw(cpux,sysx,wallx)
        if (iprlvl.ge.5) 
     &     call prtim(luout,'time in 2int pre-sort',
     &     cpux-cpu0,sysx-sys0,wallx-wall0)
        
        ! out-of-core: process presorted integrals
        call import_2el_sort(fflist,ffpre,ffchain,lenlist,len_bin)

        call atim_csw(cpu,sys,wall)

        if (iprlvl.ge.5) 
     &     call prtim(luout,'time in 2int sort',
     &     cpu-cpux,sys-sysx,wall-wallx)

        call file_close_delete(ffpre)
        call file_close_delete(ffchain)
      end if

      if (close_again) call file_close_keep(fflist)
        
      ifree=mem_flushmark()

      call atim_csw(cpu,sys,wall)

      if (iprlvl.ge.1) 
     &     call prtim(luout,'time in 2int import',
     &     cpu-cpu0,sys-sys0,wall-wall0)


      return

      end

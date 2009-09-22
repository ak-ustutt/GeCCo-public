*----------------------------------------------------------------------*
      subroutine import_cmox_dalton_special(ffcmox,orb_info)
*----------------------------------------------------------------------*
*     read MO-coefficients from SIRIFC (sirius interface file)
*----------------------------------------------------------------------*
      implicit none

c      include 'opdim.h'
      include 'stdunit.h'
      include 'def_orbinf.h'
      include 'def_filinf.h'
      include 'ifc_memman.h'
      include 'multd2h.h'
      include 'par_dalton.h'

      integer, parameter ::
     &     ntest = 100

      type(orbinf), intent(in), target ::
     &     orb_info
      type(filinf), intent(inout) ::
     &     ffcmox

      type(filinf) ::
     &     ffcmox_dalton
      logical ::
     &     closeit, lcmox
      real(8) ::
     &     cpu0,sys0,wall0,cpu,sys,wall
      integer ::
     &     lusir, isym, ifree, igas, igasr,
     &     ncmo, ncmox, imo, imo_reo, ioff, ioff_reo, len,
     &     ioff_sym(8), ioff_mo(8)

      real(8), pointer ::
     &     cmox(:)

      call atim_csw(cpu0,sys0,wall0)

      ifree = mem_setmark('import_cmox')

      ! get buffers (SAO basis -> symmetry blocked)
      ncmo = 0
      ncmox = 0
      imo = 0
      do isym = 1, orb_info%nsym
        ioff_sym(isym) = ncmox
        ioff_mo(isym) =  imo
        ncmo = ncmo +
     &       orb_info%ntoobs(isym)*orb_info%nbas(isym)
        ncmox = ncmox +
     &       orb_info%cab_orb(isym)*
     &       (orb_info%nbas(isym)+orb_info%nxbas(isym))
        imo = imo + orb_info%cab_orb(isym)
      end do
      
      ! buffer for CMOX-matrix as read from AUXBAS or CMOX
      ifree = mem_alloc_real(cmox,ncmox,'cmox')

      ! if CMOX exists, prefer this file
      inquire(file=cmox_file,exist=lcmox)
      if (lcmox) then
        call file_init(ffcmox_dalton,cmox_file,ftyp_sq_unf,0)
        call file_open(ffcmox_dalton)
        call get_cmox(cmox,ffcmox_dalton%unit)
        call file_close_keep(ffcmox_dalton)
      else
        ! too dangerous: doesn't work for canonical aux. bas.!
        call quit(0,'import_cmox_dalton_special',
     &            'CMOX file is missing!')
c        ! else read from AUXBAS
c        call file_init(ffcmox_dalton,auxbas_file,ftyp_sq_frm,0)
c        call file_open(ffcmox_dalton)
c        call get_auxbas(cmox,ffcmox_dalton%unit)
c        call file_close_keep(ffcmox_dalton)
      end if

      if (ntest.ge.100) then
        write(luout,*) 'CMOX:'
        call wr_blkmat(cmox,orb_info%nbas+orb_info%nxbas,
     &                     orb_info%cab_orb,
     &                     orb_info%nsym,0)
      end if

      if (ffcmox%unit.le.0) then
        call file_open(ffcmox)
        closeit = .true.
      else
        closeit = .false.
      end if
      
      ! write to file
      call put_vec(ffcmox,cmox,ncmo+1,ncmo+ncmox)

      if (closeit)
     &     call file_close_keep(ffcmox)

      ifree = mem_flushmark('import_cmox')

      call atim_csw(cpu,sys,wall)

      if (iprlvl.ge.10) 
     &     call prtim(luout,'time in cmox import',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      return
      end

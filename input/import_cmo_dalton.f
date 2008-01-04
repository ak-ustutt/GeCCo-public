*----------------------------------------------------------------------*
      subroutine import_cmo_dalton(ffcmo,orb_info)
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
     &     ntest = 00

      type(orbinf), intent(in), target ::
     &     orb_info
      type(filinf), intent(inout) ::
     &     ffcmo

      type(filinf) ::
     &     ffsir
      logical ::
     &     closeit
      real(8) ::
     &     cpu0,sys0,wall0,cpu,sys,wall
      integer ::
     &     lusir, isym, ifree, igas, igasr,
     &     ncmo, imo, imo_reo, ioff, ioff_reo, len,
     &     ioff_sym(8), ioff_mo(8)

      real(8), pointer ::
     &     cmo(:), cmo_reo(:)

      call atim_csw(cpu0,sys0,wall0)

      ifree = mem_setmark('import_cmo')

      ! get buffers (SAO basis -> symmetry blocked)
      ncmo = 0
      imo = 0
      do isym = 1, orb_info%nsym
        ioff_sym(isym) = ncmo
        ioff_mo(isym) =  imo
        ncmo = ncmo +
     &       orb_info%ntoobs(isym)*orb_info%nbas(isym)
        imo = imo + orb_info%ntoobs(isym)
      end do
      
      ! buffer for CMO-matrix as read from SIRIFC
      ifree = mem_alloc_real(cmo,ncmo,'cmo_ori')
      ! buffer for reordered CMO-matrix
      ifree = mem_alloc_real(cmo_reo,ncmo,'cmo_reo')

      ! open files
      call file_init(ffsir,sirifc,ftyp_sq_unf,0)
      call file_open(ffsir)

      lusir = ffsir%unit
      rewind lusir
      
      if (ffcmo%unit.le.0) then
        call file_open(ffcmo)
        closeit = .true.
      else
        closeit = .false.
      end if

      call mollab('TRCCINT ',lusir,luout)

      ! skip dimension record
      read(lusir)
      ! skip orbital eigenvalue record
      read(lusir)

      ! and finally: the inactive fock matrix in symmetry-blocked
      ! upper triangular form
      read (lusir) cmo(1:ncmo)

      if (ntest.ge.100) then
        write(luout,*) 'CMO (original):'
        call wr_blkmat(cmo,orb_info%nbas,orb_info%ntoobs,
     &                     orb_info%nsym,0)
      end if

      call file_close_keep(ffsir)

      ! and reorder
      ioff_reo = 0
      do imo_reo = 1, orb_info%ntoob
        isym = orb_info%igamorb(imo_reo)
        imo = orb_info%ireots(imo_reo) - ioff_mo(isym)
        len  = orb_info%nbas(isym)
        ioff = ioff_sym(isym) + (imo-1)*len
        cmo_reo(ioff_reo+1:ioff_reo+len) = cmo(ioff+1:ioff+len)
        ioff_reo = ioff_reo+len
      end do

      if (ntest.ge.100) then
        ioff_reo = 1
        do igas = 1, orb_info%ngas
          igasr = orb_info%gas_reo(igas)
          write(luout,*) 'orbital shell # ',
     &         igas,orb_info%igassh(1:orb_info%nsym,igasr)
          call wr_blkmat(cmo_reo(ioff_reo),orb_info%nbas,
     &                                     orb_info%igassh(1,igasr),
     &                                     orb_info%nsym,0)          
          do isym = 1, orb_info%nsym
            ioff_reo = ioff_reo + orb_info%nbas(isym)*
     &                            orb_info%igassh(isym,igasr)
          end do
        end do
      end if
      
      ! write to file
      call put_vec(ffcmo,cmo_reo,1,ncmo)

      if (closeit)
     &     call file_close_keep(ffcmo)

      ifree = mem_flushmark('import_cmo')

      call atim_csw(cpu,sys,wall)

      if (iprlvl.ge.10) 
     &     call prtim(luout,'time in cmo import',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      return
      end

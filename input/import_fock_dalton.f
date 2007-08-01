*----------------------------------------------------------------------*
      subroutine import_fock_dalton(ffham,hop,str_info,orb_info)
*----------------------------------------------------------------------*
*     read fock operator from SIRIFC (sirius interface file)
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'par_globalmarks.h'
      include 'def_operator.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'def_filinf.h'
      include 'ifc_memman.h'
      include 'multd2h.h'
      include 'par_dalton.h'

      integer, parameter ::
     &     ntest = 00

      type(operator), intent(in) ::
     &     hop
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in), target ::
     &     orb_info
      type(filinf), intent(inout) ::
     &     ffham

      type(filinf) ::
     &     ffsir
      logical ::
     &     closeit
      real(8) ::
     &     potnuc,emy,eactiv,emcscf
      ! DALTON writes integer*4, so we must take care of that
      integer*4 ::
     &     istate,ispin,nactel,lsym
      real(8) ::
     &     cpu0,sys0,wall0,cpu,sys,wall
      integer ::
     &     lusir, isym, ifree,
     &     nfock, nh1reo, iocc_cls

      real(8), pointer ::
     &     xfock(:), xh1reo(:)

      real(8), external ::
     &     dnrm2

      call atim_csw(cpu0,sys0,wall0)

      ifree = mem_setmark('import_h1')

      ! get buffers
      nfock = 0
      do isym = 1, orb_info%nsym
        nfock = nfock +
     &       (orb_info%ntoobs(isym)+1)*orb_info%ntoobs(isym)/2
      end do
      
      nh1reo = 0
      do iocc_cls = 1, hop%n_occ_cls
        if (max(hop%ica_occ(1,iocc_cls),hop%ica_occ(2,iocc_cls)).gt.1)
     &       cycle
        nh1reo = nh1reo + hop%len_op_occ(iocc_cls)
      end do

      ! buffer for fock-matrix as read from SIRIFC
      ifree = mem_alloc_real(xfock,nfock,'fock_symblk')
      ! buffer for reordered elements
      if (ffham%buffered) then
        if (ffham%nbuffer.lt.nh1reo)
     &       call quit(1,'import_fock_dalton',
     &       'incore handling for H1 is inconsistent')
        xh1reo => ffham%buffer
        ! check that indices on disc and in memory coincide
        ! (currently assumed in h1_sym2str_reo)
        do iocc_cls = 1, hop%n_occ_cls
          ! ignore R12 stuff
          if(iextr.gt.0.and.
     &       hop%ihpvca_occ(iextr,1,iocc_cls)+
     &       hop%ihpvca_occ(iextr,2,iocc_cls).gt.0 ) cycle

          if (ffham%incore(iocc_cls).gt.0.and.
     &        ffham%idxrec(iocc_cls).ne.hop%off_op_occ(iocc_cls)) then
            call quit(1,'import_fock_dalton',
     &           'non-contiguous buffering: not prepared for that')
          end if
        end do
      else
        ifree = mem_alloc_real(xh1reo,nh1reo,'fock_reo')
      end if

      ! open files
      call file_init(ffsir,sirifc,ftyp_sq_unf,0)
      call file_open(ffsir)

      lusir = ffsir%unit
      rewind lusir
      
      if (ffham%unit.le.0) then
        call file_open(ffham)
        closeit = .true.
      else
        closeit = .false.
      end if

      call mollab('SIR IPH ',lusir,luout)

      read (lusir) potnuc,emy,eactiv,emcscf,istate,ispin,nactel,lsym
      read (lusir) ! overread dimension information
      read (lusir) ! overread CMO data
      read (lusir) ! overread DV data
      read (lusir) ! overread F data
      read (lusir) ! overread PV data
      ! and finally: the inactive fock matrix in symmetry-blocked
      ! upper triangular form
      read (lusir,err=16) xfock(1:nfock)
      goto 1
      ! on error try next record
 16   write(luout,*) 'Trying new DALTON format ...'
      read (lusir) xfock(1:nfock) 

  1   call file_close_keep(ffsir)

      if (dnrm2(nfock,xfock,1).lt.1d-12)
     &   call quit(0,'import_fock_dalton',
     &               'No sensible fock matrix found!')

      ! and reorder
      call h1_sym2str_reo(emcscf,xfock,xh1reo,hop,str_info,orb_info)

      call put_vec(ffham,xh1reo,1,nh1reo)

c      call
      if (closeit)
     &     call file_close_keep(ffham)

      ifree = mem_flushmark('import_h1')

      call atim_csw(cpu,sys,wall)

      if (iprlvl.ge.5) 
     &     call prtim(luout,'time in 1int import',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      return
      end

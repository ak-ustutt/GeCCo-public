*----------------------------------------------------------------------*
      subroutine import_exchange_dalton(klist,kname,str_info,orb_info)
*----------------------------------------------------------------------*
*     Read the exchange matrix from the file MO_K.
*     This is required for approximation B and beyond in R12 models.
*     Modified by GWR from import_fock_dalton.f (Feb 2008)
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
      include 'def_me_list.h'
      include 'ifc_memman.h'
      include 'multd2h.h'
      include 'par_dalton.h'

      integer, parameter ::
     &     ntest = 00

      type(me_list), intent(in) ::
     &     klist
      character*256, intent(in) ::
     &     kname
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in), target ::
     &     orb_info

      type(filinf) ::
     &     ffmok
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
     &     lumok, isym, ifree, irec,
     &     nexch, nh1reo, iocc_cls

      type(operator), pointer ::
     &     kop
      type(filinf), pointer ::
     &     ffexc

      real(8), pointer ::
     &     kexch(:), xh1reo(:)

      real(8), external ::
     &     dnrm2

      call atim_csw(cpu0,sys0,wall0)

      ffexc => klist%fhand
      kop   => klist%op

      ifree = mem_setmark('import_k1')

      ! get buffers
      nexch = 0
      do isym = 1, orb_info%nsym
        nexch = nexch +
     &       (orb_info%ntoobs(isym)+1)*orb_info%ntoobs(isym)/2
      end do
      
      nh1reo = 0
      do iocc_cls = 1, kop%n_occ_cls
        if (max(kop%ica_occ(1,iocc_cls),kop%ica_occ(2,iocc_cls)).gt.1)
     &       cycle
        if(kop%formal_blk(iocc_cls))cycle
        nh1reo = nh1reo + klist%len_op_occ(iocc_cls)
      end do

      ! buffer for exchange-matrix as read from MO_K
      ifree = mem_alloc_real(kexch,nexch,'exchange_symblk')
      ! buffer for reordered elements
      if (ffexc%buffered) then
        if (ffexc%nbuffer.lt.nh1reo)
     &       call quit(1,'import_exchange_dalton',
     &       'incore handling for K1 is inconsistent')
        xh1reo => ffexc%buffer
        ! check that indices on disc and in memory coincide
        ! (currently assumed in h1_sym2str_reo)
        do iocc_cls = 1, kop%n_occ_cls
          if (ffexc%incore(iocc_cls).gt.0.and.
     &        ffexc%idxrec(iocc_cls).ne.klist%off_op_occ(iocc_cls)) then
            call quit(1,'import_exchange_dalton',
     &           'non-contiguous buffering: not prepared for that')
          end if
        end do
      else
        ifree = mem_alloc_real(xh1reo,nh1reo,'exch_reo')
      end if

      ! open files
      call file_init(ffmok,trim(kname),ftyp_sq_frm,0)
      call file_open(ffmok)

      lumok = ffmok%unit
      rewind lumok
      
      if (ffexc%unit.le.0) then
        call file_open(ffexc)
        closeit = .true.
      else
        closeit = .false.
      end if

      ! and finally: the inactive exchange matrix in symmetry-blocked
      ! upper triangular form
      read (lumok,err=16) kexch(1:nexch)
      goto 1
      ! on error try next record
 16   write(lulog,*) 'Trying new DALTON format ...'
      read (lumok) kexch(1:nexch) 

  1   call file_close_keep(ffmok)

      if (dnrm2(nexch,kexch,1).lt.1d-12)
     &   call quit(0,'import_exchange_dalton',
     &               'No sensible fock matrix found!')

      ! and reorder
      call h1_sym2str_reo(emcscf,kexch,xh1reo,klist,str_info,orb_info)

      call put_vec(ffexc,xh1reo,1,nh1reo)

      if (closeit)
     &     call file_close_keep(ffexc)

      ifree = mem_flushmark('import_k1')

      call atim_csw(cpu,sys,wall)

      if (iprlvl.ge.5) 
     &     call prtim(lulog,'time in exchange integral import',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      return
      end

*----------------------------------------------------------------------*
      subroutine import_fock_dalton(hlist,str_info,orb_info,
     &     use_file,name)
*----------------------------------------------------------------------*
*     read fock operator from SIRIFC (sirius interface file)
*     or from file "name", if use_file is set to .true.
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
      include 'par_gamess.h'

      integer, parameter ::
     &     ntest = 00

      type(me_list), intent(in) ::
     &     hlist
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in), target ::
     &     orb_info
      logical, intent(in) ::
     &     use_file
      character*(*) ::
     &     name

      logical ::
     &     closeit
      real(8) ::
     &     eref
      real(8) ::
     &     cpu0,sys0,wall0,cpu,sys,wall
      integer ::
     &     lusir, isym, ifree, irec,
     &     nfock, nh1reo, iocc_cls,
     &     idxst, idxnd, idxbufst, idxbufnd

      type(operator), pointer ::
     &     hop
      type(filinf), pointer ::
     &     ffham

      real(8), pointer ::
     &     xfock(:), xh1reo(:)

      if (ntest.ge.100)
     &     call write_title(luout,wst_dbg_subr,'import_fock_dalton')

      call atim_csw(cpu0,sys0,wall0)

      ffham => hlist%fhand
      hop   => hlist%op

      ifree = mem_setmark('import_h1')

      ! get buffers
      if (.not.use_file) then
        nfock = 0
        do isym = 1, orb_info%nsym
          nfock = nfock +
     &       (orb_info%ntoobs(isym)+1)*orb_info%ntoobs(isym)/2
        end do
      else
        nfock = (orb_info%ntoob+orb_info%caborb)*
     &          (orb_info%ntoob+orb_info%caborb+1)/2
      end if
      
      nh1reo = 0
      do iocc_cls = 1, hop%n_occ_cls
        if (max(hop%ica_occ(1,iocc_cls),hop%ica_occ(2,iocc_cls)).gt.1)
     &       cycle
        if(hop%formal_blk(iocc_cls))cycle
        ! Quick fix to ignore Fock operators with an external index.
        if(.not.use_file.and.
     &      (hop%ihpvca_occ(iextr,1,iocc_cls).gt.0.or.
     &       hop%ihpvca_occ(iextr,2,iocc_cls).gt.0))cycle
        nh1reo = nh1reo + hlist%len_op_occ(iocc_cls)
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
          if (hop%formal_blk(iocc_cls)) cycle
          ! ignore R12 stuff unless use_file is selected
          if(.not.use_file.and.
     &       iextr.gt.0.and.
     &       hop%ihpvca_occ(iextr,1,iocc_cls)+
     &       hop%ihpvca_occ(iextr,2,iocc_cls).gt.0 ) cycle

          if (ffham%incore(iocc_cls).gt.0.and.
     &        ffham%idxrec(iocc_cls).ne.hlist%off_op_occ(iocc_cls)) then
            call quit(1,'import_fock_dalton',
     &           'non-contiguous buffering: not prepared for that')
          end if
        end do
      else
        ifree = mem_alloc_real(xh1reo,nh1reo,'fock_reo')
      end if

      ! read fock matrix and reorder
      if (.not.use_file) then
        call read_fock_from_sirifc(eref,xfock,nfock)
        call h1_sym2str_reo(eref,xfock,xh1reo,hlist,str_info,orb_info)
      else
        select case(trim(name))
        case (moints) ! this option is for GAMESS input
          call read_fock_from_file_gms(eref,xfock,nfock)
        case default
          call read_fock_from_file(eref,xfock,nfock,name)
        end select
        call h1_full2str_reo(eref,xfock,xh1reo,hlist,str_info,orb_info)
      end if

      if (ffham%unit.le.0) then
        call file_open(ffham)
        closeit = .true.
      else
        closeit = .false.
      end if

      ! well, on file the h1-blocks may be non-contiguous, so
      ! we have to collect the adjacent blocks and write, as 
      ! soon as a different buffer comes in between
      idxbufst = 1 ! counter for buffer
      idxbufnd = 0
      idxst = 1    ! counter for list on disc
      idxnd = 0
      do iocc_cls = 1, hop%n_occ_cls+1
c dbg fix by mh
        if (iocc_cls.le.hop%n_occ_cls) then
c dbg original
        if (iocc_cls.le.hop%n_occ_cls.and.hop%formal_blk(iocc_cls)) 
     &     cycle
c dbg resume fix
        end if
c dbg end fix
        ! on any of these conditions, we have to write the present
        ! buffer slice to disc
c dbg fix by mh
        if (iocc_cls.gt.hop%n_occ_cls) then
          if (idxbufnd.ge.idxbufst)
     &         call put_vec(ffham,xh1reo(idxbufst),idxst,idxnd)
          idxbufst = idxbufnd+1
          idxst = idxnd+1
        else if (iocc_cls.gt.hop%n_occ_cls .or.
c dbg original line        if (iocc_cls.gt.hop%n_occ_cls .or.
c dbg end fix
     &      max(hop%ica_occ(1,iocc_cls),hop%ica_occ(2,iocc_cls)).gt.1
     &       .or. (.not.use_file.and.
     &      (hop%ihpvca_occ(iextr,1,iocc_cls).gt.0.or.
     &       hop%ihpvca_occ(iextr,2,iocc_cls).gt.0)) ) then
          ! write, if applicable
          if (idxbufnd.ge.idxbufst)
     &         call put_vec(ffham,xh1reo(idxbufst),idxst,idxnd)
          idxbufst = idxbufnd+1
c dbg fix by mh
          if (iocc_cls.le.hop%n_occ_cls) then
c dbg original
          idxnd = idxnd+hlist%len_op_occ(iocc_cls)
c dbg resume fix
          end if
c dbg end fix
          idxst = idxnd+1
        else
c dbg fix by mh
          if (iocc_cls.le.hop%n_occ_cls) then
c dbg original
          idxnd    = idxnd    + hlist%len_op_occ(iocc_cls)
          idxbufnd = idxbufnd + hlist%len_op_occ(iocc_cls)
c dbg resume fix
          end if
c dbg end fix
        end if
      end do

      if (closeit)
     &     call file_close_keep(ffham)

      ifree = mem_flushmark('import_h1')

      call atim_csw(cpu,sys,wall)

      if (iprlvl.ge.5) 
     &     call prtim(luout,'time in 1int import',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      return
      end

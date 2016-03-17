*----------------------------------------------------------------------*
      subroutine import_fock_cfour(hlist,str_info,orb_info)
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
      include 'def_me_list.h'
      include 'ifc_memman.h'
      include 'multd2h.h'

      character(len=17), parameter ::
     &     i_am = 'import_fock_cfour'
 
      integer, parameter ::
     &     ntest = 00

      type(me_list), intent(in) ::
     &     hlist
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in), target ::
     &     orb_info

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
     &     call write_title(lulog,wst_dbg_subr,i_am)

      call atim_csw(cpu0,sys0,wall0)

      ffham => hlist%fhand
      hop   => hlist%op

      ifree = mem_setmark('import_h1')

      nh1reo = 0
      do iocc_cls = 1, hop%n_occ_cls
        if (max(hop%ica_occ(1,iocc_cls),hop%ica_occ(2,iocc_cls)).gt.1)
     &       cycle
        if(hop%formal_blk(iocc_cls))cycle
        ! Quick fix to ignore Fock operators with an external index.
        if(
     &      (hop%ihpvca_occ(iextr,1,iocc_cls).gt.0.or.
     &       hop%ihpvca_occ(iextr,2,iocc_cls).gt.0))cycle
        nh1reo = nh1reo + hlist%len_op_occ(iocc_cls)
      end do

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
          if(
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

      nfock = (orb_info%ntoob+orb_info%caborb)*
     &        (orb_info%ntoob+orb_info%caborb+1)/2

      ifree = mem_alloc_real(xfock,nfock,'fock_symblk')

      ! read fock matrix and reorder
      call read_fock_from_file_cfour(eref,xfock,nfock,'fockin.dat')
      call h1_full2str_reo(eref,xfock,xh1reo,hlist,str_info,orb_info)

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
     &       .or. (
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
     &     call prtim(lulog,'time in 1int import',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      return
      end

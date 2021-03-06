*----------------------------------------------------------------------*
      subroutine join_operator(op1,op2,orb_info)
*----------------------------------------------------------------------*
*     add occupations from op2 to op1 (if not already present in op1)
*     to be called in stage 0 of operator initialization only
*     (i.e. before setting dimensions etc.)
*----------------------------------------------------------------------*

      implicit none
      include 'opdim.h'
      include 'stdunit.h'
      include 'def_operator.h'
      include 'def_orbinf.h'
      include 'ifc_memman.h'
      include 'par_globalmarks.h'

      integer, parameter ::
     &     ntest = 00

      type(operator), intent(inout) ::
     &     op1
      type(operator), intent(in) ::
     &     op2
      type(orbinf), intent(in) ::
     &     orb_info

      integer ::
     &     njoined, nblk, nblk_old, ioffblk2, ioffblk,
     &     iblk, iblk2, idxblk, nadd, ngas, nspin, ijoin, ifree
      integer, pointer ::
     &     ihpvca_occ_new(:,:,:),
     &     ica_occ_new(:,:),
     &     igasca_restr_new(:,:,:,:,:,:),
     &     idx_graph_new(:,:,:),
     &     blk_version_new(:)
      logical, pointer ::
     &     formal_blk_new(:)
      
      integer, external ::
     &     iblk_occ

      ngas = orb_info%ngas
      nspin = orb_info%nspin

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'join_operator')
        write(lulog,*) 'joining ',trim(op1%name),' and ',trim(op2%name)
      end if

      if (op1%type.ne.op2%type)
     &     call quit(1,'join_operator',
     &     'cannot join operators with different type')

      if (op1%njoined.ne.op2%njoined)
     &     call quit(1,'join_operator',
     &     'cannot join operators with different number'//
     &     ' of primitive vertices')
      
      njoined = op1%njoined

      nadd = 0
      do iblk2 = 1, op2%n_occ_cls
        idxblk = (iblk2-1)*njoined+1
        if (iblk_occ(op2%ihpvca_occ(1,1,idxblk),.false.,op1,
     &      op2%blk_version(iblk2)).le.0)
     &       nadd = nadd+1
      end do

      if (nadd.gt.0) then
        nblk_old = op1%n_occ_cls
        nblk = nblk_old + nadd
        allocate(ihpvca_occ_new(ngastp,2,nblk*njoined),
     &         igasca_restr_new(2,orb_info%ngas,2,2,nspin,nblk*njoined),
     &           ica_occ_new(2,nblk),
     &           formal_blk_new(nblk),blk_version_new(nblk))
        ! save old info
        ihpvca_occ_new(1:ngastp,1:2,1:nblk_old*njoined) =
     &       op1%ihpvca_occ(1:ngastp,1:2,1:nblk_old*njoined)
        igasca_restr_new(1:2,1:ngas,1:2,1:2,1:nspin,1:nblk_old*njoined)=
     &       op1%igasca_restr(1:2,1:ngas,1:2,1:2,1:nspin,
     &                                              1:nblk_old*njoined)
        ica_occ_new(1:2,1:nblk_old) =
     &       op1%ica_occ(1:2,1:nblk_old)
        formal_blk_new(1:nblk_old) =
     &       op1%formal_blk(1:nblk_old)
        blk_version_new(1:nblk_old) =
     &       op1%blk_version(1:nblk_old)

        ! add new blocks
        iblk = nblk_old
        do iblk2 = 1, op2%n_occ_cls
          ioffblk2 = (iblk2-1)*njoined
          if (iblk_occ(op2%ihpvca_occ(1,1,ioffblk2+1),
     &                 .false.,op1,op2%blk_version(iblk2)).le.0) then
            iblk = iblk+1
            if (iblk.gt.nblk)
     &           call quit(1,'join_operator','something''s wrong...')
            ioffblk = (iblk-1)*njoined
            ihpvca_occ_new(1:ngastp,1:2,ioffblk+1:ioffblk+njoined) =
     &          op2%ihpvca_occ(1:ngastp,1:2,ioffblk2+1:ioffblk2+njoined)
            igasca_restr_new(1:2,1:ngas,1:2,1:2,1:nspin,
     &                                  ioffblk+1:ioffblk+njoined) =
     &           op2%igasca_restr(1:2,1:ngas,1:2,1:2,1:nspin,
     &                                  ioffblk2+1:ioffblk2+njoined)
            ica_occ_new(1:2,iblk) =
     &           op2%ica_occ(1:2,iblk2)
            formal_blk_new(iblk) = op2%formal_blk(iblk2)
            blk_version_new(iblk) = op2%blk_version(iblk2)
          end if
        end do

        op1%n_occ_cls = nblk
        
        deallocate(op1%ihpvca_occ,op1%igasca_restr,
     &       op1%ica_occ,op1%formal_blk,op1%blk_version)

        op1%ihpvca_occ => ihpvca_occ_new
        op1%igasca_restr => igasca_restr_new
        op1%ica_occ => ica_occ_new
        op1%formal_blk => formal_blk_new
        op1%blk_version => blk_version_new

        call mem_pushmark()
        ifree = mem_gotomark(operator_def)

        ifree = mem_dealloc(trim(op1%name))
        ifree = mem_register(2*ngastp*nblk*njoined
     &                      +2*nblk
     &                      +8*orb_info%ngas*nblk*njoined
     &                      +nblk,
     &       trim(op1%name))

        call mem_popmark()

      end if

      if (ntest.ge.100) then
        write(lulog,*) 'generated ',op1%n_occ_cls,' blocks'
        do iblk = 1, op1%n_occ_cls
          write(lulog,'(/x,a,i4)') 'Occupation Nr. ',iblk
          ioffblk = (iblk-1)*njoined
          call wrt_occ_n(lulog,op1%ihpvca_occ(1,1,ioffblk+1),njoined)
          do ijoin = 1, njoined
            call wrt_rstr(lulog,
     &           op1%igasca_restr(1,1,1,1,1,ioffblk+ijoin),ngas)
          end do
        end do
      end if
      return
      end

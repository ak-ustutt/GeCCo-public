*----------------------------------------------------------------------*
      subroutine set_opti_info_signs(opti_info,mode,nopt,
     &               me_trv,me_mvp,me_met,me_rhs,use_s)
*----------------------------------------------------------------------*
*     evaluates the sign for the contraction vector*(matrix*vector),
*     which is later (optc_update_redsp3, optc_update_redsp4,
*     optc_orthvec) done by a simple scalar product.
*     Relevant, when the operator assigned to (matrix*vector) has two
*     vertices --> sign can be different for each occupation class.
*
*     matthias, fall 2009 (adopted from set_opti_info)
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     maxsec = 30
      
      include 'opdim.h'
      include 'stdunit.h'
      include 'def_optimize_info.h'
      include 'ifc_input.h'
      include 'ifc_memman.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'def_me_list_array.h'
      include 'hpvxseq.h'
      include 'ifc_operators.h'

      type(optimize_info), intent(inout) ::
     &     opti_info
      integer, intent(in) ::
     &     nopt, mode
      type(me_list_array), intent(in) ::
     &     me_trv(*), me_mvp(*), me_met(*), me_rhs(*)
      logical, intent(in) ::
     &     use_s(*)

      integer ::
     &     ifree, iopt, ioff, nsect, nwfpsec(maxsec), idstsec(maxsec),
     &     sign_cur, sign_old, iblk, ncblk, nablk, ncblk2, nablk2,
     &     iocc_dag(ngastp,2), iocc0(ngastp,2),iocc_cntr(ngastp,2),
     &     i0(ngastp,2),j0(ngastp,2)
      integer(8) ::
     &     base, topo_trv, topo1, topo2
      integer, pointer ::
     &     nsec, iocc_trv(:,:), iocc_mvp(:,:), iocc1(:,:), iocc2(:,:),
     &     idx_graph(:,:,:), hpvx_csec(:), hpvx_asec(:), graph_csec(:),
     &     graph_asec(:), graph_csec2(:), graph_asec2(:),
     &     iocc3(:,:), iocc4(:,:)
      real(8) ::
     &     signsec(maxsec)
      type(operator), pointer ::
     &     op_trv, op_mvp, op_met, op_rhs

      integer, external ::
     &     sign_hpvx, sign_contr
!     integer, external :: 
!    &     iocc_overlap(ngastp,2)
!     logical, external ::
!    &     iocc_equal
      integer(8), external ::
     &     int8_pack

c      if (mode.ne.2.and.mode.ne.3) call quit(1,'set_opti_info_signs',
c     &                             'only for LEQ and EVP solvers yet')

      base = pack_base
      ifree = mem_alloc_int (opti_info%nsec,nopt,'nsec')

      nsect = 0
      do iopt = 1, nopt
        op_trv => me_trv(iopt)%mel%op
        op_mvp => me_mvp(iopt)%mel%op
        if (use_s(iopt)) op_met => me_met(iopt)%mel%op
        if (.not.use_s(iopt)) op_met => me_mvp(iopt)%mel%op  ! do not leave undefined
        if (mode.eq.2) op_rhs => me_rhs(iopt)%mel%op
        if (mode.ne.2) op_rhs => me_mvp(iopt)%mel%op
        nsec => opti_info%nsec(iopt)
        nsec = 0
        sign_old = 0

        if (op_trv%n_occ_cls.ne.op_mvp%n_occ_cls.or.
     &      use_s(iopt).and.op_trv%n_occ_cls.ne.op_met%n_occ_cls.or.
     &      mode.eq.2.and.op_trv%n_occ_cls.ne.op_rhs%n_occ_cls)
     &      call quit(1,'set_opti_info_signs','inconsistent n_occ_cls!')
        if ((use_s(iopt).and.op_mvp%njoined.ne.op_met%njoined).or.
     &      (mode.eq.2.and.op_mvp%njoined.ne.op_rhs%njoined))
     &      call quit(1,'set_opti_info_signs','inconsistent njoined!')

c dbg
c        print *,'njoined for trv, mvp: ',op_trv%njoined, op_mvp%njoined
c dbgend
        if (op_trv%njoined.eq.1.and.op_mvp%njoined.eq.1) then
          do iblk = 1, op_trv%n_occ_cls
            iocc_trv => op_trv%ihpvca_occ(1:ngastp,1:2,iblk)
            iocc1 => op_mvp%ihpvca_occ(1:ngastp,1:2,iblk)
            if (.not.iocc_equal(iocc_trv,.false.,iocc1,.false.).and.
     &          .not.iocc_equal(iocc_trv,.false.,iocc1,.true.))
     &            call quit(1,'set_opti_info_signs',
     &                      'operators do not match (1a)!')
            if (use_s(iopt)) then
              iocc1 => op_met%ihpvca_occ(1:ngastp,1:2,iblk)
              if (.not.iocc_equal(iocc_trv,.false.,iocc1,.false.).and.
     &            .not.iocc_equal(iocc_trv,.false.,iocc1,.true.))
     &             call quit(1,'set_opti_info_signs',
     &                       'operators do not match (1b)!')
            end if
            if (mode.eq.2) then
              iocc1 => op_rhs%ihpvca_occ(1:ngastp,1:2,iblk)
              if (.not.iocc_equal(iocc_trv,.false.,iocc1,.false.).and.
     &            .not.iocc_equal(iocc_trv,.false.,iocc1,.true.))
     &             call quit(1,'set_opti_info_signs',
     &                       'operators do not match (1c)!')
            end if
          end do
          nsec = 1
          nsect = nsect + 1
          if (nsect.gt.maxsec) call quit(1,'set_opti_info_signs',
     &                'increase maxsec!')
          nwfpsec(nsect) = me_trv(iopt)%mel%len_op
          idstsec(nsect) = 1
          signsec(nsect) = 1d0

        else if (op_trv%njoined.eq.1.and.op_mvp%njoined.eq.2) then
          do iblk = 1, op_trv%n_occ_cls
c dbg
c            print *,'current occupation class: ',iblk
c dbgend
            iocc_trv => op_trv%ihpvca_occ(1:ngastp,1:2,iblk)
            iocc1 => op_mvp%ihpvca_occ(1:ngastp,1:2,iblk*2-1)
            iocc2 => op_mvp%ihpvca_occ(1:ngastp,1:2,iblk*2)
c dbg
c        print *, 'test'
c        call wrt_occ(6,iocc1)
c        call wrt_occ(6,iocc2)
c        call wrt_occ(6,iocc_trv)
c dbgend
            iocc_dag(1:ngastp,1) = iocc_trv(1:ngastp,2)
            iocc_dag(1:ngastp,2) = iocc_trv(1:ngastp,1)
            topo_trv = int8_pack(iocc_trv,ngastp*2,base)
            topo1 = int8_pack(iocc1,ngastp*2,base)
            topo2 = int8_pack(iocc2,ngastp*2,base)
            if (topo_trv.ne.topo1+topo2) call quit(1,
     &            'set_opti_info_signs','operators do not match (2a)!')
            if (use_s(iopt)) then
              iocc3 => op_met%ihpvca_occ(1:ngastp,1:2,iblk*2-1)
              iocc4 => op_met%ihpvca_occ(1:ngastp,1:2,iblk*2)
              if (.not.iocc_equal(iocc1,.false.,iocc3,.false.).or.
     &            .not.iocc_equal(iocc2,.false.,iocc4,.false.))
     &             call quit(1,'set_opti_info_signs',
     &                       'operators do not match (2b)!')
            end if
            if (mode.eq.2) then
              iocc3 => op_rhs%ihpvca_occ(1:ngastp,1:2,iblk*2-1)
              iocc4 => op_rhs%ihpvca_occ(1:ngastp,1:2,iblk*2)
              if (.not.iocc_equal(iocc1,.false.,iocc3,.false.).or.
     &            .not.iocc_equal(iocc2,.false.,iocc4,.false.))
     &             call quit(1,'set_opti_info_signs',
     &                       'operators do not match (2c)!')
            end if

            ! check if graph indices match --> lists in same order
            call get_num_subblk(ncblk,nablk,
     &             op_trv%ihpvca_occ(1,1,iblk),1)
            call get_num_subblk(ncblk2,nablk2,
     &             op_mvp%ihpvca_occ(1,1,iblk*2-1),2)
            if (ncblk.ne.ncblk2.or.nablk.ne.nablk2)
     &            call quit(1,'set_opti_info_signs',
     &                      'operators do not match (3)!')
            allocate(hpvx_csec(ncblk),hpvx_asec(nablk),
     &               graph_csec(ncblk), graph_asec(nablk),
     &               graph_csec2(ncblk), graph_asec2(nablk))
            idx_graph => me_trv(iopt)%mel%idx_graph(1:ngastp,1:2,
     &                                              iblk:iblk)
            call condense_occ(graph_csec,graph_asec,
     &                  hpvx_csec,hpvx_asec,
     &                  idx_graph,1,hpvxblkseq)
            idx_graph => me_mvp(iopt)%mel%idx_graph(1:ngastp,1:2,
     &                                              iblk*2-1:iblk*2)
            call condense_occ(graph_csec2,graph_asec2,
     &                  hpvx_csec,hpvx_asec,
     &                  idx_graph,2,hpvxblkseq)
            if (.not.all(graph_csec(1:ncblk)-graph_csec2(1:ncblk).eq.0)
     &          .or.
     &          .not.all(graph_asec(1:nablk)-graph_asec2(1:nablk).eq.0))
     &            call quit(1,'set_opti_info_signs',
     &                      'graph indices do not match!')
            deallocate(hpvx_csec,hpvx_asec,
     &                 graph_csec, graph_asec,
     &                 graph_csec2, graph_asec2)

            ! contraction sign: consider contraction occ1 - occ_trv
            ! the contraction is given by occ1
!           iocc_dag(1:ngastp,1)=iocc_dag(1:ngastp,1)-iocc1(1:ngastp,2)
!           iocc_dag(1:ngastp,2)=iocc_dag(1:ngastp,2)-iocc1(1:ngastp,1)
!           iocc0 = 0
!           ! get occ_trv in contraction form: KA^+ J0C J0A KC^+
!           sign_cur = sign_hpvx(1,iocc1,.true.,iocc_dag,.false.)
!           ! now the contraction:
!           sign_cur = sign_cur 
!    &               * sign_contr(iocc1,iocc0,iocc_dag,0,.false.)

!!try to get the sign in a generalised way
!!          contraction between first two vertices:
!           print*, 'debugging set_opti_info_sign for', iopt, iblk
!           call wrt_occ(6,iocc1)
!           call wrt_occ(6,iocc_dag)
            iocc_cntr = iocc_overlap(iocc1,.false.,iocc_dag,.true.)
!           call wrt_occ(6,iocc_cntr)

            i0(1:ngastp,1)=iocc1(1:ngastp,1)-iocc_cntr(1:ngastp,1)
            i0(1:ngastp,2)=iocc1(1:ngastp,2)-iocc_cntr(1:ngastp,2)
!           call wrt_occ(6,i0)
             
            j0(1:ngastp,1)=iocc_dag(1:ngastp,1)-iocc_cntr(1:ngastp,2)
            j0(1:ngastp,2)=iocc_dag(1:ngastp,2)-iocc_cntr(1:ngastp,1)
!           call wrt_occ(6,j0)
      
            sign_cur = sign_hpvx(1,iocc_cntr,.false.,i0,.false.)
!           print*, 'sign after first sign_hpvx:',sign_cur
            sign_cur = sign_cur*sign_hpvx(1,iocc_cntr,.true.,j0,.false.)
!           print*, 'sign after second sign_hpvx:',sign_cur
            sign_cur = sign_cur*
     &           sign_contr(iocc_cntr,i0,j0,0,.false.)
!           print*, 'sign after first sign_contr:',sign_cur
           
!!          contraction between second and third vertices:
            iocc_dag = j0
!           call wrt_occ(6,iocc_dag)
!           call wrt_occ(6,iocc2)
            iocc_cntr = iocc_overlap(iocc_dag,.false.,iocc2,.true.)
!           call wrt_occ(6,iocc_cntr)

            i0(1:ngastp,1)=iocc_dag(1:ngastp,1)-iocc_cntr(1:ngastp,1)
            i0(1:ngastp,2)=iocc_dag(1:ngastp,2)-iocc_cntr(1:ngastp,2)
!           call wrt_occ(6,i0)
             
            j0(1:ngastp,1)=iocc2(1:ngastp,1)-iocc_cntr(1:ngastp,2)
            j0(1:ngastp,2)=iocc2(1:ngastp,2)-iocc_cntr(1:ngastp,1)
!           call wrt_occ(6,j0)
      
            sign_cur = sign_cur*sign_hpvx(1,iocc_cntr,.false.,
     &                                    i0,.false.)
!           print*, 'sign after third sign_hpvx:',sign_cur
            sign_cur = sign_cur*sign_hpvx(1,iocc_cntr,.true.,j0,.false.)
!           print*, 'sign after fourth sign_hpvx:',sign_cur
            sign_cur = sign_cur*
     &           sign_contr(iocc_cntr,i0,j0,0,.false.)
!           print*, 'sign after second sign_contr:',sign_cur

!!done
            if (sign_cur.eq.sign_old) then
              ! same section: just add length of current block
              nwfpsec(nsect) = nwfpsec(nsect)
     &                       + me_trv(iopt)%mel%len_op_occ(iblk)
            else
              ! create new section
              nsec = nsec + 1
              nsect = nsect + 1
              if (nsect.gt.maxsec) call quit(1,'set_opti_info_signs',
     &                    'increase maxsec!')
              nwfpsec(nsect) = me_trv(iopt)%mel%len_op_occ(iblk)
              idstsec(nsect) = me_trv(iopt)%mel%off_op_occ(iblk) + 1
              signsec(nsect) = dble(sign_cur)
c              signsec(nsect) = 1d0
              sign_old = sign_cur
            end if
          end do
        else
          call quit(1,'set_opti_info_signs',
     &              'not yet adapted for this case')
        end if
      end do

      ifree = mem_alloc_int(opti_info%nwfpsec,nsect,'nwfpsec')
      ifree = mem_alloc_int(opti_info%idstsec,nsect,'idstsec')
      ifree = mem_alloc_real(opti_info%signsec,nsect,'signsec')

      opti_info%nwfpsec(1:nsect) = nwfpsec(1:nsect)
      opti_info%idstsec(1:nsect) = idstsec(1:nsect)
      opti_info%signsec(1:nsect) = signsec(1:nsect)

c dbg
c      write(lulog,*) 'set_opti_info_signs:'
c      write(lulog,*) '--------------------'
c      write(lulog,'(a,10i8)') 'nsec: ', opti_info%nsec(1:nopt)
c      write(lulog,'(a,30i8)') 'nwfpsec: ',opti_info%nwfpsec(1:nsect)
c      write(lulog,'(a,30i8)') 'idstsec: ',opti_info%idstsec(1:nsect)
c      write(lulog,'(a,30f8.1)') 'signsec : ',opti_info%signsec(1:nsect)
c dbgend

      return
      end

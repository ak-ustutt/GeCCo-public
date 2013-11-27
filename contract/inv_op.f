      subroutine inv_op(ninp,label_inp,nlist,label_inv,mode,
     &     op_info,orb_info,str_info,strmap_info)
*----------------------------------------------------------------------*
*     Wrapper subroutine used in the inversion of the matrix 
*     representation of an operator, op_inp.
*     The resultant is op_inv.
*     Optionally, the unitary matrix is returned on a second output list
*     GWR November 2007.
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 00
      
      include 'stdunit.h'
      include 'opdim.h'
      include 'def_orbinf.h'
      include 'mdef_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'ifc_memman.h'
      include 'par_opnames_gen.h'

      integer, intent(in) ::
     &     ninp, nlist
      character(*), intent(in) ::
     &     label_inp(ninp), label_inv(nlist)
      type(operator_info), intent(inout) ::
     &     op_info
      type(orbinf), intent(in) ::
     &     orb_info
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf), intent(in) ::
     &     strmap_info
      character(len=*), intent(in) ::
     &     mode

      type(me_list), pointer ::
     &     me_inp, me_inv, me_u, me_spc

      integer ::
     &     idx_inp, idx_inv, idx_u, idx_spc,
     &     njoined, join_off, idx, nocc_cls, iocc_cls
      integer ::
     &     opinp_temp(ngastp,2), opinv_temp(ngastp,2)
      logical ::
     &     open_close_inv, open_close_inp, open_close_u, open_close_spc
      real(8) ::
     &     cpu, sys, wall, cpu0, sys0, wall0

      logical, external ::
     &     iocc_equal
      integer, external ::
     &     idx_mel_list

      call atim_csw(cpu0,sys0,wall0)

      idx_inp = idx_mel_list(label_inp(1),op_info)
      idx_inv = idx_mel_list(label_inv(1),op_info)
      if (ninp.gt.2.or.nlist.gt.2)
     &   call quit(1,'inv_op','too many lists given!')
      if (ninp.eq.2) 
     &   idx_spc = idx_mel_list(label_inp(2),op_info)
      if (nlist.eq.2)
     &   idx_u = idx_mel_list(label_inv(2),op_info)

      if (idx_inp.lt.0.or.idx_inv.lt.0
     &    .or.nlist.eq.2.and.idx_u.lt.0
     &    .or.ninp.eq.2.and.idx_spc.lt.0) then
        write(luout,*) '"',label_inv,'" "',label_inp,'"'
        write(luout,*) idx_inv, idx_inp, idx_u, idx_spc
        call quit(1,'inv_op','label not on list')
      end if

      ! Point to the relevant operators and their associated files.
      me_inv => op_info%mel_arr(idx_inv)%mel
      me_inp => op_info%mel_arr(idx_inp)%mel
      if (nlist.eq.2) then
        me_u => op_info%mel_arr(idx_u)%mel
      else
        me_u => me_inv ! dummy, just to avoid problems
      end if
      if (ninp.eq.2) then
        me_spc => op_info%mel_arr(idx_spc)%mel
      else
        me_spc => me_inv ! dummy, just to avoid problems
      end if

      if (.not.associated(me_inv%fhand))
     &     call quit(1,'inv_op','no file handle defined for '//
     &                  trim(me_inv%label))
      if (.not.associated(me_inp%fhand))
     &     call quit(1,'inv_op','no file handle defined for '//
     &                  trim(me_inp%label))
      if (nlist.eq.2.and..not.associated(me_u%fhand))
     &     call quit(1,'inv_op','no file handle defined for '//
     &                  trim(me_u%label))

      open_close_inv = me_inv%fhand%unit.le.0
      open_close_u = .false.
      if (nlist.eq.2) open_close_u = me_u%fhand%unit.le.0
      if (ninp.eq.2) open_close_spc = me_spc%fhand%unit.le.0
      if(open_close_inv)then
        call file_open(me_inv%fhand)
      endif
      open_close_inp = me_inp%fhand%unit.le.0
      if(open_close_inp)then
        call file_open(me_inp%fhand)
      endif
      if (open_close_u) call file_open(me_u%fhand)
      if (open_close_spc) call file_open(me_spc%fhand)

      if(ntest.ge.100)then
        write(luout,*) '===================='
        write(luout,*) ' Operator inversion '
        write(luout,*) '===================='
        write(luout,*) 'To be inverted: ',trim(me_inp%label)
        write(luout,*) 'The inverse: ',trim(me_inv%label)
        if (nlist.eq.2) write(luout,*) 'Unitary mat.: ',trim(me_u%label)
        if (ninp.eq.2) write(luout,*) 'Spec. mat.: ',trim(me_spc%label)
      endif

      ! Check that the two operators have the same shape.
      ! not necessarily: in principle, the inverse is the 
      ! contravariant operator with a different shape!
      njoined = me_inp%op%njoined
      if(njoined.ne.me_inv%op%njoined
     &   .or.nlist.eq.2.and.njoined.ne.me_u%op%njoined)
     &     call quit(1,'inv_op','in and out incompatible: njoined')
      nocc_cls = me_inp%op%n_occ_cls
      if(nocc_cls.ne.me_inv%op%n_occ_cls
     &   .or.nlist.eq.2.and.nocc_cls.ne.me_u%op%n_occ_cls)
     &     call quit(1,'inv_op','in and out incompatible: nocc_cls')
      

      ! Compare the block structures of the operators. The assumption is 
      ! that the comparable blocks are in the same order. 
      opinp_temp(1:ngastp,1:2)=0
      opinv_temp(1:ngastp,1:2)=0

      do iocc_cls = 1, nocc_cls
        join_off = (iocc_cls-1)*njoined
        do idx=1,njoined
          opinp_temp(1:ngastp,1:2) =
     &         me_inp%op%ihpvca_occ(1:ngastp,1:2,join_off+idx)
          opinv_temp(1:ngastp,1:2) = 
     &         me_inv%op%ihpvca_occ(1:ngastp,1:2,join_off+idx)

          if (.not.iocc_equal(opinp_temp,me_inp%op%dagger,
     &         opinv_temp,me_inv%op%dagger)) then
            call quit(1,'inv_op','in and out incompatible: occs.')
          endif  
          if (nlist.eq.2) then
            opinv_temp(1:ngastp,1:2) =
     &           me_u%op%ihpvca_occ(1:ngastp,1:2,join_off+idx)
            if (.not.iocc_equal(opinp_temp,me_inp%op%dagger,
     &         opinv_temp,me_u%op%dagger))
     &         call quit(1,'inv_op','in and out incompatible: occs.(2)')
          end if

        enddo
      enddo

      ! NOTE:
      ! actually, we should check whether the input operator consists
      ! of diagonal blocks only (N(A)==N(C) in each of HPVX)
      ! (or that any off-diagonal block is zero)

      ! Call the actual inversion routine.
      if (mode(1:7).ne.'invsqrt') then
        if (nlist.eq.2) call warn('inv_op','Unitary matrix unavailable')
        if (mode(1:9).eq.'pseudoinv') then
          call pseudoinv(me_inp,me_inv,nocc_cls,
     &                   op_info,orb_info)
        else
          call invert(me_inp,me_inv,nocc_cls,
     &         op_info,orb_info)
        end if
      else
        write(luout,*) 'Calculating square root of inverse'
        call invsqrt(me_inp,me_inv,nocc_cls,mode(8:11).eq.'half',
     &       nlist.eq.2,me_u,ninp.eq.2,me_spc,
     &       op_info,orb_info,str_info,strmap_info)
      end if

      if(ntest.ge.1000)then
        call wrt_mel_file(luout,5,me_inv,1,
     &       me_inv%op%n_occ_cls,str_info,orb_info)
      endif

      if (open_close_u)
     &     call file_close_keep(me_u%fhand)
      if (open_close_spc)
     &     call file_close_keep(me_spc%fhand)
      if (open_close_inv)
     &     call file_close_keep(me_inv%fhand)
      if (open_close_inp)
     &     call file_close_keep(me_inp%fhand)

      call atim_csw(cpu,sys,wall)
      call prtim(luout,'time for inversion',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      return
      end

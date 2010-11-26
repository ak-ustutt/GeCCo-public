      subroutine inv_op(label_inp,label_inv,mode,
     &     op_info,orb_info,str_info,strmap_info)
*----------------------------------------------------------------------*
*     Wrapper subroutine used in the inversion of the matrix 
*     representation of an operator, op_inp.
*     The resultant is op_inv.
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

      character(*), intent(in) ::
     &     label_inp, label_inv
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
     &     me_inp, me_inv

      integer ::
     &     idx_inp, idx_inv,
     &     njoined, join_off, idx, nocc_cls, iocc_cls
      integer ::
     &     opinp_temp(ngastp,2), opinv_temp(ngastp,2)
      logical ::
     &     open_close_inv, open_close_inp
      real(8) ::
     &     cpu, sys, wall, cpu0, sys0, wall0

      logical, external ::
     &     iocc_equal
      integer, external ::
     &     idx_mel_list

      call atim_csw(cpu0,sys0,wall0)

      idx_inp = idx_mel_list(label_inp,op_info)
      idx_inv = idx_mel_list(label_inv,op_info)

      if (idx_inp.lt.0.or.idx_inv.lt.0) then
        write(luout,*) '"',trim(label_inv),'" "',trim(label_inp),'"'
        write(luout,*) idx_inv, idx_inp
        call quit(1,'inv_op','label not on list')
      end if

      ! Point to the relevant operators and their associated files.
      me_inv => op_info%mel_arr(idx_inv)%mel
      me_inp => op_info%mel_arr(idx_inp)%mel

      if (.not.associated(me_inv%fhand))
     &     call quit(1,'inv_op','no file handle defined for '//
     &                  trim(me_inv%label))
      if (.not.associated(me_inp%fhand))
     &     call quit(1,'inv_op','no file handle defined for '//
     &                  trim(me_inp%label))

      open_close_inv = me_inv%fhand%unit.le.0
      open_close_inp = me_inp%fhand%unit.le.0
      if(open_close_inv)then
        call file_open(me_inv%fhand)
      endif
      if(open_close_inp)then
        call file_open(me_inp%fhand)
      endif

      if(ntest.ge.100)then
        write(luout,*) '===================='
        write(luout,*) ' Operator inversion '
        write(luout,*) '===================='
        write(luout,*) 'To be inverted: ',trim(me_inp%label)
        write(luout,*) 'The inverse: ',trim(me_inv%label)
      endif

      ! Check that the two operators have the same shape.
      ! not necessarily: in principle, the inverse is the 
      ! contravariant operator with a different shape!
      njoined = me_inp%op%njoined
      if(njoined.ne.me_inv%op%njoined)
     &     call quit(1,'inv_op','in and out incompatible: njoined')
      nocc_cls = me_inp%op%n_occ_cls
      if(nocc_cls.ne.me_inv%op%n_occ_cls)
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

        enddo
      enddo

      ! NOTE:
      ! actually, we should check whether the input operator consists
      ! of diagonal blocks only (N(A)==N(C) in each of HPVX)
      ! (or that any off-diagonal block is zero)

      ! Call the actual inversion routine.
      if (mode(1:7).ne.'invsqrt') then
        call invert(me_inp,me_inv,nocc_cls,
     &       op_info,orb_info)
      else
        write(luout,*) 'Calculating square root of inverse'
        call invsqrt(me_inp,me_inv,nocc_cls,
     &       mode(8:11).eq.'half',mode(12:15).eq.'sgrm',
     &       op_info,orb_info,str_info,strmap_info)
      end if

      if(ntest.ge.1000)then
        call wrt_mel_file(luout,5,me_inv,1,
     &       me_inv%op%n_occ_cls,str_info,orb_info)
      endif

      if (open_close_inv)
     &     call file_close_keep(me_inv%fhand)
      if (open_close_inp)
     &     call file_close_keep(me_inp%fhand)

      call atim_csw(cpu,sys,wall)
      call prtim(luout,'time for inversion',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      return
      end

*----------------------------------------------------------------------*
      subroutine set_restr_prel(irestr_res,
     &     contr,op_info,ihpvgas,ngas)
*----------------------------------------------------------------------*
*     a preliminary fix for restrictions on final result of contraction
*     valid for standard CC only
*----------------------------------------------------------------------*

      implicit none
      
      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'

      integer, parameter ::
     &     ntest = 00

      type(contraction), intent(in) ::
     &     contr
      type(operator_info), intent(in) ::
     &     op_info
      integer, intent(in) ::
     &     ngas, ihpvgas(ngas)
      integer, intent(out) ::
     &     irestr_res(2,ngas,2,2)

      integer ::
     &     idxop, iblkop, maxex, ica, igas, igastp

      type(operator), pointer ::
     &     cur_op

      integer, external ::
     &     maxxlvl_op

      idxop = contr%idx_res
      iblkop = contr%iblk_res
c quick fix:
c      print *,'idxop,iblkop',idxop,iblkop
      if (idxop.eq.0) then
        maxex = 98
        irestr_res(1,1:ngas,1:2,1) = maxex
        irestr_res(2,1:ngas,1:2,1) = maxex
        irestr_res(1:2,1:ngas,1:2,2) = 0
      else
        cur_op => op_info%op_arr(idxop)%op
        maxex = maxxlvl_op(cur_op)
        irestr_res = cur_op%igasca_restr(1:2,1:ngas,1:2,1:2,iblkop)

      do ica = 1, 2
        do igas = 1, ngas
          igastp = ihpvgas(igas)
          if (cur_op%ihpvca_occ(igastp,ica,iblkop).eq.0) then
! ???
            irestr_res(1,igas,ica,1) = maxex+2
c            irestr_res(1,igas,ica,1) = 0
            irestr_res(2,igas,ica,1) = maxex+2
          end if
          ! check for quasi-inactive spaces:
          if (
     &       (cur_op%igasca_restr(1,igas,1,1,iblkop).eq.
     &        cur_op%igasca_restr(2,igas,1,1,iblkop).and.
     &        cur_op%igasca_restr(2,igas,1,1,iblkop).gt.0).or.
     &       (cur_op%igasca_restr(1,igas,2,1,iblkop).eq.
     &        cur_op%igasca_restr(2,igas,2,1,iblkop).and.
     &        cur_op%igasca_restr(2,igas,2,1,iblkop).gt.0)) then
            irestr_res(1,igas,ica,1) = maxex+2
            irestr_res(2,igas,ica,1) = maxex+2
          end if
        end do
      end do

      end if

      if (ntest.ge.100) then
        write(luout,*) 'restriction used for intermediates:'
        call wrt_rstr(luout,irestr_res,ngas)
      end if

      return
      end

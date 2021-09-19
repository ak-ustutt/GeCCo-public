*----------------------------------------------------------------------*
      subroutine prt_contr_export(lulog,contr,op_info)
*----------------------------------------------------------------------*
*     write export info on contraction onto unit lulog
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'

      integer, intent(in) ::
     &     lulog
      type(contraction), intent(in) ::
     &     contr
      type(operator_info), intent(in) ::
     &     op_info

      type(operator_array), pointer ::
     &     ops(:)
      type(operator), pointer ::
     &     op

      integer, parameter ::
     &     maxscr = 32
      
      character(len=maxscr) ::
     &     scr
      
      character(2) ::
     &     cdag
      character(256) ::
     &     opstr
      integer ::
     &     idx, idx1, idx2, nj, iblk

      ops => op_info%op_arr
      write(lulog,'(2x,a)') '/RESULT/'
      if (contr%idx_res.gt.0) then
        op => ops(contr%idx_res)%op
        iblk = contr%iblk_res
        nj   = op%njoined
        idx1 = (iblk-1)*nj+1
        idx2 = iblk*nj
        call occ2descr(scr,maxscr,
     &       op%ihpvca_occ(1:ngastp,1:2,idx1:idx2),nj)
        write(lulog,'(2x,a10,x,l,x,"[",a,"]")')trim(op%name),
     &       contr%dagger,trim(scr)
      else
        call quit(1,'prt_contr_export',
     &       'Do not call me for intermediate contractions.')
      end if

      write(lulog,'(2x,a,f25.14,i4,f25.14)') '/FACTOR/',
     &     contr%fac,contr%total_sign,contr%eqvl_fact
      write(lulog,'(2x,a,2i5)')
     &     '/#VERTICES/ ',contr%nvtx,contr%nsupvtx
      write(lulog,'(2x,a,15i4)') '/SVERTEX/',contr%svertex(1:contr%nvtx)
      write(lulog,'(2x,a,2i5)')
     &     '/#ARCS/ ',contr%narc,contr%nxarc

      write(lulog,'(2x,a)') '/VERTICES/'
      do idx = 1, contr%nvtx
        op => ops(contr%vertex(idx)%idx_op)%op
        iblk = contr%vertex(idx)%iblk_op
        nj   = op%njoined
        idx1 = iblk
        idx2 = iblk+nj-1
        call occ2descr(scr,maxscr,
     &       op%ihpvca_occ(1:ngastp,1:2,idx1:idx2),nj)
        write(lulog,'(x,a10,x,l,x,"[",a,"]")') trim(op%name),
     &       contr%vertex(idx)%dagger,trim(scr)
      end do
      write(lulog,'(2x,a)') '/ARCS/'
      do idx = 1, contr%narc
        call occ2descr(scr,maxscr,
     &       contr%arc(idx)%occ_cnt,1)
        write(lulog,'(7x,2i3,2x,"[",a,"]")') contr%arc(idx)%link(1:2),
     &       trim(scr)
      end do
      write(lulog,'(2x,a)') '/XARCS/'
      do idx = 1, contr%nxarc
        call occ2descr(scr,maxscr,
     &       contr%xarc(idx)%occ_cnt,1)
        write(lulog,'(7x,2i3,2x,"[",a,"]")') contr%xarc(idx)%link(1:2),
     &       trim(scr)
      end do

      write(lulog,'(2x,a)') '/CONTR_STRING/'
      idx = contr%nidx
      write(lulog,'(2x,20i4)') contr%contr_string(1:idx)%vtx
      write(lulog,'(2x,20i4)') contr%contr_string(1:idx)%ca
      write(lulog,'(2x,20i4)') contr%contr_string(1:idx)%hpvx
      write(lulog,'(2x,20l4)') contr%contr_string(1:idx)%ext
      write(lulog,'(2x,20i4)') contr%contr_string(1:idx)%cnt
      write(lulog,'(2x,20i4)') contr%contr_string(1:idx)%idx

      
      write(lulog,'(2x,a)') '/RESULT_STRING/'
      idx = contr%nxidx
      write(lulog,'(2x,20i4)') contr%result_string(1:idx)%vtx
      write(lulog,'(2x,20i4)') contr%result_string(1:idx)%ca
      write(lulog,'(2x,20i4)') contr%result_string(1:idx)%hpvx
      !write(lulog,'(2x,20l4)') contr%result_string(1:idx)%ext
      write(lulog,'(2x,20i4)') contr%result_string(1:idx)%cnt
      write(lulog,'(2x,20i4)') contr%result_string(1:idx)%idx
      
      return
      end

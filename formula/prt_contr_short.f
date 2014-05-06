*----------------------------------------------------------------------*
      subroutine prt_contr_short(lulog,idx,contr,op_info)
*----------------------------------------------------------------------*
*     write info on contraction onto unit lulog
*     one line printout version
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'

      integer, intent(in) ::
     &     lulog, idx
      type(contraction), intent(in) ::
     &     contr
      type(operator_info), intent(in) ::
     &     op_info

      type(operator_array), pointer ::
     &     ops(:)

      integer, parameter ::
     &     maxstr = 512, maxscr = 32

      character(2) ::
     &     cdag
      character(maxstr) ::
     &     str
      character(maxscr) ::
     &     scr
      integer ::
     &     ii, idxph, ipos

      ops => op_info%op_arr

      str(1:maxstr) = ' '
      ! write term number and result
      cdag = '  '
      if (contr%dagger) cdag = '^+'
      write(str,'(i6,1x,a,"(",i3,") <- ")') idx,
     &     trim(ops(contr%idx_res)%op%name)//cdag,contr%iblk_res
      ipos = 19+len_trim(ops(contr%idx_res)%op%name) + 1

      ! write factor
      write(str(ipos:),'(f10.6,1x)') contr%fac
      ipos = ipos+11

      do ii = 1, contr%nvtx
        cdag = '  '
        if (contr%vertex(ii)%dagger) then
          idxph=2
          cdag = '^+'
        end if
        if (contr%vertex(ii)%idx_op.eq.0) then
          write(str(ipos:),'("0")')
          ipos = ipos + 4
        else
          write(str(ipos:),'(a)')
     &       trim(ops(contr%vertex(ii)%idx_op)%op%name)//cdag
          ipos = ipos+len_trim(ops(contr%vertex(ii)%idx_op)%op%name)
          if (contr%vertex(ii)%dagger) ipos = ipos+2
          if (contr%vertex(ii)%iblk_op.lt.10) then
            write(str(ipos:),'("(",i1,") ")')
     &           contr%vertex(ii)%iblk_op 
            ipos = ipos+4
          else if (contr%vertex(ii)%iblk_op.lt.100) then
            write(str(ipos:),'("(",i2,") ")')
     &           contr%vertex(ii)%iblk_op 
            ipos = ipos+5
          else
            write(str(ipos:),'("(",i3,") ")')
     &         contr%vertex(ii)%iblk_op 
            ipos = ipos+6
          end if
        end if
        if (ipos.gt.maxstr)
     &       call quit(1,'prt_contr_short','string too short !')
      end do

!      write(lulog,*) trim(str)

      do ii = 1, contr%narc
        if (contr%arc(ii)%occ_cnt(1,1).lt.0) then
          ! prototype connection:
         write(str(ipos:),'("[",i1,",",i1,",?]")')
     &         contr%arc(ii)%link(1:2)
         ipos = ipos+7
        else
          call occ2descr(scr,maxscr,
     &         contr%arc(ii)%occ_cnt(1:ngastp,1),1)
          write(str(ipos:),'("[",i1,",",i1,",",a,"]")')
     &         contr%arc(ii)%link(1:2),trim(scr)
          ipos = ipos+6+len_trim(scr)
        end if
      end do
      do ii = 1, contr%nxarc
        if (contr%xarc(ii)%occ_cnt(1,1).lt.0) then
          ! prototype connection:
         write(str(ipos:),'("{",i1,",",i1,",?}")')
     &         contr%xarc(ii)%link(1:2)
         ipos = ipos+7
        else
          call occ2descr(scr,maxscr,
     &         contr%xarc(ii)%occ_cnt(1:ngastp,1),1)
          write(str(ipos:),'("{",i1,",",i1,",",a,"}")')
     &         contr%xarc(ii)%link(1:2),trim(scr)
          ipos = ipos+6+len_trim(scr)
        end if
      end do

      write(lulog,'(a)') trim(str)

      return
      end

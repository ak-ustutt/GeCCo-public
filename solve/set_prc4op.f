*----------------------------------------------------------------------*
      subroutine set_prc4op(label_prc,label_ham,
     &                      op_info,
     &                      str_info,orb_info)
*----------------------------------------------------------------------*
*     driver routine for setting up a diagonal preconditioner
*
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 000

      include 'stdunit.h'
      include 'mdef_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'ifc_memman.h'

      character(*), intent(in) ::
     &     label_ham, label_prc
      type(operator_info), intent(inout) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info

      logical ::
     &     open_close_ham,
     &     open_close_prc
      integer ::
     &     ifree, nops,
     &     idxprc, idxham
      real(8) ::
     &     cpu, sys, wall, cpu0, sys0, wall0

      type(me_list), pointer ::
     &     me_ham, me_prc
      type(filinf), pointer ::
     &     ffopham, ffopprc

      real(8), pointer ::
     &     x1dia(:)

      integer, external ::
     &     idx_mel_list

      call atim_csw(cpu0,sys0,wall0)

      ifree = mem_setmark('prc4op')
      ! this assumption is probably not too bad:
      ifree = mem_alloc_real(x1dia,2*(orb_info%ntoob+orb_info%caborb)
     &     ,'x1dia')

      idxham = idx_mel_list(label_ham,op_info)
      idxprc = idx_mel_list(label_prc,op_info)

      if (idxham.lt.0.or.idxprc.lt.0) then
        write(luout,*) '"',trim(label_ham),'" "',trim(label_prc),'"'
        write(luout,*) idxham, idxprc
        call quit(1,'set_prc4op','label not on list')
      end if

      me_ham => op_info%mel_arr(idxham)%mel
      me_prc => op_info%mel_arr(idxprc)%mel

      if (.not.associated(me_ham%fhand))
     &     call quit(1,'set_prc4op','no file handle defined for '//
     &                  trim(me_ham%label))
      if (.not.associated(me_prc%fhand))
     &     call quit(1,'set_prc4op','no file handle defined for '//
     &                  trim(me_prc%label))

      if (iprlvl.ge.1)
     &     write(luout,*) 'set up diagonal'//
     &     ' from rank 1 part of ',trim(me_ham%op%name)

      open_close_ham = me_ham%fhand%unit.le.0
      open_close_prc = me_prc%fhand%unit.le.0
      if (open_close_ham)
     &     call file_open(me_ham%fhand)
      if (open_close_prc)
     &     call file_open(me_prc%fhand)

      ! extract the fock-matrix diagonal
      call onedia_from_op(x1dia,me_ham,orb_info)

      ! set up preconditioner
      call dia4op(me_prc,x1dia,str_info,orb_info)      

      if (ntest.ge.1000) then
        call wrt_mel_file(luout,5,me_prc,1,
     &       me_prc%op%n_occ_cls,str_info,orb_info)
      end if

      if (open_close_ham)
     &     call file_close_keep(me_ham%fhand)
      if (open_close_prc)
     &     call file_close_keep(me_prc%fhand)

      ifree = mem_flushmark()

      call atim_csw(cpu,sys,wall)

      call prtim(luout,'time for diagonal ',
     &                cpu-cpu0,sys-sys0,wall-wall0)

      return
      end

*----------------------------------------------------------------------*
      subroutine spin_prj_list_drv(label,s2,op_info,
     &                      str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
* driver for spin_prj_list
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'stdunit.h'
      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'def_strmapinf.h'
      include 'ifc_memman.h'
      include 'par_opnames_gen.h'

      character(len=*), intent(in) ::
     &     label
      integer, intent(in) ::
     &     s2
      type(operator_info), intent(inout) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf), intent(in) ::
     &     strmap_info
      type(orbinf), intent(in) ::
     &     orb_info

      integer ::
     &     ifree, idx
      real(8) ::
     &     cpu, sys, wall, cpu0, sys0, wall0, xdum

      type(me_list), pointer ::
     &     me
      integer, external ::
     &     idx_mel_list

      call atim_csw(cpu0,sys0,wall0)

      ifree = mem_setmark('spin_prj_mel')

      idx = idx_mel_list(label,op_info)
      if (idx.lt.0) call quit(1,'spin_prj_list_drv',
     &                           'label not on list: '//trim(label))
      me => op_info%mel_arr(idx)%mel

      call spin_prj_list(1.0,me,me,s2,
     &     xdum,.false.,
     &     op_info,str_info,strmap_info,orb_info)

      ifree = mem_flushmark('spin_prj_mel')

      call atim_csw(cpu,sys,wall)

      call prtim(luout,'time for spin-projecting ME-list ',
     &                cpu-cpu0,sys-sys0,wall-wall0)

      return
      end

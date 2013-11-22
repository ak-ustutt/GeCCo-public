*----------------------------------------------------------------------*
      subroutine mod_op_for_ge_test(label_mel,
     &     iRdef,norb,icase,icaseF,
     &     op_info,str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*     modify operator elements (most likely the Hamiltonian)
*     for generalized extensivity test
*     iRdef(norb) defines system R, the rest is assumed system S
*     icaseF: 0 = Fock is part of E0, else = Fock is also split into R/S
*     case 1: set R and S to zero
*     case 2: only R
*     case 3: only S
*     case 4: R and S (coupling elements are 0 nontheless)
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 1000

      include 'stdunit.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'par_opnames_gen.h'
      include 'mdef_operator_info.h'
      include 'def_strmapinf.h'

      character(*), intent(in) ::
     &     label_mel
      type(operator_info), intent(in) ::
     &     op_info
      integer, intent(in) ::
     &     norb, iRdef(norb), icase, icaseF
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf) ::
     &     strmap_info
      type(orbinf), intent(in) ::
     &     orb_info

      integer ::
     &     ipri, mode, scaling, idx_mel
      real(8) ::
     &     xdum
      logical ::
     &     list_exists
      type(me_list), pointer ::
     &     mel_target

      integer, external ::
     &     idx_mel_list


      if (ntest.ge.10) then
        call write_title(lulog,wst_dbg_subr,'mod_op_for_ge_test')
        write(lulog,*) 'Rdef:  ',iRdef(1:norb)
        write(lulog,*) 'case:  ',icase
        write(lulog,*) 'icaseF:',icaseF
      end if

      if (icase.lt.1.and.icase.gt.4)
     &     call quit(1,'mod_op_for_ge_test',
     &     'illegal value for icase')

      idx_mel = idx_mel_list(label_mel,op_info)
      if (idx_mel.lt.0)
     &     call quit(1,'mod_op_for_ge_test','Label not on list: "'//
     &     trim(label_mel)//'"')

      mel_target => op_info%mel_arr(idx_mel)%mel

      call touch_file_rec(mel_target%fhand)

      call getest_mod_mel(mel_target,
     &     iRdef,norb,icase,icaseF,
     &     str_info,orb_info)

      if (ntest.ge.10.and.(.not.mel_target%op%formal)) then
        write(lulog,*)
        write(lulog,*) 'modified list: ',trim(mel_target%label)
        if (ntest.ge.10) ipri = 1
        if (ntest.ge.50) ipri = 2
        if (ntest.ge.100) ipri = 3
        if (ntest.ge.500) ipri = 4
        if (ntest.ge.1000) ipri = 5
        call wrt_mel_file(lulog,ipri,mel_target,
     &       1,mel_target%op%n_occ_cls,
     &       str_info,orb_info)
      end if

      return
      end

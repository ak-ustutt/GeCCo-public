*----------------------------------------------------------------------*
      subroutine define_me_list(label_mel,label_op,
     &     absym,casym,gamma,s2,ms,ms_fix,
     &     rec_lo,rec_hi,diag_type,gamdiag,msdiag,
     &     op_info,orb_info,str_info,strmap_info)
*----------------------------------------------------------------------*
*     define a new ME-list with label "label_mel"
*     it is assigned to operator "label_op" (must exist)
*     and has the symmetry properties as given in line 2 above
*     the range of active records is given in line 3 
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'mdef_operator_info.h'
      include 'def_orbinf.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'

      integer, intent(in) ::
     &     absym, casym, gamma, s2, ms, rec_lo, rec_hi,
     &     diag_type, gamdiag, msdiag
      logical, intent(in) ::
     &     ms_fix
      character(*), intent(in) ::
     &     label_mel, label_op
      type(operator_info), intent(inout) ::
     &     op_info
      type(orbinf), intent(in) ::
     &     orb_info
      type(strinf), intent(inout) ::
     &     str_info
      type(strmapinf), intent(inout) ::
     &     strmap_info

      integer ::
     &     idx
      type(me_list), pointer ::
     &     mel

      integer, external ::
     &     idx_mel_list, idx_oplist2

      ! the label should not exist previously ...
      if (idx_mel_list(label_mel,op_info).gt.0)
     &     call quit(1,'def_me_list',
     &     'list exists already: "'//trim(label_mel)//'"')
      call add_me_list(label_mel,op_info)
      idx = idx_mel_list(label_mel,op_info)
      mel => op_info%mel_arr(idx)%mel

      ! ... but the operator-label must be defined
      idx = idx_oplist2(label_op,op_info)
      if (idx.le.0)
     &       call quit(1,'def_me_list',
     &     'operator not found: "'//trim(label_op)//'"')
      mel%op => op_info%op_arr(idx)%op
      ! set associated list on operator:
      mel%op%assoc_list(1:mxlen_melabel) = ' '
      mel%op%assoc_list = trim(label_mel)
      call init_me_list(1,mel,orb_info)

      ! set symmetry info
      mel%absym = absym
      mel%casym = casym
      mel%gamt  = gamma
      mel%mst   = ms
      mel%s2    = s2
      mel%fix_vertex_ms = ms_fix
      mel%diag_type = diag_type
      if (diag_type.eq.1) then
        write(luout,*) 'defining diagonal blocks only!'
        mel%gamdiag = gamdiag
        mel%msdiag = msdiag
      else if (diag_type.ne.0) then
        call quit(1,'define_me_list',
     &            'requested diag_type not yet available')
      end if

      if (abs(absym).gt.1.or.abs(casym).gt.1.or.
     &     gamma.le.0.or.gamma.gt.orb_info%nsym) then
        write(luout,*) 'absym: ',absym,abs(absym).gt.1
        write(luout,*) 'casym: ',casym,abs(casym).gt.1
        write(luout,*) 'gamma: ',gamma,
     &       gamma.le.0.or.gamma.gt.orb_info%nsym
        write(luout,*)
     &       'the symmetry specifiers marked by T are not correct'
        call quit(1,'define_me_list',
     &       'incorrect symmetry specifier for ME list')
      end if

c dbg
c      print *,'name ',trim(mel%op%name)
c      print *,'ms_fix = ', mel%fix_vertex_ms
c dbg

      ! make sure that all graphs exist to address the ME-list
      call update_graphs(str_info,mel,orb_info)
      call update_strmap(str_info,strmap_info)

      ! set dimensions
      call set_op_dim2(1,mel,str_info,orb_info%nsym)

      call init_me_list(2,mel,orb_info)

      call set_op_dim2(2,mel,str_info,orb_info%nsym)
      
      ! initialize file
      call init_mel_file(mel,-1,rec_lo,rec_hi,0)

      ! update op_list array in order to set up the lookup-table
      call update_op_arr(op_info)

      return
      end

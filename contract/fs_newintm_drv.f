*----------------------------------------------------------------------*
      subroutine fs_newintm_drv(fl_item,
     &     op_info,str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*     driver for [NEW INTERMEDIATE] operation
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'multd2h.h'
      include 'mdef_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'def_orbinf.h'

      integer, parameter ::
     &     ntest = 00

      type(formula_item), intent(in) ::
     &     fl_item
      type(strinf) ::
     &     str_info
      type(strmapinf) ::
     &     strmap_info
      type(operator_info), intent(inout) ::
     &     op_info
      type(orbinf), intent(in) ::
     &     orb_info

      type(operator), pointer ::
     &     newint_info, op_new
      type(me_list), pointer ::
     &     me_new, me_ma, me_pa

      integer ::
     &     idx, itra, itra1, itra2

      integer, external ::
     &     idx_oplist2, idx_mel_list

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'fs_newintm_drv')
        write(luout,'(2x,a)') trim(fl_item%interm%name)
        write(luout,'(2x,"attribute parentage: ",a," ",a)')
     &                        trim(fl_item%parent1),
     &                        trim(fl_item%parent2)
        call print_op_occ(luout,fl_item%interm)
      end if

      newint_info => fl_item%interm

      call add_operator(newint_info%name,op_info)
      idx = idx_oplist2(newint_info%name,op_info)
      
      op_new => op_info%op_arr(idx)%op

      call clone_operator(op_new,newint_info,.false.,orb_info)

      call add_me_list(trim(newint_info%name)//'LIST',op_info)
      idx = idx_mel_list(trim(newint_info%name)//'LIST',op_info)

      me_new => op_info%mel_arr(idx)%mel

      me_new%op => op_new

      me_new%op%assoc_list(1:mxlen_melabel) = ' '
      me_new%op%assoc_list = trim(newint_info%name)//'LIST'

      call init_me_list(1,me_new,orb_info)

      ! determine list attributes
      idx = idx_oplist2(fl_item%parent1,op_info)
      if (idx.lt.1)
     &     call err_label('fs_newintm_drv','operator',fl_item%parent1)
      idx = op_info%op2list(idx)
      if (idx.lt.1)
     &     call quit(1,'fs_newintm_drv',
     &     'no list associated with "'//trim(fl_item%parent1)//'"')
      me_ma => op_info%mel_arr(idx)%mel

      ! transpositions: sign change in mst
      itra = 1
      itra1 = 1
      itra2 = 1
      if (fl_item%tra) itra = -1
      if (fl_item%tra1) itra1 = -1
      if (fl_item%tra2) itra2 = -1

      if (len_trim(fl_item%parent2).eq.0.or.
     &        trim(fl_item%parent2).eq.'-'.or.
     &        trim(fl_item%parent2).eq.'---') then
        ! single parent only
        me_new%absym = me_ma%absym
        me_new%casym = me_ma%casym
        me_new%mst   = itra*itra1*me_ma%mst  
        me_new%gamt  = me_ma%gamt 
        me_new%s2    = me_ma%s2   
      else
        idx = idx_oplist2(fl_item%parent2,op_info)
        if (idx.lt.1)
     &       call err_label('fs_newintm_drv','operator',fl_item%parent2)
        idx = op_info%op2list(idx)
        if (idx.lt.1)
     &       call quit(1,'fs_newintm_drv',
     &       'no list associated with "'//trim(fl_item%parent2)//'"')
        me_pa => op_info%mel_arr(idx)%mel
        
        me_new%absym = me_ma%absym * me_pa%absym
        me_new%casym = 0 ! no general statement possible
cmh        me_new%mst   = me_ma%mst + me_pa%mst
cmh        if (me_ma%mst.ne.0.and.me_pa%mst.ne.0)
cmh     &       call quit(1,'fs_newintm_drv',
cmh     &       'MS handling is still not correct in general!!!!')
        me_new%mst = itra*(itra1*me_ma%mst + itra2*me_pa%mst)
        me_new%gamt = multd2h(me_ma%gamt,me_pa%gamt)
        me_new%s2    = me_ma%s2
        if (me_ma%s2.ne.0.or.me_pa%s2.ne.0)
     &       call quit(1,'fs_newintm_drv',
     &       'S2 handling is not implemented!!!!')
      end if
      me_new%diag_type = 0

      ! make sure that all graphs exist to address the ME-list
      call update_graphs(str_info,me_new,orb_info)
      call update_strmap(str_info,strmap_info)

      ! set dimensions
      call set_op_dim2(1,me_new,str_info,orb_info%nsym)

      call init_me_list(2,me_new,orb_info)

      call set_op_dim2(2,me_new,str_info,orb_info%nsym)

      ! initialize file
      call init_mel_file(me_new,-1,-1,-1,0)

      ! update op_list array in order to set up the lookup-table
      call update_op_arr(op_info)

      return
      end

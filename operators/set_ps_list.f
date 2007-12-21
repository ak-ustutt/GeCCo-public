*----------------------------------------------------------------------*
      subroutine set_ps_list(list,label,
     &     absym,casym,mst,gamt,s2,
     &     str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*     set a list (outside op_info), used for scratch files
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'def_strmapinf.h'

      type(me_list), intent(inout) ::
     &     list
      integer, intent(in) ::
     &     absym, casym, mst, gamt, s2
      character(*), intent(in) ::
     &     label
      type(strinf) ::
     &     str_info
      type(strmapinf) ::
     &     strmap_info
      type(orbinf) ::
     &     orb_info

      ! to be set outside:
      ! list%op => op
      if (.not.associated(list%op))
     &     call quit(1,'set_ps_list','list%op is not set!')
      list%label(1:len_opname) = ' '
      list%label = trim(label)
      ! set associated list on operator:
      list%op%assoc_list(1:len_opname) = ' '
      list%op%assoc_list(1:len_opname) = trim(label)
      call init_me_list(1,list,orb_info)

      ! set symmetry info
      list%absym = absym
      list%casym = casym
      list%gamt  = gamt
      list%mst   = mst
      list%s2    = s2

      ! make sure that all graphs exist to address the ME-list
      call update_graphs(str_info,list,orb_info)
      call update_strmap(str_info,strmap_info)

      ! set dimensions
      call set_op_dim2(1,list,str_info,orb_info%nsym)

      call init_me_list(2,list,orb_info)

      call set_op_dim2(2,list,str_info,orb_info%nsym)
      
      ! initialize file
      call init_mel_file(list,-1,-1,-1,0)

      return
      end

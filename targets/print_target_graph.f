*----------------------------------------------------------------------*
      subroutine print_target_graph(tgt_info,only_processed)
*----------------------------------------------------------------------*
*     print target list and dependencies in the dot language, to be 
*     interpreted by programs like dot.
*     If only_processed is set to true, only precessed targets are 
*     printed (%last_mod)
*
*     Dot is a part of graphviz (http://graphviz.org/), a software 
*     to visualize graphs. A possible usage of the dot program is as
*     follow:
*
*     dot -T eps dependency_tree.gv -Gratio="1x1" > tgt_gecco.eps
*
*     Check the man page for more options. The wikipedia gives an
*     interesting explanation on the dot laguage:
*     http://en.wikipedia.org/wiki/DOT_(graph_description_language)
*
*----------------------------------------------------------------------*
      implicit none
      
      include 'mdef_target_info.h'
      include 'def_filinf.h'
      
      type(target_info), intent(in) ::
     &     tgt_info
      logical, intent(in) ::
     &     only_processed

      type(filinf) ::
     &     fdot
      
      integer ::
     &     itgt, idx, jdx, ndim
      type(target), pointer ::
     &     tgt
      
      call file_init(fdot,"dependency_tree.gv",ftyp_sq_frm,0)
      call file_open(fdot)
      
      write(fdot%unit, '("digraph dependency_tree {")')
      do itgt = 1, tgt_info%ntargets
       tgt => tgt_info%array(itgt)%tgt
       if (tgt%last_mod.ge.0.or..not.only_processed) then
        do idx = 1, tgt%n_depends_on
         write(fdot%unit,'("""",a,"""",x,"->",x,"""",a,"""",";")')
     &        trim(tgt%name), trim(tgt%depends_on(idx))
        end do
       end if
      end do
      write(fdot%unit, '("}")')
      
      call file_close_keep(fdot)
      
      return
      end

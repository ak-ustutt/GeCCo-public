*----------------------------------------------------------------------*
      subroutine optk_move_bcontr(node_move,node_append)
*----------------------------------------------------------------------*
*
*     auxiliary routine for optimize_kernel
*
*     move the binary contraction (on item at node_move)
*     and append after node_append
*     make sure that all intermediate on which node_move depends
*     are moving too (and recursively)
*
*     andreas, oct. 2012
*  
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'def_formula_item.h'
      include 'def_formula_item_list.h'
      include 'ifc_fl_aux.h'

      character(len=16), parameter ::
     &     i_am = 'optk_move_bcontr'
      integer, parameter ::
     &     ntest = 00
      integer, parameter ::
     &     maxiter = 10 000

      type(formula_item), intent(inout), target ::
     &     node_move,
     &     node_append

      type(formula_item_list), target ::
     &     fpl_to_process, fpl_to_move
      type(formula_item_list), pointer ::
     &     fpl_to_process_pnt, fpl_to_move_pnt, loop_pnt
      type(formula_item), pointer ::
     &     fl_pnt

      type(binary_contr), pointer ::
     &     bcontr
      character(len=len_opname) ::
     &     label_op
      integer ::
     &     isave, iop


      if (ntest.gt.0) then
        call write_title(luout,wst_dbg_subr,i_am)
      end if 

      ! list of binary contractions which have to be
      ! investigated for dependencies
      call init_formula_plist(fpl_to_process)
      fpl_to_process_pnt => fpl_to_process
      ! push input node here
      fpl_to_process_pnt%item => node_move

      ! list of items to move
      call init_formula_plist(fpl_to_move)
      fpl_to_move_pnt => fpl_to_move
      ! push input node here, too
      fpl_to_move_pnt%item => node_move

      loop_pnt => fpl_to_process
      isave = 0
      main_loop: do
        isave = isave+1
        if (isave.gt.maxiter) 
     &       call quit(1,i_am,'too many iterations, infinite loop?')
        ! fetch pointer from list of unprocessed items
        bcontr => loop_pnt%item%bcontr
        ! does the present contraction depend on intermediates?
        op_loop: do iop = 1, 2        
          ! ---------------------
          ! for each intermediate
          ! ---------------------
          if (iop.eq.1) label_op = bcontr%label_op1
          if (iop.eq.2) label_op = bcontr%label_op2
          if (label_op(1:5).eq.'_LTIN')
     &              call quit(1,i_am,'_LTINs not yet considered')
          if (label_op(1:5).ne.'_STIN') cycle op_loop
          ! ---------------------------------------------------------
          ! find the previous contraction that gives the intermediate
          ! ---------------------------------------------------------
          fl_pnt => find_fl_item(loop_pnt%item,
     &                 command_list=(/command_bc,command_bc_reo/),
     &                           nlist=2,
     &                           label_res=label_op,
     &                           backward=.true.)
          if (associated(fl_pnt)) then
            ! put it on list to move
            if (ntest.ge.100) then
              write(luout,*) 'for iop = ',iop
              write(luout,*) 'put on lists: definition of '//
     &                       trim(label_op)
            end if

            call new_formula_plist_entry(fpl_to_move_pnt)
            fpl_to_move_pnt => fpl_to_move_pnt%next
            fpl_to_move_pnt%item => fl_pnt
            ! put it on list to process further
            call new_formula_plist_entry(fpl_to_process_pnt)
            fpl_to_process_pnt => fpl_to_process_pnt%next
            fpl_to_process_pnt%item => fl_pnt
          else
            call quit(1,i_am,'not found (1)')
          end if

          ! ---------------------------------------
          ! find the definition of the intermediate
          ! ---------------------------------------
          fl_pnt => find_fl_item(fpl_to_move_pnt%item,
     &                           command=command_new_intermediate,
     &                           label=label_op,
     &                           backward=.true.)
          if (associated(fl_pnt)) then
            ! put it on list to move
            call new_formula_plist_entry(fpl_to_move_pnt)
            fpl_to_move_pnt => fpl_to_move_pnt%next
            fpl_to_move_pnt%item => fl_pnt
          else
            call quit(1,i_am,'not found (2)')
          end if
        end do op_loop

        ! exit?
        if (.not.associated(loop_pnt%next)) exit main_loop
        loop_pnt => loop_pnt%next
      end do main_loop

      loop_pnt => fpl_to_move

      do
c dbg
c        print *,'relink 1'
c dbg
        if (ntest.ge.100) then
          if (loop_pnt%item%command.eq.command_new_intermediate) then
            write(luout,*) 'now moving: [NEW] '//
     &                      trim(loop_pnt%item%interm%name)
          else
            write(luout,*) 'now moving: [CONTR] '//
     &                  trim(loop_pnt%item%bcontr%label_op1)//'*'//
     &                  trim(loop_pnt%item%bcontr%label_op2)//'->'//
     &                  trim(loop_pnt%item%bcontr%label_res)
          end if
        end if
        call relink_fl_node(loop_pnt%item,node_append)
        if (.not.associated(loop_pnt%next)) exit
        loop_pnt => loop_pnt%next
      end do
 
      call dealloc_formula_plist(fpl_to_process)
      call dealloc_formula_plist(fpl_to_move)

      
      return
      end


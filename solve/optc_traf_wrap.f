*-------------------------------------------------------------------------------*
      subroutine optc_traf_wrap(me_out,me_in,ctype,
     &     ftraf,label_special,nspcl,
     &     orb_info,op_info,str_info,strmap_info)
*-------------------------------------------------------------------------------*
*
*     wrapper for optc_traf, to use it as a rule
*
*     NOTE(andreas): THE ROUTINE optc_traf IS VERY UNSAFE, DO NOT PROCEED
*           USING THIS AS A "PUBLIC" RULE UNLESS THIS IS CHANGED
*           ACTUALLY, A SIMPLE "EVALUATE" IS SUFFICIENT IN A "RULE-CONTEXT"
*     
*     Pradipta, Jan. 2015
*-------------------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'mdef_operator_info.h'
      include 'def_orbinf.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'def_formula_item.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_file_array.h'
!     include 'def_me_list.h'
!     include 'def_me_list_array.h'

      character(*), intent(in) ::
     &     ctype,label_special(nspcl)
      integer, intent(in) :: 
     &     nspcl
      type(me_list), intent(in) ::
     &     me_in,me_out
      type(formula_item), intent(in) ::
     &     ftraf

      type(orbinf), intent(in) ::
     &     orb_info
      type(strinf),intent(in) ::
     &     str_info
      type(strmapinf),intent(in) ::
     &     strmap_info
      type(operator_info), intent(inout) ::
     &     op_info


      integer :: 
     &     nwfpar, irec_in, irec_out,idx,jopt,ierr,idxmel
      real(8), pointer ::
     &     xbuf1(:)
      real(8) ::
     &     xnrm
      character(32) ::
     &     label
      type(me_list_array), pointer ::
     &     me_special(:)
      type(file_array), pointer ::
     &     ffspecial(:)
      integer, external :: 
     &     idx_mel_list


      allocate(me_special(nspcl))
      allocate(ffspecial(nspcl))

      xnrm = 0.0d0
      irec_in = 1
      irec_out = 1

      nwfpar = me_in%len_op

      allocate(xbuf1(nwfpar))

      if (me_in%fhand%unit.le.0)then
        call file_open(me_in%fhand)
      endif

      if (me_out%fhand%unit.le.0)then
        call file_open(me_out%fhand)
      endif

      ierr = 0
      do idx = 1, nspcl
        jopt = idx
        idxmel = idx_mel_list(label_special(idx),op_info)
        ierr = 1
        if (idxmel.le.0) exit
        me_special(idx)%mel  => op_info%mel_arr(idxmel)%mel
        ffspecial(idx)%fhand => op_info%mel_arr(idxmel)%mel%fhand
        ierr = 2
        if (.not.associated(ffspecial(idx)%fhand)) exit
        ierr = 0
        if (me_special(idx)%mel%fhand%unit.le.0)then
          call file_open(me_special(idx)%mel%fhand)
        endif
      end do


      if (ierr.gt.0) then
        if (ierr.eq.1.or.ierr.eq.2) label = label_special(jopt)
        if (ierr.eq.1)
     &       call quit(1,'solve_leq',
     &       'did not find list '//trim(label))
        if (ierr.eq.2)
     &       call quit(1,'solve_leq',
     &       'no file associated to list '//trim(label))
      end if

      call optc_traf(me_out,irec_out,xnrm,me_in,irec_in,
     &     ftraf,ctype,me_special,nspcl,
     &     nwfpar,xbuf1,
     &     orb_info,op_info,str_info,strmap_info) 

      return
      end
    

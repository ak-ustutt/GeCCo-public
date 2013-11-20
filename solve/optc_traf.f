*----------------------------------------------------------------------*
      subroutine optc_traf(me_out,irec_out,xnrm,me_in,irec_in,
     &     ftraf,ctype,me_traf,n_me_traf,
     &     nwfpar,xbuf1,
     &     orb_info,op_info,str_info,strmap_info)
*----------------------------------------------------------------------*
*
*     transforms me_in to me_out using formula ftraf
*
*     me_out          : list of actual output-vector (record irec_out)
*     me_in           : list of actual input-vector (record irec_in)
*     xnrm            : norm of output-vector
*     ftraf           : transformation formula
*     me_traf(1)      : list (information only) of output vector in ftraf
*     me_traf(2)      : list (information only) of input vector in ftraf
*     me_traf(3,4)    : list of transformation matrices
*     me_traf(5)      : projector
*
*     andreas, oct. 2012 from matthias's routines
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'mdef_operator_info.h'
      include 'def_orbinf.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'

      integer, parameter ::
     &     ntest = 00
      character(len=9), parameter ::
     &     i_am = 'optc_traf'

      character, intent(in) ::
     &     ctype
      real(8), intent(out) ::
     &     xnrm
      integer, intent(in) ::
     &     n_me_traf, nwfpar, irec_in, irec_out
      type(me_list_array), intent(inout) ::
     &     me_traf(n_me_traf)
      type(me_list), intent(in) ::
     &     me_in,me_out
      real(8), intent(inout) ::
     &     xbuf1(*)

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
     &     idx, ioff,
     &     nblk, iblk, nj, iblkoff, jblk
      character(len_opname) ::
     &     op_in_name, op_out_name, op_trf_name,
     &     op_in_mel, op_out_mel, op_trf_mel,
     &     op_in_name2, op_out_name2

      integer, pointer ::
     &     occ(:,:,:)

      real(8) ::
     &     xdum

      integer, external ::
     &     idx_mel_list, iblk_occ

      real(8), external :: dnrm2


      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,i_am)
        write(luout,*)'me_in:  ',trim(me_in%label), ' record: ',irec_in
        write(luout,*)'me_out: ',trim(me_out%label),' record: ',irec_out
        write(luout,*)'type : ',ctype
        write(luout,*)'n_me_traf: ',n_me_traf
        do idx = 1, n_me_traf
          write(luout,*)'   ',trim(me_traf(idx)%mel%label)
        end do
c        write(luout,*) 'formula:'
c        call print_form_list(luout,ftraf,op_info)
      end if

      ! remember names of originally assigned operators
      op_in_name = trim(me_in%op%name)
      op_out_name = trim(me_out%op%name)
c      ! also remember the *present primary* list to which these 
c      ! operators were assigend (not necessarily the ME lists passed
c      ! here!)
c      op_in_mel  = trim(me_in%op%assoc_list)
c      op_out_mel = trim(me_out%op%assoc_list)
      ! also remember operator names and the lists assigned to them
      op_in_name2 = trim(me_traf(2)%mel%op%name)
      op_out_name2 = trim(me_traf(1)%mel%op%name)
      op_in_mel = trim(me_traf(2)%mel%op%assoc_list)
      op_out_mel = trim(me_traf(1)%mel%op%assoc_list)
      ! assign input list to input operator of the trafo-formula
      call assign_me_list(me_in%label,
     &                       me_traf(2)%mel%op%name,op_info)
      ! assign output list to output operator of the trafo-formula
      call assign_me_list(me_out%label,
     &                       me_traf(1)%mel%op%name,op_info)

      select case(ctype)
      case ('F','f')
        ! do nothing
      case ('B','b')
        ! backtrafo: use daggered transformation matrix
        if (n_me_traf.lt.4) call quit(1,i_am,'B: missing me_traf(4)')
        ! the error is most likely due to too few special lists in the
        ! call to, e.g., solve_leq, solve_evp
        op_trf_name = trim(me_traf(4)%mel%op%name)
        op_trf_mel  = trim(me_traf(4)%mel%op%assoc_list)
        call assign_me_list(me_traf(4)%mel%label,
     &                         me_traf(3)%mel%op%name,op_info)
      case ('P','p')
        ! use projection matrix
        if (n_me_traf.lt.5) call quit(1,i_am,'P: missing me_traf(5)')
        op_trf_name = trim(me_traf(5)%mel%op%name)
        op_trf_mel  = trim(me_traf(5)%mel%op%assoc_list)
        call assign_me_list(me_traf(5)%mel%label,
     &                         me_traf(3)%mel%op%name,op_info)
      case default
        call quit(1,'optc_traf','unknown type: '//ctype)
      end select

      call switch_mel_record(me_in,irec_in)
      call switch_mel_record(me_out,irec_out)
      ! evaluate projected vector
      call evaluate2(ftraf,.true.,.true.,
     &         op_info,str_info,strmap_info,orb_info,xnrm,.true.)
c dbg
c      print *,'after evaluate2: xnrm = ',xnrm
c dbg

      ! restore law and order: switch back association of lists and operators
      ! only restore the broken assignments in an unilateral way,
      ! as a bilateral re-assignment can mess up the previous assignments
      ! of op_*_name and op_*_mel
      call assign_op2me(me_in%label,
     &                       op_in_name,op_info)
      call assign_me2op(op_in_mel,
     &                       op_in_name2,op_info)
      call assign_op2me(me_out%label,
     &                       op_out_name,op_info)
      call assign_me2op(op_out_mel,
     &                       op_out_name2,op_info)
      if (ctype.eq.'b'.or.ctype.eq.'B') then
        call assign_me_list(me_traf(4)%mel%label,
     &                       op_trf_name,op_info)
        call assign_me_list(op_trf_mel,
     &                       op_trf_name,op_info)
      end if
      if (ctype.eq.'p'.or.ctype.eq.'P') then
        call assign_me_list(me_traf(5)%mel%label,
     &                       op_trf_name,op_info)
        call assign_me_list(op_trf_mel,
     &                       op_trf_name,op_info)
      end if

      if (ntest.ge.100) write(luout,*) 'end of optc_traf'

      return
      end

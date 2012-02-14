*----------------------------------------------------------------------*
      subroutine form_select_line(f_output,f_input,
     &                      title,label_opres,
     &                      ncmpnd,label_op,
     &                      igastp,mode_str,
     &                      op_info)
*----------------------------------------------------------------------*
*     collect those contractions on ffoutput in which two operators
*     on the list label_op(1:ncmpnd) are contracted via a line given
*     by igastp.
*     mode_str = 'delete': delete all these terms
*     mode_str = 'keep'  : delete all other terms
*     mode_str = 'no_ext': delete terms where a listed op has external
*                          lines in the given space
*     mode_str = 'ext'   : keep only terms where some listed op has
*                          external lines in the given space
*
*     matthias, nov 2009 (adopted from form_invariant)
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'def_formula.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     ncmpnd, igastp
      character(*), intent(in) ::
     &     title, mode_str,
     &     label_opres, label_op(ncmpnd)
      type(formula), intent(inout) ::
     &     f_input, f_output
      type(operator_info) ::
     &     op_info
      
      logical ::
     &     same, transpose
      integer ::
     &     icmpnd, idxres, idxop(ncmpnd), len
      character ::
     &     name*(form_maxlen_label*2)

      type(formula_item), target ::
     &     flist

      integer, external ::
     &     idx_oplist2

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,
     &                   'here speaks form_select_line')
        write(luout,*) ' f_input  = ',trim(f_input%label)
        write(luout,*) ' f_output = ',trim(f_output%label)
        write(luout,*) ' op_res  = ',trim(label_opres)
        do icmpnd = 1, ncmpnd
          write(luout,*) 'compound # ',icmpnd
          write(luout,*) ' op  = ',trim(label_op(icmpnd))
        end do
        write(luout,*) 'mode_str, igastp: ',trim(mode_str),igastp
      end if

      same = trim(f_input%label).eq.trim(f_output%label)

      ! get indices
      idxres = idx_oplist2(label_opres,op_info)
      if (idxres.lt.0)
     &     call quit(1,'form_select_line',
     &     'required operators are not yet defined? '//
     &       trim(label_opres))

      do icmpnd = 1, ncmpnd

        ! look for transposition label
        len = len_trim(label_op(icmpnd))
        transpose = (label_op(icmpnd)(len-1:len).eq.'^+') 
        if (transpose) len = len-2

        idxop(icmpnd) = idx_oplist2(label_op(icmpnd)(1:len),op_info)

        if (idxop(icmpnd).lt.0)
     &     call quit(1,'form_select_line',
     &     'required operators are not yet defined? '//
     &       label_op(icmpnd)(1:len))

        if (transpose) idxop(icmpnd) = -idxop(icmpnd)

      end do
c dbg
c      print *,idxop
c      print *,idxres
c dbg

      ! read in input formula
      call init_formula(flist)
      call read_form_list(f_input%fhand,flist,.true.)

      call select_line(flist,idxres,idxop,ncmpnd,igastp,mode_str)

      ! write result
      if (.not.same) then
        write(name,'(a,".fml")') trim(f_output%label)
        call file_init(f_output%fhand,name,ftyp_sq_unf,0)      
      end if
      f_output%comment = trim(title)
      call write_form_list(f_output%fhand,flist,title)

c dbg
c      write(luout,*)'TeX list'
c      call tex_form_list(luout,flist,op_info)
c dbg

      call dealloc_formula_list(flist)
      
      return
      end

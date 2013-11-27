*----------------------------------------------------------------------*
      subroutine form_insert_op(f_output,f_input,
     &                      title,label_opres,label_opins,
     &                      ncmpnd,label_op,
     &                      op_info)
*----------------------------------------------------------------------*
*
*     matthias, may 10
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
     &     ncmpnd
      character(*), intent(in) ::
     &     title,
     &     label_opres, label_opins, label_op(ncmpnd)
      type(formula), intent(inout) ::
     &     f_input, f_output
      type(operator_info) ::
     &     op_info
      
      logical ::
     &     ok, same, transpose, del_mode
      integer ::
     &     nterms, iarc, idx, icmpnd, idxres, idxop(ncmpnd), len,
     &     idxins, nvtx, iblk, nblk
      character ::
     &     name*(form_maxlen_label*2)
      logical, allocatable ::
     &     vtxinlist(:)

      type(formula_item), target ::
     &     flist
      type(formula_item), pointer ::
     &     fl_pnt, fl_pnt_next
      type(contraction), pointer ::
     &     contr

      integer, external ::
     &     idx_oplist2

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,
     &                   'here speaks form_insert_op')
        write(lulog,*) ' f_input  = ',trim(f_input%label)
        write(lulog,*) ' f_output = ',trim(f_output%label)
        write(lulog,*) ' op_res  = ',trim(label_opres)
        write(lulog,*) ' op_ins  = ',trim(label_opins)
        do icmpnd = 1, ncmpnd
          write(lulog,*) 'compound # ',icmpnd
          write(lulog,*) ' op  = ',trim(label_op(icmpnd))
        end do
      end if

      same = trim(f_input%label).eq.trim(f_output%label)

      ! get indices
      idxres = idx_oplist2(label_opres,op_info)
      if (idxres.lt.0)
     &     call quit(1,'form_insert_op',
     &     'required operators are not yet defined? '//
     &       trim(label_opres))
      idxins = idx_oplist2(label_opins,op_info)
      if (idxins.lt.0)
     &     call quit(1,'form_insert_op',
     &     'required operators are not yet defined? '//
     &       trim(label_opins))

      nblk = op_info%op_arr(idxins)%op%n_occ_cls

      do icmpnd = 1, ncmpnd

        ! look for transposition label
        len = len_trim(label_op(icmpnd))
        transpose = (label_op(icmpnd)(len-1:len).eq.'^+') 
        if (transpose) len = len-2

        idxop(icmpnd) = idx_oplist2(label_op(icmpnd)(1:len),op_info)

        if (idxop(icmpnd).lt.0)
     &     call quit(1,'form_insert_op',
     &     'required operators are not yet defined? '//
     &       label_op(icmpnd)(1:len))

        if (transpose) idxop(icmpnd) = -idxop(icmpnd)

      end do
c dbg
c      print *,idxop
c      print *,idxres
c      print *,idxins
c dbg

      ! read in input formula
      call init_formula(flist)
      call read_form_list(f_input%fhand,flist,.true.)

      fl_pnt => flist

      if (.not.associated(fl_pnt))
     &     call quit(1,'form_insert_op',
     &     'empty formula list? something is buggy')

      nterms = 0
      fl_loop: do

        fl_pnt_next => fl_pnt%next
        if (fl_pnt%command.eq.command_end_of_formula)
     &       exit fl_loop

        fl_pnt%target = idxres

        if (fl_pnt%command.eq.command_add_contribution) then

          contr => fl_pnt%contr

          nvtx = 0
          do iblk = nblk, 1, -1

            if (contr%nvtx.ne.nvtx) then
              nvtx = contr%nvtx
              ! find which vertex belongs to an operator from the given list
              allocate(vtxinlist(nvtx))
              vtxinlist = .false.
              cmp_loop: do idx = 1, nvtx
                do icmpnd = 1, ncmpnd
                 if (contr%vertex(idx)%idx_op.eq.abs(idxop(icmpnd)).and.
     &              (contr%vertex(idx)%dagger.eqv.(idxop(icmpnd).lt.0)))
     &           then
                   vtxinlist(idx) = .true.
                   cycle cmp_loop
                 end if
                end do
              end do cmp_loop
            end if

            ! insert operator if possible
            call contr_insert(contr,op_info,nvtx,vtxinlist,idxins,iblk,
     &                        idxres)

            if (contr%nvtx.ne.nvtx.or.iblk.eq.1) deallocate(vtxinlist)
          end do
          nterms = nterms+1
        end if

        fl_pnt => fl_pnt_next
        if (.not.associated(fl_pnt))
     &       call quit(1,'form_insert_op',
     &       'unexpected end of formula list')

      end do fl_loop

      ! write result
      if (.not.same) then
        write(name,'(a,".fml")') trim(f_output%label)
        call file_init(f_output%fhand,name,ftyp_sq_unf,0)      
      end if
      f_output%comment = trim(title)
      call write_form_list(f_output%fhand,flist,title)

c dbg
c      write(lulog,*)'TeX list'
c      call tex_form_list(lulog,flist,op_info)
c dbg

      call dealloc_formula_list(flist)

      if (ntest.ge.10) then
        write(lulog,*) 'generated terms: ',nterms
      end if
      
      return
      end

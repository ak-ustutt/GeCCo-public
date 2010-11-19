*----------------------------------------------------------------------*
      subroutine form_mod_factorization(f_output,f_input,
c     &                      title,
     &                      nterms,idxterms,
     &                      op_info)
*----------------------------------------------------------------------*
*     modify factorization (to user given values)
*     idxterms(1:nterms) = (/idx_term,n,iarc(1),iarc(2),..,iarc(n),.../)
*     where iarc(i) is the principle contraction arc defining the 
*     binary contraction (if a second arc links the same super-vertices,
*     this arc is contracted as well....)
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
     &     nterms, idxterms(nterms)
      type(formula), intent(inout) ::
     &     f_input, f_output
      type(operator_info) ::
     &     op_info
      
      logical ::
     &     same, this_term
      integer ::
     &     idx, jdx, nfac
      character ::
     &     name*(form_maxlen_label*2)

      type(formula_item), target ::
     &     flist
      type(formula_item), pointer ::
     &     fl_pnt, fl_pnt_next
      type(contraction), pointer ::
     &     contr

      integer, external ::
     &     idx_oplist2, imltlist

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,
     &       'here speaks form_mod_factorization')
        write(luout,*) ' f_input  = ',trim(f_input%label)
        write(luout,*) ' f_output = ',trim(f_output%label)
      end if

      same = trim(f_input%label).eq.trim(f_output%label)


      ! read in input formula
      call init_formula(flist)
      call read_form_list(f_input%fhand,flist,.true.)

      fl_pnt => flist

      if (.not.associated(fl_pnt))
     &     call quit(1,'form_mod_factorization',
     &     'empty formula list? something is buggy')

      idx = 0
      jdx = 1
      fl_loop: do

        if (jdx.gt.nterms) exit fl_loop
        if (nterms.lt.jdx+1 .or.
     &      nterms.lt.jdx+1+idxterms(jdx+1)) then
          write(luout,*) 'idxterms:'
          write(luout,'(x,5i4,x,5i4)') idxterms(1:nterms)
          write(luout,*) 'present index: ',jdx
          call quit(1,'form_mod_factorization',
     &         'not enough entries on idxterms')
        end if

        fl_pnt_next => fl_pnt%next
        if (fl_pnt%command.eq.command_end_of_formula)
     &       exit fl_loop

        if (fl_pnt%command.eq.command_add_contribution) then

          idx = idx+1
c dbg          
          print  *,'present term ',idx
          call prt_contr2(luout,fl_pnt%contr,op_info)        
c dbg

          this_term = idxterms(jdx).eq.idx

          if (this_term) then

            jdx = jdx+1
            nfac = idxterms(jdx)
c dbg
            print *,'nfac = ',nfac
c dbg
            call resize_contr(fl_pnt%contr,
     &                        0, !fl_pnt%contr%nvtx,
     &                        0, !fl_pnt%contr%narc,
     &                        0, !fl_pnt%contr%nxarc,
     &                        nfac)
c dbg          
          print  *,'after resize ',idx
          call prt_contr2(luout,fl_pnt%contr,op_info)        
c dbg
            fl_pnt%contr%nfac = nfac
            if (nfac.gt.0) then
              fl_pnt%contr%inffac(4,1:nfac) = idxterms(jdx+1:jdx+nfac)
              fl_pnt%contr%inffac(5,1:nfac) = idxterms(jdx+1:jdx+nfac)
            end if
            
            jdx = jdx+nfac+1

            if (iprlvl.ge.10) then
              write(luout,*) 'modifying term # ',idx
              call prt_contr2(luout,fl_pnt%contr,op_info)
            end if
            
          end if
        end if

        fl_pnt => fl_pnt_next
        if (.not.associated(fl_pnt))
     &       call quit(1,'form_mod_factorization',
     &       'unexpected end of formula list')

      end do fl_loop

      ! write result
      if (.not.same) then
        write(name,'(a,".fml")') trim(f_output%label)
        call file_init(f_output%fhand,name,ftyp_sq_unf,0)      
      end if
      f_output%comment = 'XXX'
      call write_form_list(f_output%fhand,flist,'XXX')

      call dealloc_formula_list(flist)

      return
      end

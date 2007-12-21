*----------------------------------------------------------------------*
      subroutine form_invariant(f_output,f_input,
     &                      title,label_opres,
     &                      ncmpnd,label_op,
     &                      op_info)
*----------------------------------------------------------------------*
*     collect those contractions on ffoutput that do not depend on
*     any operator on list label_op(1:ncmpnd)
*
*     andreas, april 2007 (new version: december 2007)
*
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
c      include 'def_filinf.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'
      include 'def_contraction_list.h'
      include 'def_formula.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     ncmpnd
      character(*), intent(in) ::
     &     title,
     &     label_opres, label_op(ncmpnd)
      type(formula), intent(inout) ::
     &     f_input, f_output
      type(operator_info) ::
     &     op_info
      
      type(contraction) ::
     &     contr

      logical ::
     &     ok
      integer ::
     &     luinput, luoutput, len,
     &     nterms, idum, idxinp, idx, icmpnd, idxres, idxop(ncmpnd)
      character ::
     &     name*(form_maxlen_label*2)

      integer, external ::
     &     idx_oplist2
      logical, external ::
     &     rd_contr

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'here speaks form_indep')
        write(luout,*) ' f_input  = ',trim(f_input%label)
        write(luout,*) ' f_output = ',trim(f_output%label)
        write(luout,*) ' op_res  = ',trim(label_opres)
        do icmpnd = 1, ncmpnd
          write(luout,*) 'compound # ',icmpnd
          write(luout,*) ' op  = ',trim(label_op(icmpnd))
        end do
      end if

      ! get indices
      idxres = idx_oplist2(label_opres,op_info)
      do icmpnd = 1, ncmpnd
        idxop(icmpnd) = idx_oplist2(label_op(icmpnd),op_info)
      end do
c dbg
c      print *,idxop
c      print *,idxres
c dbg
      if (idxop(1).lt.0.or.idxres.lt.0)
     &     call quit(1,'form_invariant',
     &     'required operators are not yet defined')

      write(name,'(a,".fml")') trim(f_output%label)
      call file_init(f_output%fhand,name,ftyp_sq_unf,0)      
      f_output%comment = trim(title)

      call file_open(f_input%fhand)
      call file_open(f_output%fhand)
      luinput = f_input%fhand%unit
      luoutput = f_output%fhand%unit
      rewind luinput
      rewind luoutput

      read(luinput)
      read(luinput) idum,idxinp

      len = len_trim(title)
      write(luoutput) len,trim(title)
      write(luoutput) idum,idxres

      ! signal, that still nothing is allocated
      contr%mxvtx = 0
      contr%mxarc = 0
      contr%mxfac = 0

      nterms = 0
      do while(rd_contr(luinput,contr,idxres))
        
        ok = .true.
        cmp_loop: do idx = 1, contr%nvtx
          do icmpnd = 1, ncmpnd
            if (contr%vertex(idx)%idx_op.eq.idxop(icmpnd)) then
              ok = .false.
              exit cmp_loop
            end if
          end do
        end do cmp_loop

        if (ok) then
          nterms = nterms+1
          contr%idx_res = idxres ! not completely OK
          call wrt_contr(luoutput,contr)
          if (ntest.ge.100) then
            call prt_contr2(luout,contr,op_info)
          end if
        end if

      end do

      call file_close_keep(f_output%fhand)
      call file_close_keep(f_input%fhand)

      call dealloc_contr(contr)

      if (ntest.ge.10) then
        write(luout,*) 'generated terms: ',nterms
      end if
      
      return
      end

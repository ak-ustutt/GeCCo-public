*----------------------------------------------------------------------*
      subroutine form_deriv(f_deriv,f_input,
     &                      label_deriv,title_deriv,
     &                      ncmpnd,idxder,idxmlt,idxres,
     &                      ops,nops)
*----------------------------------------------------------------------*
*     get the derivatives of all contraction on ffinput
*
*     ncmpnd: number of compounds in multi-compound operator -
*       e.g. for R12 : 'T' and 'C' are a compound op. (ncmpnd=2) 
*     idxder(1:ncmpnd) is the index of the operator with respect to which the
*     derivative has to be taken
*     idxmlt(1:ncmpnd) is the index of the operator which the derivative is 
*     multiplied with (0 if only the derivative is taken)
*     idxres is the index of the resulting operator (0, if scalar)     
*
*     andreas, end of 2006, put to gecco april 2007
*
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_contraction.h'
      include 'def_contraction_list.h'
      include 'def_formula.h'

      character, intent(in) ::
     &     title_deriv*(*),label_deriv*(*)
      type(formula), intent(inout) ::
     &     f_input, f_deriv
      integer, intent(in) ::
     &     nops, ncmpnd, idxder(ncmpnd), idxmlt(ncmpnd), idxres
      type(operator), intent(in) ::
     &     ops(nops)
      
      type(contraction) ::
     &     contr
      type(contraction_list), target ::
     &     conder
      type(contraction_list), pointer ::
     &     cur_conder

      integer ::
     &     luinput, luderiv,
     &     nterms, nder, idx, idum, idxinp, len, icmpnd
      character ::
     &     name*(form_maxlen_label*2)

      logical, external ::
     &     rd_contr

      write(name,'(a,".fml")') label_deriv
      call file_init(f_deriv%fhand,name,ftyp_sq_unf,0)      
      f_deriv%label = label_deriv
      f_deriv%comment = title_deriv

      call file_open(f_input%fhand)
      call file_open(f_deriv%fhand)
      luinput = f_input%fhand%unit
      luderiv = f_deriv%fhand%unit
      rewind luinput
      rewind luderiv

      read(luinput)
      read(luinput) idum,idxinp

      len = len_trim(title_deriv)
      write(luderiv) len,title_deriv
      write(luderiv) idum,idxres

      ! signal, that still nothing is allocated
      contr%mxvtx = 0
      contr%mxarc = 0
      contr%mxfac = 0
      nullify(conder%contr)

      nterms = 0
      do while(rd_contr(luinput,contr,idxinp))
        
        do icmpnd = 1, ncmpnd
          call contr_deriv(conder,nder,contr,ops,
     &         idxder(icmpnd),idxmlt(icmpnd),idxres)

          cur_conder => conder
          do idx = 1, nder
            nterms = nterms+1
            call wrt_contr(luderiv,cur_conder%contr)
            if (idx.lt.nder) cur_conder => cur_conder%next
          end do
        end do

      end do

      call file_close_keep(f_deriv%fhand)
      call file_close_keep(f_input%fhand)

      call dealloc_contr(contr)
      call dealloc_contr_list(conder)
      
      return
      end

*----------------------------------------------------------------------*
      subroutine form_deriv2(f_deriv,f_input,
     &                      title,
     &                      ncmpnd,
     &                      label_opres,label_opder,label_opmlt,
     &                      op_info)
*----------------------------------------------------------------------*
*     get the derivatives of all contraction on ffinput
*
*     ncmpnd: number of compounds in multi-compound operator -
*       e.g. for R12 : 'T' and 'C' are a compound op. (ncmpnd=2) 
*     idxder(1:ncmpnd) is the index of the operator with respect to which the
*     derivative has to be taken
*     idxmlt(1:ncmpnd) is the index of the operator which the derivative is 
*     multiplied with (array of 0s if only the derivative is taken)
*     idxres is the index of the resulting operator (0, if scalar)     
*
*     andreas, end of 2006, put to gecco april 2007
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
     &     label_opder(ncmpnd),label_opmlt(ncmpnd),label_opres
      type(formula), intent(inout) ::
     &     f_input, f_deriv
      type(operator_info) ::
     &     op_info

      type(contraction) ::
     &     contr
      type(contraction_list), target ::
     &     conder
      type(contraction_list), pointer ::
     &     cur_conder

      logical ::
     &     reo, transpose
c dbg
     &     ,first
c dbg
      integer ::
     &     idxder(ncmpnd), idxmlt(ncmpnd), idxres,
     &     luinput, luderiv, nvtx,
     &     nterms, nder, idx, idum, idxinp, len, icmpnd, ieqvfac
      character ::
     &     name*(form_maxlen_label*2)

      integer, pointer ::
     &     ivtx_reo(:),occ_vtx(:,:,:)
      logical, pointer ::
     &     fix_vtx(:)

      integer, external ::
     &     idx_oplist2
      logical, external ::
     &     rd_contr

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,
     &     'form_deriv2 at work')
        write(luout,*) ' f_input = ',trim(f_input%label)
        write(luout,*) ' f_deriv = ',trim(f_deriv%label)
        write(luout,*) ' op_res  = ',trim(label_opres)
        do icmpnd = 1, ncmpnd
          write(luout,*) 'compound # ',icmpnd
          write(luout,*) ' op_der  = ',trim(label_opder(icmpnd))
          write(luout,*) ' op_mlt  = ',trim(label_opmlt(icmpnd))
        end do
      end if

      ! get indices
      idxres = idx_oplist2(label_opres,op_info)
      do icmpnd = 1, ncmpnd
        ! look for transposition label
        len = len_trim(label_opder(icmpnd))
        transpose = (label_opder(icmpnd)(len-1:len).eq.'^+') 
        if (transpose) len = len-2

        idxder(icmpnd) = idx_oplist2(label_opder(icmpnd)(1:len),op_info)

        if (idxder(icmpnd).lt.0)
     &     call quit(1,'form_deriv2',
     &     'required operators are not yet defined? '//
     &       label_opder(icmpnd)(1:len))

        if (transpose) idxder(icmpnd) = -idxder(icmpnd)

      end do
      if (len_trim(label_opmlt(1)).eq.0) then
        if (ncmpnd.gt.1)
     &       call quit(1,'form_deriv',
     &       'no multi-compound derivative without multiplication')
        idxmlt(1) = 0
      else
        do icmpnd = 1, ncmpnd
          idxmlt(icmpnd) = idx_oplist2(label_opmlt(icmpnd),op_info)          
          if (idxmlt(icmpnd).lt.0)
     &     call quit(1,'form_deriv2',
     &         'required operators are not yet defined? '//
     &         label_opmlt(icmpnd)(1:len))
        end do
      end if
      if (idxres.lt.0) then
        write(luout,*) 'idxder:', idxder
        write(luout,*) 'idxmlt:', idxmlt
        write(luout,*) 'idxres:', idxres
        call print_op_info(luout,'ops',op_info)
        call quit(1,'form_deriv2',
     &     'required operators are not yet defined')
      end if

      write(name,'(a,".fml")') trim(f_deriv%label)
      call file_init(f_deriv%fhand,name,ftyp_sq_unf,0)      
      f_deriv%comment = trim(title)

      call file_open(f_input%fhand)
      call file_open(f_deriv%fhand)
      luinput = f_input%fhand%unit
      luderiv = f_deriv%fhand%unit
      rewind luinput
      rewind luderiv

      read(luinput)
      read(luinput) idum,idxinp

      len = len_trim(title)
      write(luderiv) len,trim(title)
      write(luderiv) idum,idxres

      ! signal, that still nothing is allocated
      call init_contr(contr)
      nullify(conder%contr)

c dbg
c      print *,'allow?',rd_contr(luinput,contr,idxinp)
c dbg

      nterms = 0
      do while(rd_contr(luinput,contr,idxinp))

        do icmpnd = 1, ncmpnd
          call contr_deriv3(conder,nder,contr,op_info,
     &         idxder(icmpnd),idxmlt(icmpnd),idxres)

          cur_conder => conder
          do idx = 1, nder
            nterms = nterms+1

            ! enforce law and order
            nvtx = cur_conder%contr%nvtx
            allocate(ivtx_reo(nvtx),fix_vtx(nvtx),
     &           occ_vtx(ngastp,2,nvtx))
            fix_vtx = .true.    ! "fix" all vertices -> ieqvfac will be 1
            call occvtx4contr(1,occ_vtx,cur_conder%contr,op_info)

            call topo_contr2(ieqvfac,reo,ivtx_reo,
     &           cur_conder%contr,fix_vtx,op_info)
c            call topo_contr(ieqvfac,reo,ivtx_reo,
c     &           cur_conder%contr,occ_vtx,fix_vtx)
            ! ieqvfac is ignored
            call canon_contr(cur_conder%contr,reo,ivtx_reo)
            deallocate(ivtx_reo,fix_vtx,occ_vtx)

            call wrt_contr(luderiv,cur_conder%contr)

            if(ntest.ge.100)then
              call prt_contr2(luout,cur_conder%contr,op_info)
            endif

c dbg
c            first = .true.
c            call tex_contr(luout,first,cur_conder%contr,op_info)
c            if(first) first = .false.
c dbg

            if (idx.lt.nder) cur_conder => cur_conder%next
          end do

          call dealloc_contr_list(conder)
        end do

      end do
c dbg
c      print *,'nterms = ',nterms
c dbg

      call file_close_keep(f_deriv%fhand)
      call file_close_keep(f_input%fhand)

      call dealloc_contr(contr)
      
      return
      end

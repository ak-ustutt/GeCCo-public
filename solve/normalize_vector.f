*-----------------------------------------------------------*
      subroutine normalize_vector(
     &     me_vec,me_met,use_s,nopt,mode)
*
*subroutine for normalizing vector/group of vectors and then 
*scaling all the vector all the vector in the list accordingly
*
*me_vec: The list of vectors need to be scaled by the norm.
*me_met: Metric-vector product corresponding to the vectors.
*use_s: A logical to specify the requirement of overlap 
*       matrix for all the vectors.
*nopt: Total number of vectors in the list
*mode: specify the vector(s) for whom the norm is calculated.
*      index of the vector in the list. If one wants the norm
*      to be calculated using all the vectors, then mode 
*      should be greater than nopt.
* written by Pradipta and Yuri, April, 2015
*-----------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'mdef_operator_info.h'
      include 'def_optimize_info.h'
      include 'ifc_memman.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     nopt,mode
      logical ::
     &     use_s(nopt)
      type(me_list_array), intent(in) ::
     &     me_vec(*), me_met(*)


      integer ::
     &     ifree, nwfpar(nopt), iopt,irecv,irecm
      real(8) ::
     &     xnrm
      real(8), pointer ::   
     &     xbuf1(:), xbuf2(:)

      real(8), external ::
     &     ddot, da_ddot

      if (ntest.ge.100) then
        write(lulog,*) '-------------------------'
        write(lulog,*) ' normalize_vector at work '
        write(lulog,*) '-------------------------'
      end if

      if (mode.le.0) 
     &       call quit(1,'normalize_vector',
     &       'mode should be greater than zero')

      do iopt = 1, nopt
        nwfpar(iopt) = me_vec(iopt)%mel%len_op
      enddo

      ifree = mem_setmark('normalize_vector')
      ifree = mem_alloc_real(xbuf1,maxval(nwfpar(1:nopt)),'xbuf1')
      ifree = mem_alloc_real(xbuf2,maxval(nwfpar(1:nopt)),'xbuf2')

      xnrm = 0d0
      if (mode.le.nopt) then
        iopt = mode
        irecv = me_vec(iopt)%mel%fhand%current_record
        if (use_s(iopt)) then
          irecm = me_met(iopt)%mel%fhand%current_record
          call vec_from_da(me_vec(iopt)%mel%fhand,irecv,
     &                     xbuf1,nwfpar(iopt))
          call vec_from_da(me_met(iopt)%mel%fhand,irecm,
     &                     xbuf2,nwfpar(iopt))
          xnrm = ddot(nwfpar(iopt),xbuf1,1,xbuf2,1)
        else
          call vec_from_da(me_vec(iopt)%mel%fhand,1,
     &                     xbuf1,nwfpar(iopt))
          xnrm = ddot(nwfpar(iopt),xbuf1,1,xbuf1,1)
        end if
      else 
         do iopt = 1, nopt 
           irecv = me_vec(iopt)%mel%fhand%current_record
           if (use_s(iopt)) then
             irecm = me_met(iopt)%mel%fhand%current_record
             call vec_from_da(me_vec(iopt)%mel%fhand,irecv,
     &                        xbuf1,nwfpar(iopt))
             call vec_from_da(me_met(iopt)%mel%fhand,irecm,
     &                        xbuf2,nwfpar(iopt))
             xnrm = xnrm + ddot(nwfpar(iopt),xbuf1,1,xbuf2,1)
           else
             xnrm = xnrm + da_ddot(me_vec(iopt)%mel%fhand,irecv,
     &                      me_vec(iopt)%mel%fhand,irecv,
     &                      nwfpar(iopt),
     &                      xbuf1,xbuf2,
     &                      nwfpar(iopt))          
           end if
        end do
      end if

      xnrm = sqrt(xnrm)
     
      print*, 'norm for mode',mode,'is:',xnrm

      if (xnrm.gt.1d-12) then
        do iopt = 1, nopt
          irecv = me_vec(iopt)%mel%fhand%current_record
          call da_sccpvec(me_vec(iopt)%mel%fhand,irecv,
     &                    me_vec(iopt)%mel%fhand,irecv,
     &                    1d0/xnrm,nwfpar(iopt),
     &                    xbuf1,nwfpar(iopt))
         end do 
      else 
        call quit(1,'normalize_vector','norm is zero')
      end if
 
      ifree = mem_flushmark()

      return

      end

*-----------------------------------------------------------*
      subroutine normalize_vector_wrap(
     &     label_mel,label_met,
     &     nopt,mode,op_info)
*-----------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'mdef_operator_info.h'
      include 'def_file_array.h'

      character(*), intent(in) ::
     &     label_mel(nopt),label_met(nopt)
      integer, intent(in) ::
     &     nopt, mode
      type(operator_info), intent(inout) ::
     &     op_info

      type(me_list_array), pointer ::
     &     me_opt(:),me_met(:)
      character(len_opname) ::
     &     label
      integer ::
     &     iopt, jopt, ierr, idxmel
      logical ::
     &     use_s(nopt)

      integer, external ::
     &     idx_mel_list

      allocate(me_opt(nopt),me_met(nopt))

      do iopt = 1, nopt
        ! pointer array for operators:
        ierr = 1
        jopt = iopt
        idxmel = idx_mel_list(label_mel(iopt),op_info)
        if (idxmel.le.0) exit
        ierr = 2
        me_opt(iopt)%mel   => op_info%mel_arr(idxmel)%mel
        jopt = iopt
        idxmel = idx_mel_list(label_met(iopt),op_info)
        if (idxmel.le.0) exit
        ierr = 0
        me_met(iopt)%mel   => op_info%mel_arr(idxmel)%mel
        use_s(iopt) = label_mel(iopt).eq.label_met(iopt)
      end do

      ! error handling
      if (ierr.gt.0) then
        if (ierr.eq.1) label = label_mel(jopt)
        if (ierr.eq.2) label = label_met(jopt)
        call quit(1,'solve_evp',
     &           'did not find list '//trim(label))
      end if

      call normalize_vector(
     &     me_opt,me_met,use_s,nopt,mode)

      return

      end

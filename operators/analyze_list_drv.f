      subroutine analyze_list_drv(label_cov,ncov,label_contrv,ncontrv,
     &     mode,records,orb_info,str_info,op_info)

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_orbinf.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'mdef_operator_info.h'

      integer, parameter ::
     &     ntest = 00
      character(len=16), parameter ::
     &     i_am='analyze_list_drv'

      integer, intent(in) ::
     &     ncov, ncontrv, records(2)
      character(len=*), intent(in) ::
     &     label_cov(ncov), label_contrv(ncontrv), mode

      type(orbinf), intent(inout) ::
     &     orb_info
      type(strinf), intent(inout) ::
     &     str_info
      type(operator_info), intent(inout) ::
     &     op_info

      integer ::
     &     idx, jdx, ierr, idxmel, nrec
      character(len=mxlen_melabel) ::
     &     label

      type(me_list_array), pointer ::
     &     me_cov(:), me_contrv(:)

      integer, external ::
     &     idx_mel_list

      if (ntest.ge.10) then
        call write_title(lulog,wst_dbg_subr,i_am)
        write(lulog,*) ' mode = "',trim(mode),'"'
      end if

      if (ncov.ne.ncontrv) then
        write(luout,*) 'ncov,ncontrv: ',ncov,ncontrv
        call quit(0,i_am,
     &      'Number of covariant and contravariant list must be equal!')
      end if

      allocate(me_cov(ncov),me_contrv(ncontrv))

      ! identify all lists
      do idx = 1, ncov
        jdx = idx
        ierr = 1
        idxmel = idx_mel_list(label_cov(idx),op_info)
        if (idxmel.le.0) exit
        me_cov(idx)%mel => op_info%mel_arr(idxmel)%mel
        ierr = 2
        idxmel = idx_mel_list(label_contrv(idx),op_info)
        if (idxmel.le.0) exit
        me_contrv(idx)%mel => op_info%mel_arr(idxmel)%mel
        ierr = 0
      end do

      if (ierr.ne.0) then
        if (ierr.eq.1) label = label_cov(jdx)
        if (ierr.eq.2) label = label_contrv(jdx)
        call quit(1,i_am,'List not found: '//trim(label))
      end if

      ! inquire number of records
      nrec = me_cov(1)%mel%fhand%active_records(2)
      do idx = 1, ncov
        jdx = idx
        ierr = 1
        if (me_cov(idx)%mel%fhand%active_records(2).ne.nrec) exit
        ierr = 0
      end do

      if (ierr.ne.0) then
        if (ierr.eq.1) label = label_cov(jdx)
        if (ierr.eq.2) label = label_contrv(jdx)
        call quit(1,i_am,
     &       'Incompatible number of records: '//trim(label))
      end if

      ! loop over all records
      do idx = 1, nrec
        do jdx = 1, ncov
          call switch_mel_record(me_cov(jdx)%mel,idx)
          call switch_mel_record(me_contrv(jdx)%mel,idx)
        end do

        call analyze_list_core(me_cov,me_contrv,ncov,idx,mode,
     &       orb_info,str_info)

      end do

      deallocate(me_cov,me_contrv)

      end

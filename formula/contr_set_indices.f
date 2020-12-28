*-----------------------------------------------------------------------*
      subroutine contr_set_indices(contr,op_info)
*-----------------------------------------------------------------------*
*
*     Assign explicit index labels to all arcs of a contraction and
*     determine the corresponding total sign. Some heuristics are
*     employed to arrive at typical topologies (and sign) for standard
*     diagrams. The signs do not fit for (present) GeCCo internal
*     contractions but are suited for generation of input for external
*     tensor engines (in particular ITF).
*
*     input/output: contr (the index information is held in a special
*                          set of records of this structure
*                          see explanations for type string_element
*                          in def_contraction.h)
*     auxiliary input: op_info
*
*     andreas, dec. 2020
*
*-----------------------------------------------------------------------*
      Implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'

      integer, parameter ::
     &     ntest = 100
      character(len=17), parameter ::
     &     i_am = 'contr_set_indices'
      
      type(contraction), intent(inout) ::
     &     contr
      type(operator_info), intent(in) ::
     &     op_info

      integer ::
     &     nvtx, narc, nxarc,
     &     iarc, ivtx, ivtx1, ivtx2, ivtxleft, ivtxright,
     &     ihpvx, ica, jhpvx, jca,
     &     iseq, jseq, kseq,
     &     nencl, count, ncnt, nex, hpvxend, sign
      integer ::
     &     nidx, nxidx, idx, jdx, ii, jj, inc, ind, ist,
     &     ioff1, nidx1, ioff2, nidx2,
     &     idx1, idx2, jdx1, jdx2, nexta, nextc, icnt,
     &     npair, nextra, ipair, iextra, nxvtx, ntransp,
     &     ixvtx_last
      integer ::
     &     idxoff(ngastp)

      type(string_element), pointer ::
     &     string(:), result(:)
      
      integer, pointer ::
     &     occ_vtx(:,:,:), occ_cnt(:,:), vtxseq(:,:)
      integer, pointer ::
     &     ivtxoff(:), nidxvtx(:), ixvtxoff(:), nidxxvtx(:)
      type(cntr_arc), pointer ::
     &     arc(:), xarc(:)

      if (ntest.gt.0) call write_title(lulog,wst_dbg_subr,i_am)

      if (ntest.ge.100) then
        write(lulog,*) 'contraction: '
        call prt_contr2(lulog,contr,op_info)
      end if

      nvtx = contr%nvtx
      narc = contr%narc
      nxarc = contr%nxarc

      allocate(occ_vtx(ngastp,2,nvtx),vtxseq(2,2*ngastp))

      ! get occupations of operators (transposed operators are handled here)
      call occvtx4contr(1,occ_vtx,contr,op_info)

!     set sequence for looping over vertex:
!     C: H P V X, A: X V P H
      iseq=0
      do ihpvx = 1, ngastp
        iseq=iseq+1
        vtxseq(1,iseq) = 1
        vtxseq(2,iseq) = ihpvx
      end do
      do ihpvx = ngastp, 1, -1
        iseq=iseq+1
        vtxseq(1,iseq) = 2
        vtxseq(2,iseq) = ihpvx
      end do

! total length of string
      nidx = sum(occ_vtx(1:ngastp,1:2,1:nvtx))

      allocate(string(nidx),ivtxoff(nvtx),nidxvtx(nvtx))

      idx = 0
      do ivtx = 1, nvtx
        ivtxoff(ivtx) = idx
        nidxvtx(ivtx) = 0
        do iseq = 1, 2*ngastp
          ica = vtxseq(1,iseq)
          ihpvx = vtxseq(2,iseq)
!          if (occ_vtx(ihpvx,ica,ivtx)==0) cycle
          do ii = 1, occ_vtx(ihpvx,ica,ivtx)
            idx = idx+1
            nidxvtx(ivtx) = nidxvtx(ivtx)+1
            string(idx)%vtx = ivtx
            string(idx)%ca = ica
            string(idx)%hpvx = ihpvx
            string(idx)%cnt = 0
            string(idx)%idx = 0
            string(idx)%ext = .false.
            string(idx)%del = .false.
          end do
        end do
      end do

      if (ntest.ge.100) then
        call print_string(string,nidx,nvtx)
      end if

!     go through contractions and assign:
!     'c' or 'x' for contraction or external line
!     number of contraction (increasing index)
!     index number (increasing index for each h p v x)

      arc => contr%arc
      do iarc = 1, narc

        ! get contraction info
        ivtx1 = arc(iarc)%link(1)
        ivtx2 = arc(iarc)%link(2)
        occ_cnt => arc(iarc)%occ_cnt

! distribute info to string
        ioff1 = ivtxoff(ivtx1)
        nidx1 = nidxvtx(ivtx1)
        ioff2 = ivtxoff(ivtx2)
        nidx2 = nidxvtx(ivtx2)
        do ica = 1, 2
          do ihpvx = 1, ngastp
            if (occ_cnt(ihpvx,ica)==0) cycle
            ncnt = occ_cnt(ihpvx,ica)
! search first vertex
            if (ica==1) then
              ist = ioff1+1
              ind = ioff1+nidx1
              inc = +1
            else
              ist = ioff1+nidx1
              ind = ioff1+1
              inc = -1
            end if
            do idx = ist, ind, inc
              if ( string(idx)%ca==ica.and.
     &             string(idx)%hpvx==ihpvx.and.
     &             string(idx)%cnt==0 ) then
                string(idx)%cnt = iarc
                ncnt = ncnt-1
              end if
              if (ncnt==0) exit
            end do
            ncnt = occ_cnt(ihpvx,ica)
! search second vertex
            if (ica==2) then
              ist = ioff2+1
              ind = ioff2+nidx2
              inc = +1
            else
              ist = ioff2+nidx2
              ind = ioff2+1
              inc = -1
            end if            
            do idx = ist, ind, inc
              if ( string(idx)%ca==3-ica.and.  ! search for A if CNT has a C and vice versa
     &             string(idx)%hpvx==ihpvx.and.
     &             string(idx)%cnt==0 ) then
                string(idx)%cnt = iarc
                ncnt = ncnt-1
              end if
              if (ncnt==0) exit
            end do
          end do
        end do
      end do

      !call print_string_cnt(string,nidx,nxvtx)

      if (nxarc>0) xarc => contr%xarc
      do iarc = 1, nxarc

        ! get contraction info
        ivtx1 = xarc(iarc)%link(1)
        ivtx2 = xarc(iarc)%link(2)  ! vertex of result
        occ_cnt => xarc(iarc)%occ_cnt

! distribute info to string
        ioff1 = ivtxoff(ivtx1)
        nidx1 = nidxvtx(ivtx1)
        do ica = 1, 2
          do ihpvx = 1, ngastp
            if (occ_cnt(ihpvx,ica)==0) cycle
            ncnt = occ_cnt(ihpvx,ica)
            
! search first vertex
            if (ica==1) then
              ist = ioff1+1
              ind = ioff1+nidx1
              inc = +1
            else
              ist = ioff1+nidx1
              ind = ioff1+1
              inc = -1
            end if
            do idx = ist, ind, inc
              if ( string(idx)%ca==ica.and.
     &             string(idx)%hpvx==ihpvx.and.
     &             string(idx)%cnt==0 ) then
                string(idx)%cnt = ivtx2
                string(idx)%ext = .true.
                ncnt = ncnt-1
              end if
              if (ncnt==0) exit
            end do
          end do
        end do
      end do

      if (ntest.ge.1000) then
        call print_string_cnt(string,nidx,nvtx)
      end if

!     now: distribute C A pairs over external indices
      nextc = 0
      nexta = 0
      do idx = 1, nidx
        if (string(idx)%ext) then
          if (string(idx)%ca==1) nextc = nextc+1
          if (string(idx)%ca==2) nexta = nexta+1
        end if
      end do

      npair = min(nextc,nexta)
      nextra = max(nextc,nexta)-npair
!     if (nextra>0)
!     & call quit(1,i_am,'have a look what to do for unpaired externals')

      idxoff(1:ngastp) = 0
      
      pair_loop: do ipair = 1, npair
! 1. try to distribute within single vertices
        do ivtx = 1, nvtx
          ioff1 = ivtxoff(ivtx)
          nidx1 = nidxvtx(ivtx)
          do idx = 1, nidx1/2
            idx1 = ioff1 + idx
            idx2 = ioff1 + nidx1 + 1 - idx
            ! look for pair of externals which still have no index
            if ( string(idx1)%ext.and.string(idx2)%ext.and.
     &           string(idx1)%idx==0.and.string(idx2)%idx==0 ) then
              ihpvx = string(idx1)%hpvx
              string(idx1)%idx = idxoff(ihpvx)+1
              idxoff(ihpvx) = idxoff(ihpvx)+1
              ihpvx = string(idx2)%hpvx
              string(idx2)%idx = idxoff(ihpvx)+1
              idxoff(ihpvx) = idxoff(ihpvx)+1
              cycle pair_loop
            end if
          end do
        end do
!     2. try do distribute over entire string
        do idx1 = 1, nidx
          ! look for external which still has no index
          if (string(idx1)%ext.and.string(idx1)%idx==0) then
            do idx2 = nidx, idx1+1, -1
              if (string(idx2)%ext.and.string(idx2)%idx==0) then
                ihpvx = string(idx1)%hpvx
                string(idx1)%idx = idxoff(ihpvx)+1
                idxoff(ihpvx) = idxoff(ihpvx)+1
                ihpvx = string(idx2)%hpvx
                string(idx2)%idx = idxoff(ihpvx)+1
                idxoff(ihpvx) = idxoff(ihpvx)+1
                cycle pair_loop
              end if
            end do
          end if
        end do
          
      end do pair_loop

      extra_loop: do iextra = 1, nextra
        do idx1 = 1, nidx
          if (string(idx1)%ext.and.string(idx1)%idx==0) then
            ihpvx = string(idx1)%hpvx
            string(idx1)%idx = idxoff(ihpvx)+1
            idxoff(ihpvx) = idxoff(ihpvx)+1
            cycle extra_loop
          end if
        end do
      end do extra_loop

!     align p/h indices
      idx = max(idxoff(IHOLE),idxoff(IPART))
      idxoff(IHOLE) = idx
      idxoff(IPART) = idx
      
!     now distribute contraction indices
      do iarc = 1, narc
!     get contraction info
        ivtx1 = arc(iarc)%link(1)
        ivtx2 = arc(iarc)%link(2)
        occ_cnt => arc(iarc)%occ_cnt

        ncnt = sum(occ_cnt(1:ngastp,1:2))
     
!     distribute info to string
        ioff1 = ivtxoff(ivtx1)
        nidx1 = nidxvtx(ivtx1)
        ioff2 = ivtxoff(ivtx2)
        nidx2 = nidxvtx(ivtx2)
       
        cnt_loop: do icnt = 1, ncnt
          do idx = 1, nidx1
            idx1 = ioff1 + idx
            if ( string(idx1)%cnt==iarc.and.
     &           string(idx1)%idx==0) then
              ihpvx = string(idx1)%hpvx
              ica = 3-string(idx1)%ca    ! exchange C<->A for the partner index of contraction
              do jdx = 1, nidx2
                idx2 = ioff2 + nidx2 + 1 - jdx  ! search backwards on partner vertex
                if ( string(idx2)%cnt==iarc.and.
     &               string(idx2)%hpvx == ihpvx .and.
     &               string(idx2)%ca == ica .and.
     &               string(idx2)%idx==0 ) then
                  idxoff(ihpvx) = idxoff(ihpvx)+1
                  string(idx1)%idx = idxoff(ihpvx)
                  string(idx2)%idx = idxoff(ihpvx)
                  cycle cnt_loop
                end if
              end do
            end if
          end do 
          
        end do cnt_loop
        
      end do

!     test: enforce some standard order here:
!     as this is our defined starting point, we ignore the sign
      sign = 0
      do ivtx = 1, nvtx
        ioff1 = ivtxoff(ivtx)
        nidx1 = nidxvtx(ivtx)
        if (nidx1==0) cycle
        call sort_string(sign,string(ioff1+1),nidx1)
      end do

      if (ntest.ge.100) then
        write(lulog,*) 'generated contraction index string:'
        call print_string(string,nidx,nvtx)
        call print_string_cnt(string,nidx,nvtx)
        call print_string_idx(string,nidx,nvtx)
      end if
        
!     now: determine the sign
      sign = 1
      do idx = 1, nidx
        if (.not.string(idx)%ext.and..not.string(idx)%del) then
          ihpvx = string(idx)%hpvx  ! get index to look for
          idx1 = string(idx)%idx
          string(idx)%del = .true.  ! inactivate this one
          nencl = 0
          do jdx = idx+1, nidx  ! search for partner index and count on the way
            if (string(jdx)%del) cycle ! skip deleted vertices
            if (string(jdx)%hpvx==ihpvx.and.string(jdx)%idx==idx1) then
              string(jdx)%del = .true. ! inactivate
              if (mod(nencl,2)>0) sign = -sign
              exit
            else
              nencl = nencl+1
            end if
          end do
          if (ntest.ge.1000) then
            write(lulog,'(1x,a,i4,a,i4)') 'cnt removed, nencl = ',nencl,
     &           ' -> new sign = ',sign
            call print_string_idx(string,nidx,nvtx)
          end if
        end if
      end do

      ! store info on diagram labels
      contr%index_info = .true.
      allocate(contr%contr_string(nidx))
      contr%nidx = nidx
      string(1:nidx)%del = .false.  ! set all active again
      contr%contr_string = string
      
! assemble result
      if (nxarc>0) then
        nxidx = 0
        nxvtx = 0
        do idx = 1, nidx
          if (string(idx)%ext) then
            nxidx = nxidx+1
            nxvtx = max(nxvtx,string(idx)%cnt)
          end if
        end do

        allocate(result(nxidx),ixvtxoff(nxvtx),nidxxvtx(nxvtx))

        ixvtxoff(1:nxvtx) = 0
        nidxxvtx(1:nxvtx) = 0
        
        idx1 = 0
        ixvtx_last = 0
        do idx = 1, nidx
          if (string(idx)%ext) then
            idx1 = idx1+1
            result(idx1)%vtx  = string(idx)%cnt
            result(idx1)%ca   = string(idx)%ca
            result(idx1)%hpvx = string(idx)%hpvx
            result(idx1)%cnt  = string(idx)%cnt
            result(idx1)%idx  = string(idx)%idx
            result(idx1)%ext = .true.
            result(idx1)%del = .false.
            do while (ixvtx_last<result(idx1)%vtx) 
              ixvtx_last = ixvtx_last+1
              ixvtxoff(ixvtx_last) = idx1-1
            end do
            nidxxvtx(ixvtx_last) = nidxxvtx(ixvtx_last)+1
          end if
        end do

        if (ntest.ge.1000) then
          write(lulog,*) 'initial result: sign = ',sign
          call print_string(result,nxidx,nxvtx)
          call print_string_cnt(result,nxidx,nxvtx)
          call print_string_idx(result,nxidx,nxvtx)
        end if

! run over vertices and sort creations and annihilations to proper order
        do ivtx = 1, nxvtx
          ioff1 = ixvtxoff(ivtx)
          nidx1 = nidxxvtx(ivtx)
          call sort_string(sign,result(ioff1+1),nidx1)
        end do
        
        if (ntest.ge.100) then
          write(lulog,*) 'final result: sign = ',sign
          call print_string(result,nxidx,nxvtx)
          call print_string_cnt(result,nxidx,nxvtx)
          call print_string_idx(result,nxidx,nxvtx)
        end if

        ! store information on result string
        allocate(contr%result_string(nxidx))
        contr%nxidx = nxidx
        contr%result_string = result
        
        deallocate(ixvtxoff,nidxxvtx)
        deallocate(result)
        
      end if

      if (ntest.ge.100) then
        write(lulog,*) 'final sign = ',sign
      end if
      ! store sign
      contr%total_sign = sign
      
      deallocate(string,ivtxoff,nidxvtx)
      deallocate(occ_vtx,vtxseq)

      return
      
      end subroutine

      
      

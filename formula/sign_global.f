*----------------------------------------------------------------------*
      integer function sign_global(arc_seq,contr,op_info)
*----------------------------------------------------------------------*
*     determine global sign of contraction (for printing purposes only)
*     arc_seq determines the sequence to consider the contraction arcs
*     (phase factors may depend on that)
*     In the sequence of contraction, the outer indices of each
*     vertex are contracted:
*
*        standard sequence: VTX(CH,CP,CV,CX;AX,AV,AP,AH)
*
*     the sign (hopefully) conforms with the automatic indexing done
*     by tex_contr
*
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     arc_seq(*)
      type(contraction), intent(in) ::
     &     contr
      type(operator_info), intent(in) ::
     &     op_info

      integer ::
     &     nvtx, narc, nxarc,
     &     count, iarc0, iarc, ivtx, ivtx1, ivtx2,
     &     hpvx, hpvx2, nencl, iocc, n_vtx1_arem, n_vtx2_crem, nrem

      integer, pointer ::
     &     occ_vtx(:,:,:), occ_cnt(:,:)
      type(cntr_arc), pointer ::
     &     arc(:)

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'sign_global')
      end if

      nvtx = contr%nvtx
      narc = contr%narc
      nxarc = contr%nxarc

      arc => contr%arc

      allocate(occ_vtx(ngastp,2,nvtx))

      call occvtx4contr(1,occ_vtx,contr,op_info)

      count = 0

      ! loop over arcs and get sign changes from successive contractions
      do iarc0 = 1, narc
        iarc = arc_seq(iarc0)
        if (ntest.ge.100) then
          write(luout,*) 'iarc = ',iarc
          write(luout,*) 'present occ_vtx:'
          call wrt_occ_n(luout,occ_vtx,nvtx)
        end if

        count = mod(count,2)

        ivtx1 = arc(iarc)%link(1)
        ivtx2 = arc(iarc)%link(2)
        occ_cnt => arc(iarc)%occ_cnt

        if (ntest.ge.100) then
          write(luout,*) 'ivtx1/2: ',ivtx1,ivtx2
          write(luout,*) 'cnt='
          call wrt_occ(luout,occ_cnt)
        end if

        ! count enclosed indices
        nencl = 0
        do ivtx = ivtx1+1, ivtx2-1
          nencl = nencl + sum(occ_vtx(1:ngastp,1:2,ivtx))
        end do

        if (ntest.ge.100)
     &       write(luout,*) 'nencl = ',nencl

        ! remove CNT from vertices
        occ_vtx(1:,1:,ivtx1) =
     &       occ_vtx(1:,1:,ivtx1) - occ_cnt(1:,1:)
        occ_vtx(1:,1,ivtx2) =
     &       occ_vtx(1:,1,ivtx2) - occ_cnt(1:,2)
        occ_vtx(1:,2,ivtx2) =
     &       occ_vtx(1:,2,ivtx2) - occ_cnt(1:,1)

        if (ntest.ge.100) then
          write(luout,*) 'reduced occ_vtx:'
          call wrt_occ_n(luout,occ_vtx,nvtx)
        end if

        ! count remaining A of vtx1 and C of vtx2
        n_vtx1_arem = sum(occ_vtx(1:ngastp,2,ivtx1))
        n_vtx2_crem = sum(occ_vtx(1:ngastp,1,ivtx2))

        if (ntest.ge.100)
     &       write(luout,*) 'n_vtx1_arem,n_vtx2_crem: ',
     &                       n_vtx1_arem,n_vtx2_crem

        ! handle C part of contraction
        ! loop over HPVX/CA of contraction
        do hpvx = 1, ngastp
          iocc = occ_cnt(hpvx,1)
          if (iocc.eq.0) cycle
          nrem = n_vtx1_arem+n_vtx2_crem
          do hpvx2 = hpvx, ngastp     ! include remaining indices of same HPVX
            nrem = nrem + occ_vtx(hpvx2,1,ivtx1)
          end do
          do hpvx2 = hpvx, ngastp   ! include remaining indices of same HPVX
            nrem = nrem + occ_vtx(hpvx2,2,ivtx2)
          end do
          count = count + (nrem+nencl)*iocc
          if (ntest.eq.1000)
     &         write(luout,*) 'C,',hpvx,': iocc, nrem, nencl: ',
     &                                     iocc, nrem, nencl
        end do

        if (ntest.ge.100) write(luout,*) 'present count (C):   ',count

        ! handle A part of contraction
        do hpvx = 1, ngastp
          iocc = occ_cnt(hpvx,2)
          if (iocc.eq.0) cycle
          nrem = 0
          do hpvx2 = 1, hpvx-1
            nrem = nrem + occ_vtx(hpvx2,2,ivtx1)
          end do
          do hpvx2 = 1, hpvx-1
            nrem = nrem + occ_vtx(hpvx2,1,ivtx2)
          end do
          count = count + (nrem+nencl)*iocc
          if (ntest.eq.1000)
     &         write(luout,*) 'A,',hpvx,': iocc, nrem, nencl: ',
     &                                     iocc, nrem, nencl
        end do

        if (ntest.ge.100) write(luout,*) 'present count (C+A): ',count

      end do

      count = mod(count,2)

      if (nxarc.gt.0)
     &     call quit(1,'sign_global','adapt me for xarcs!')

      sign_global = 1
      if (count.eq.1) sign_global = -1

      deallocate(occ_vtx)

      if (ntest.ge.100) then
        write(luout,*) 'final sign: ',sign_global
      end if

      return
      end


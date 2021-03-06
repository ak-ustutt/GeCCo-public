*----------------------------------------------------------------------*
      subroutine get_reo_info(reo_op1op2,reo_other,
     &     iocc_op1op2,iocc_op1op2tmp,
     &     irst_op1op2,irst_op1op2tmp,
     &     njoined_op1op2,
     &     reo_info,str_info,orb_info)
*----------------------------------------------------------------------*
*     process raw reordering info from reduce_contr:
*
*      raise flag reo_op1op2, if op1op2 is reordered after contraction
*      raise flag reo_other,  if other operators need reordering
*
*      iocc_op1op2 is on entry the unreordered result
*      on exit, replace iocc_op1op2 by the actual result (after reo)
*      and set iocc_op1op2tmp to the temporary result
*      if no reordering of op1op2 occurs, both occupations are the same
*      (the same holds for irst_op1op2, irst_op1op2tmp)
*
*      
*
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'def_orbinf.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_reorder_info.h'
      include 'hpvxseq.h'

      integer, parameter ::
     &     ntest = 00

      type(reorder_info), intent(inout) ::
     &     reo_info
      type(orbinf), intent(in) ::
     &     orb_info
      type(strinf), intent(in) ::
     &     str_info
      logical, intent(out) ::
     &     reo_op1op2, reo_other
      integer, intent(in) ::
     &     njoined_op1op2
      integer, intent(inout) ::
     &     iocc_op1op2(ngastp,2,njoined_op1op2),
     &     iocc_op1op2tmp(ngastp,2,njoined_op1op2),
     &     irst_op1op2(2,orb_info%ngas,2,2,njoined_op1op2),
     &     irst_op1op2tmp(2,orb_info%ngas,2,2,njoined_op1op2)
      
      integer ::
     &     nreo, ireo, idx, nreo_op1op2, ireo_op1op2, nmap, len_map
      type(reorder_list), pointer ::
     &     reo(:)
      integer, pointer ::
     &     merge_map_stp1(:,:,:), merge_map_stp2(:,:,:),
     &     merge_stp1(:), merge_stp1inv(:),
     &     merge_stp2(:), merge_stp2inv(:),
     &     igrph(:,:,:), irst(:,:,:,:,:),
     &     from_to_vtx(:,:), is_op1op2(:)

      integer, external ::
     &     imltlist, sign_reo

c dbg
c      if (reo_info%nreo.gt.0) print *,'nreo = ',reo_info%nreo
c dbg
      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'get_reo_info')
        write(lulog,*) 'nreo = ',reo_info%nreo
        reo => reo_info%reo
        do ireo = 1, reo_info%nreo
          write(lulog,*) 'set # ',ireo
          write(lulog,*) 'is_bc_result: ',reo(ireo)%is_bc_result
          write(lulog,*) 'idxsuper:  ',reo(ireo)%idxsuper
          write(lulog,*) 'from/to:   ',reo(ireo)%from,reo(ireo)%to
          call wrt_occ(lulog,reo(ireo)%occ_shift)
        end do
      end if

      reo_op1op2 = .false.
      reo_other  = .false.

      ! copy input op1op2 occupation to op1op2tmp occ.
      iocc_op1op2tmp = iocc_op1op2
      irst_op1op2tmp = irst_op1op2

      ! fast exit, if nothing to do
      if (reo_info%nreo.eq.0) return

      nreo = reo_info%nreo
      reo => reo_info%reo

      ! scan over reordering info to find the number of operators
      ! to be reordered
      
      do ireo = 1, nreo
        reo_op1op2 = reo_op1op2.or.reo(ireo)%is_bc_result
        reo_other  = reo_other.or..not.reo(ireo)%is_bc_result
      end do

      ! current emergency exit
c      if (reo_other) then
c        call quit(1,'get_reo_info',
c     &       'reordering of operator other than OP1OP2 requested')
c      end if

      if (reo_op1op2) reo_info%n_op_reo = 1

      if (reo_op1op2) then
        ! make a merge-map for the primitive vertices in long form
        allocate(merge_map_stp1(njoined_op1op2,2,njoined_op1op2),
     &           merge_map_stp2(njoined_op1op2,2,njoined_op1op2))
        
        merge_map_stp1 = 0
        merge_map_stp2 = 0

        ! first node to merge for each target OP1OP2
        ! is the corresponding node in OP1OP2tmp
        do idx = 1, njoined_op1op2
          merge_map_stp1(1,1,idx) = idx
          merge_map_stp2(1,1,idx) = idx
        end do

        nreo_op1op2 = 0
        do ireo = 1, nreo
          if (.not.reo(ireo)%is_bc_result) cycle
          nreo_op1op2 = nreo_op1op2+1

          ! update op1op2, part I: remove shift occupation
          iocc_op1op2(1:ngastp,1:2,reo(ireo)%from) 
     &         = iocc_op1op2(1:ngastp,1:2,reo(ireo)%from)
     &         - reo(ireo)%occ_shift

          ! update op1op2, part II: add shift occupation
          iocc_op1op2(1:ngastp,1:2,reo(ireo)%to) 
     &         = iocc_op1op2(1:ngastp,1:2,reo(ireo)%to)
     &         + reo(ireo)%occ_shift

          ! update mapping info
          nmap = njoined_op1op2
     &         - imltlist(0,merge_map_stp1(1,2,reo(ireo)%from),
     &                    njoined_op1op2,1)+1
          merge_map_stp1(nmap,2,reo(ireo)%from) = nreo_op1op2
          nmap = njoined_op1op2
     &         - imltlist(0,merge_map_stp2(1,2,reo(ireo)%to),
     &                    njoined_op1op2,1)+1
          merge_map_stp2(nmap,2,reo(ireo)%to) = nreo_op1op2
        end do

        ! store shift occupations in iocc_reo array
        allocate(reo_info%iocc_reo(ngastp,2,nreo_op1op2))
        allocate(reo_info%from_to(2,nreo_op1op2))
        allocate(from_to_vtx(2,nreo_op1op2),
     &           is_op1op2(reo_info%nvtx_contr))
        ! ... and the remainder of the operator during reo:
        allocate(reo_info%iocc_opreo0(ngastp,2,njoined_op1op2))

        reo_info%iocc_opreo0 = iocc_op1op2tmp

        ireo_op1op2 = 0
        from_to_vtx = 0
        is_op1op2   = 0
        do ireo = 1, nreo
          if (.not.reo(ireo)%is_bc_result) cycle
          ireo_op1op2 = ireo_op1op2+1
          reo_info%iocc_reo(1:ngastp,1:2,ireo_op1op2) =
     &       reo(ireo)%occ_shift
          reo_info%from_to(1,ireo_op1op2) = reo(ireo)%from
          reo_info%from_to(2,ireo_op1op2) = reo(ireo)%to
          from_to_vtx(1,ireo_op1op2) = reo(ireo)%from_vtx
          from_to_vtx(2,ireo_op1op2) = reo(ireo)%to_vtx
          is_op1op2(reo(ireo)%from_vtx) = 1
          is_op1op2(reo(ireo)%to_vtx)   = 1
          reo_info%iocc_opreo0(1:ngastp,1:2,reo(ireo)%from) =
     &         reo_info%iocc_opreo0(1:ngastp,1:2,reo(ireo)%from) -
     &         reo(ireo)%occ_shift
        end do

        if (ntest.ge.100) then
          write(lulog,*) 'OP1OP2 tmp:'
          call wrt_occ_n(lulog,iocc_op1op2tmp,njoined_op1op2)
          write(lulog,*) 'OP1OP2:'
          call wrt_occ_n(lulog,iocc_op1op2,njoined_op1op2)
          write(lulog,*) 'REO:'
          call wrt_occ_n(lulog,reo_info%iocc_reo,nreo_op1op2)
          write(lulog,*) 'OPREO_0:'
          call wrt_occ_n(lulog,reo_info%iocc_opreo0,njoined_op1op2)
          write(lulog,*) 'is_op1op2:'
          write(lulog,'(1x,20i3)') is_op1op2(1:reo_info%nvtx_contr)
        end if

        reo_info%sign_reo = sign_reo(
     &       iocc_op1op2tmp,reo_info%iocc_opreo0,
     &       njoined_op1op2,reo_info%iocc_reo,
     &       reo_info%from_to,nreo_op1op2,0,
     &       from_to_vtx,reo_info%nca_vtx,is_op1op2,reo_info%nvtx_contr)

        call dummy_restr(irst_op1op2,
     &       iocc_op1op2,njoined_op1op2,orb_info)
        call dummy_restr(irst_op1op2tmp,
     &       iocc_op1op2tmp,njoined_op1op2,orb_info)

        ! transform merge-map to condensed representation
        ! length of map: 1 entry for each target vertex
        !  + number of non-zero entries in merge_map_xxxx
        len_map = njoined_op1op2*2 + 
     &       njoined_op1op2*njoined_op1op2*2-imltlist(0,merge_map_stp1,
     &                   njoined_op1op2*njoined_op1op2*2,1)
        allocate(merge_stp1(len_map),merge_stp1inv(len_map))
        len_map = njoined_op1op2*2 + 
     &       njoined_op1op2*njoined_op1op2*2-imltlist(0,merge_map_stp2,
     &       njoined_op1op2*njoined_op1op2*2,1)
        allocate(merge_stp2(len_map),merge_stp2inv(len_map))

        call condense_merge_map(merge_stp1,
     &       merge_map_stp1,njoined_op1op2,njoined_op1op2,.false.)
        call condense_merge_map(merge_stp1inv,
     &       merge_map_stp1,njoined_op1op2,njoined_op1op2,.true.)
        call condense_merge_map(merge_stp2,
     &       merge_map_stp2,njoined_op1op2,njoined_op1op2,.false.)
        call condense_merge_map(merge_stp2inv,
     &       merge_map_stp2,njoined_op1op2,njoined_op1op2,.true.)

        deallocate(merge_map_stp1,merge_map_stp2)

        ! transform to "condensed" info:

        call get_num_subblk(reo_info%ncblk_reo,reo_info%nablk_reo,
     &       reo_info%iocc_reo,nreo_op1op2)
        call get_num_subblk(reo_info%ncblk_reo0,reo_info%nablk_reo0,
     &       reo_info%iocc_opreo0,njoined_op1op2)

        allocate(reo_info%cinfo_reo_c(max(1,reo_info%ncblk_reo),3),
     &           reo_info%cinfo_reo_a(max(1,reo_info%nablk_reo),3),
     &           reo_info%cinfo_opreo0c(max(1,reo_info%ncblk_reo0),3),
     &           reo_info%cinfo_opreo0a(max(1,reo_info%nablk_reo0),3))
        len_map = max(1,
     &       max(reo_info%ncblk_reo+reo_info%ncblk_reo0,
     &           reo_info%nablk_reo+reo_info%nablk_reo0)*
     &           2*(njoined_op1op2+nreo_op1op2))
c dbg
c        print *,'??',njoined_op1op2,nreo_op1op2
c        print *,'len_map (2) = ',len_map,
c     &       reo_info%ncblk_reo,reo_info%nablk_reo
c dbg
        allocate(reo_info%map_reo1c(len_map),
     &           reo_info%map_reo1a(len_map),
     &           reo_info%map_reo2c(len_map),
     &           reo_info%map_reo2a(len_map))

        call set_mapping_info(reo_info%map_reo1c,reo_info%map_reo1a,
     &       1,
     &       reo_info%iocc_opreo0,njoined_op1op2,.false.,
     &       reo_info%iocc_reo,nreo_op1op2,.false.,
     &       iocc_op1op2tmp,merge_stp1,njoined_op1op2,hpvxblkseq)
c     &       iocc_op1op2,merge_stp1,njoined_op1op2,hpvxblkseq)
        call set_mapping_info(reo_info%map_reo1c,reo_info%map_reo1a,
     &       2,
     &       reo_info%iocc_reo,nreo_op1op2,.false.,
     &       reo_info%iocc_opreo0,njoined_op1op2,.false.,
     &       iocc_op1op2tmp,merge_stp1inv,njoined_op1op2,hpvxblkseq)
c     &       iocc_op1op2,merge_stp1inv,njoined_op1op2,hpvxblkseq)
        call set_mapping_info(reo_info%map_reo2c,reo_info%map_reo2a,
     &       1,
     &       reo_info%iocc_opreo0,njoined_op1op2,.false.,
     &       reo_info%iocc_reo,nreo_op1op2,.false.,
     &       iocc_op1op2,merge_stp2,njoined_op1op2,hpvxblkseq)
        call set_mapping_info(reo_info%map_reo2c,reo_info%map_reo2a,
     &       2,
     &       reo_info%iocc_reo,nreo_op1op2,.false.,
     &       reo_info%iocc_opreo0,njoined_op1op2,.false.,
     &       iocc_op1op2,merge_stp2inv,njoined_op1op2,hpvxblkseq)

        allocate(
     &       igrph(ngastp,2,max(njoined_op1op2,nreo_op1op2)),
     &       irst(2,orb_info%ngas,2,2,max(njoined_op1op2,nreo_op1op2)))

        call dummy_restr(irst,
     &       reo_info%iocc_opreo0,njoined_op1op2,orb_info)
        call get_grph4occ(igrph,
     &       reo_info%iocc_opreo0,irst,njoined_op1op2,
     &       str_info,orb_info,.true.)
        call condense_occ(reo_info%cinfo_opreo0c,reo_info%cinfo_opreo0a,
     &       reo_info%cinfo_opreo0c(1,3),reo_info%cinfo_opreo0a(1,3),
     &       reo_info%iocc_opreo0,njoined_op1op2,hpvxblkseq)
        call condense_occ(reo_info%cinfo_opreo0c(1,2),
     &                                      reo_info%cinfo_opreo0a(1,2),
     &       reo_info%cinfo_opreo0c(1,3),reo_info%cinfo_opreo0a(1,3),
     &       igrph,njoined_op1op2,hpvxblkseq)

        call dummy_restr(irst,
     &       reo_info%iocc_reo,nreo_op1op2,orb_info)
        call get_grph4occ(igrph,reo_info%iocc_reo,irst,nreo_op1op2,
     &       str_info,orb_info,.true.)
        call condense_occ(reo_info%cinfo_reo_c,reo_info%cinfo_reo_a,
     &       reo_info%cinfo_reo_c(1,3),reo_info%cinfo_reo_a(1,3),
     &       reo_info%iocc_reo,nreo_op1op2,hpvxblkseq)
        call condense_occ(reo_info%cinfo_reo_c(1,2),
     &                                       reo_info%cinfo_reo_a(1,2),
     &       reo_info%cinfo_reo_c(1,3),reo_info%cinfo_reo_a(1,3),
     &       igrph,nreo_op1op2,hpvxblkseq)

        deallocate(igrph,irst,merge_stp1,merge_stp1inv,
     &                        merge_stp2,merge_stp2inv,
     &             is_op1op2,from_to_vtx)

      end if

      return
      end

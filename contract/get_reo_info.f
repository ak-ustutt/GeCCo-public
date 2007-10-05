*----------------------------------------------------------------------*
      subroutine get_reo_info(reo_op1op2,reo_other,
     &     iocc_op1op2,iocc_op1op2tmp,
     &     irst_op1op2,irst_op1op2tmp,
     &     njoined_op1op2,
     &     contr,reo_info,str_info,orb_info)
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
     &     ntest = 100

      type(contraction), intent(in) ::
     &     contr
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
     &     nreo, ireo, idx, nreo_op1op2, nmap, len_map
      type(reorder_list), pointer ::
     &     reo(:)
      integer, pointer ::
     &     merge_map_stp1(:,:,:), merge_map_stp2(:,:,:),
     &     igrph(:,:,:), irst(:,:,:,:,:)

      integer, external ::
     &     imltlist

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'get_reo_info')
        write(luout,*) 'nreo = ',reo_info%nreo
        reo => reo_info%reo
        do ireo = 1, reo_info%nreo
          write(luout,*) 'set # ',ireo
          write(luout,*) 'is_bc_result: ',reo(ireo)%is_bc_result
          write(luout,*) 'idxsuper:  ',reo(ireo)%idxsuper
          write(luout,*) 'from/to:   ',reo(ireo)%from,reo(ireo)%to
          call wrt_occ(luout,reo(ireo)%occ_shift)
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
      if (reo_other) then
        call quit(1,'get_reo_info',
     &       'reordering of operator other than OP1OP2 requested')
      end if

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
          iocc_op1op2(1:ngastp,1:2,reo(ireo)%from) 
     &         = iocc_op1op2(1:ngastp,1:2,reo(ireo)%from)
     &         - reo(ireo)%occ_shift
          iocc_op1op2(1:ngastp,1:2,reo(ireo)%to) 
     &         = iocc_op1op2(1:ngastp,1:2,reo(ireo)%to)
     &         + reo(ireo)%occ_shift
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
        ! ... and the remainder of the operator during reo:
        allocate(reo_info%iocc_opreo0(ngastp,2,njoined_op1op2))

        reo_info%iocc_opreo0 = iocc_op1op2tmp

        do ireo = 1, nreo
          if (.not.reo(ireo)%is_bc_result) cycle
          reo_info%iocc_reo(1:ngastp,1:2,nreo_op1op2) =
     &       reo(ireo)%occ_shift
          reo_info%iocc_opreo0(1:ngastp,1:2,reo(ireo)%from) =
     &         reo_info%iocc_opreo0(1:ngastp,1:2,reo(ireo)%from) -
     &         reo(ireo)%occ_shift
        end do

        if (ntest.ge.100) then
          write(luout,*) 'OP1OP2 tmp:'
          call wrt_occ_n(luout,iocc_op1op2tmp,njoined_op1op2)
          write(luout,*) 'OP1OP2:'
          call wrt_occ_n(luout,iocc_op1op2,njoined_op1op2)
          write(luout,*) 'REO:'
          call wrt_occ_n(luout,reo_info%iocc_reo,nreo_op1op2)
          write(luout,*) 'OPREO_0:'
          call wrt_occ_n(luout,reo_info%iocc_opreo0,njoined_op1op2)
        end if

        call dummy_restr(irst_op1op2,
     &       iocc_op1op2,njoined_op1op2,orb_info%ihpvgas,orb_info%ngas)
        call dummy_restr(irst_op1op2tmp,
     &       iocc_op1op2tmp,njoined_op1op2,
     &                                  orb_info%ihpvgas,orb_info%ngas)

        ! transform merge-map to condensed representation
        ! length of map: 1 entry for each target vertex
        !  + number of non-zero entries in merge_map_xxxx
        len_map = njoined_op1op2 + 
     &       njoined_op1op2*njoined_op1op2*2-imltlist(0,merge_map_stp1,
     &                   njoined_op1op2*njoined_op1op2*2,1)
c dbg
        print *,'len_map = ',len_map
c dbg
        allocate(reo_info%merge_stp1(len_map))
        len_map = njoined_op1op2 + 
     &       njoined_op1op2*njoined_op1op2*2-imltlist(0,merge_map_stp2,
     &       njoined_op1op2*njoined_op1op2*2,1)
c dbg
        print *,'len_map = ',len_map
c dbg
        allocate(reo_info%merge_stp2(len_map))

        call condense_merge_map(reo_info%merge_stp1,
     &       merge_map_stp1,njoined_op1op2,njoined_op1op2,.false.)
        call condense_merge_map(reo_info%merge_stp2,
     &       merge_map_stp2,njoined_op1op2,njoined_op1op2,.false.)

        deallocate(merge_map_stp1,merge_map_stp2)

        ! transform to "condensed" info:

        call get_num_subblk(reo_info%ncblk_reo,reo_info%nablk_reo,
     &       reo_info%iocc_reo,nreo_op1op2)
        call get_num_subblk(reo_info%ncblk_reo0,reo_info%nablk_reo0,
     &       reo_info%iocc_opreo0,njoined_op1op2)

        allocate(reo_info%cinfo_reo_c(reo_info%ncblk_reo,3),
     &           reo_info%cinfo_reo_a(reo_info%nablk_reo,3),
     &           reo_info%cinfo_opreo0c(reo_info%ncblk_reo0,3),
     &           reo_info%cinfo_opreo0a(reo_info%nablk_reo0,3))
        len_map = 2*njoined_op1op2*nreo_op1op2
        allocate(reo_info%map_reo1c(max(1,reo_info%ncblk_reo*len_map)),
     &           reo_info%map_reo1a(max(1,reo_info%nablk_reo*len_map)),
     &           reo_info%map_reo2c(max(1,reo_info%ncblk_reo*len_map)),
     &           reo_info%map_reo2a(max(1,reo_info%nablk_reo*len_map)))

        call set_mapping_info(reo_info%map_reo1c,reo_info%map_reo1a,
     &       1,
     &       reo_info%iocc_opreo0,njoined_op1op2,.false.,
     &       reo_info%iocc_reo,nreo_op1op2,.false.,
     &       iocc_op1op2,reo_info%merge_stp1,njoined_op1op2,hpvxblkseq)
        call set_mapping_info(reo_info%map_reo1c,reo_info%map_reo1a,
     &       2,
     &       reo_info%iocc_reo,nreo_op1op2,.false.,
     &       reo_info%iocc_opreo0,njoined_op1op2,.false.,
     &       iocc_op1op2,reo_info%merge_stp1,njoined_op1op2,hpvxblkseq)
        call set_mapping_info(reo_info%map_reo2c,reo_info%map_reo2a,
     &       1,
     &       reo_info%iocc_opreo0,njoined_op1op2,.false.,
     &       reo_info%iocc_reo,nreo_op1op2,.false.,
     &       iocc_op1op2,reo_info%merge_stp2,njoined_op1op2,hpvxblkseq)
        call set_mapping_info(reo_info%map_reo2c,reo_info%map_reo2a,
     &       2,
     &       reo_info%iocc_reo,nreo_op1op2,.false.,
     &       reo_info%iocc_opreo0,njoined_op1op2,.false.,
     &       iocc_op1op2,reo_info%merge_stp2,njoined_op1op2,hpvxblkseq)

        allocate(
     &       igrph(ngastp,2,max(njoined_op1op2,nreo_op1op2)),
     &       irst(2,orb_info%ngas,2,2,max(njoined_op1op2,nreo_op1op2)))

        call dummy_restr(irst,
     &       reo_info%iocc_opreo0,njoined_op1op2,
     &       orb_info%ihpvgas,orb_info%ngas)
        call get_grph4occ(igrph,
     &       reo_info%iocc_opreo0,irst,
     &       str_info,orb_info%ihpvgas,
     &       orb_info%ngas,njoined_op1op2,.true.)
        call condense_occ(reo_info%cinfo_opreo0c,reo_info%cinfo_opreo0a,
     &       reo_info%cinfo_opreo0c(1,3),reo_info%cinfo_opreo0a(1,3),
     &       reo_info%iocc_opreo0,njoined_op1op2,hpvxblkseq)
        call condense_occ(reo_info%cinfo_opreo0c(1,2),
     &                                      reo_info%cinfo_opreo0a(1,2),
     &       reo_info%cinfo_opreo0c(1,3),reo_info%cinfo_opreo0a(1,3),
     &       igrph,njoined_op1op2,hpvxblkseq)

        call dummy_restr(irst,
     &       reo_info%iocc_reo,nreo_op1op2,
     &       orb_info%ihpvgas,orb_info%ngas)
        call get_grph4occ(igrph,reo_info%iocc_reo,irst,
     &       str_info,orb_info%ihpvgas,
     &       orb_info%ngas,nreo_op1op2,.true.)
        call condense_occ(reo_info%cinfo_reo_c,reo_info%cinfo_reo_a,
     &       reo_info%cinfo_reo_c(1,3),reo_info%cinfo_reo_a(1,3),
     &       reo_info%iocc_reo,nreo_op1op2,hpvxblkseq)
        call condense_occ(reo_info%cinfo_reo_c(1,2),
     &                                       reo_info%cinfo_reo_a(1,2),
     &       reo_info%cinfo_reo_c(1,3),reo_info%cinfo_reo_a(1,3),
     &       reo_info%iocc_reo,nreo_op1op2,hpvxblkseq)

        deallocate(igrph,irst)

      end if

      return
      end

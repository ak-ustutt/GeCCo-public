*----------------------------------------------------------------------*
      subroutine get_reo_info2(
     &     mode,idxsuper,
     &     iocc_op1op2,iocc_op1op2tmp,
     &     irst_op1op2,irst_op1op2tmp,
     &     njoined_op1op2,mst,gamt,
     &     merge_stp1,merge_stp1inv,merge_stp2,merge_stp2inv,
     &     occ_vtx,rst_vtx,svertex,info_vtx,nj_res,nvtx,
     &     reo_info,str_info,orb_info)
*----------------------------------------------------------------------*
*     process raw reordering info from reduce_contr:
*
*     mode = 1: pre-processing, when generating opt. formula
*           -1:  dto. with some other assumption on occ_vtx
*     mode = 2: post-processing, when evaluating formula
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
c      logical, intent(out) ::
c     &     reo_op1op2, reo_other
      integer, intent(in) ::
     &     mode, nj_res,idxsuper,nvtx,
     &     occ_vtx(ngastp,2,*),rst_vtx(2,orb_info%ngas,2,2,*),
     &     svertex(*),info_vtx(2,*)
      integer, intent(inout) ::
     &     njoined_op1op2,gamt,mst,
     &     merge_stp1(*), merge_stp1inv(*),
     &     merge_stp2(*), merge_stp2inv(*)
      integer, intent(inout), target ::
     &     iocc_op1op2(ngastp,2,*),
     &     iocc_op1op2tmp(ngastp,2,*),
     &     irst_op1op2(2,orb_info%ngas,2,2,*),
     &     irst_op1op2tmp(2,orb_info%ngas,2,2,*)
      
      integer ::
     &     nreo, ireo, idx, nreo_op1op2, ireo_op1op2, nmap, len_map,
     &     ivtx, jvtx, idxsuper_old1, idxsuper_old2, ivtx1, ivtx2,
     &     ngas, ij, streo, ndreo
      integer ::
     &     rst_shift(2,orb_info%ngas,2,2)
      type(reorder_list), pointer ::
     &     reo(:)
      integer, pointer ::
     &     merge_map_stp1(:,:,:), merge_map_stp2(:,:,:),
     &     igrph(:,:,:), irst0(:,:,:,:,:), irst_reo(:,:,:,:,:),
     &     from_to_vtx(:,:), is_op1op2(:),
     &     occ_shift(:,:), 
     &     occ_target(:,:,:), rst_target(:,:,:,:,:)

      integer, external ::
     &     imltlist, sign_reo


      if (abs(mode).ne.1.and.mode.ne.2)
     &     call quit(1,'get_reo_info2','invalid mode')
      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'get_reo_info2')
        write(luout,*) 'nreo = ',reo_info%nreo
        write(luout,*) 'mode:  ',mode
        if (abs(mode).eq.1) then
          reo => reo_info%reo
          do ireo = 1, reo_info%nreo
            write(luout,*) 'set # ',ireo
            write(luout,*) 'is_bc_result: ',reo(ireo)%is_bc_result
            write(luout,*) 'idxsuper:  ',reo(ireo)%idxsuper
            write(luout,*) 'from/to:   ',reo(ireo)%from,reo(ireo)%to
            call wrt_occ(luout,reo(ireo)%occ_shift)
          end do
        end if
      end if

      if (reo_info%nreo.eq.0) call quit(1,'get_reo_info2','nothing?')

      ngas = orb_info%ngas
      nreo = reo_info%nreo
      reo => reo_info%reo

      ! scan over reordering info to find the number of operators
      ! to be reordered
      
c      do ireo = 1, nreo
c        reo_op1op2 = reo_op1op2.or.reo(ireo)%is_bc_result
c        reo_other  = reo_other.or..not.reo(ireo)%is_bc_result
c      end do

      ! current emergency exit
c      if (reo_other) then
c        call quit(1,'get_reo_info',
c     &       'reordering of operator other than OP1OP2 requested')
c      end if

      if (abs(mode).eq.1) then
      jvtx = 0
      do ivtx = 1, nvtx
        do ireo = 1, nreo
c dbg
c          print *,'ivtx,ireo,idxsuper(ist/soll): ',ivtx,ireo,
c     &         reo(ireo)%idxsuper,idxsuper
c dbg
          if (reo(ireo)%idxsuper.ne.idxsuper) cycle
          idxsuper_old1 = svertex(reo(ireo)%from_vtx)
          idxsuper_old2 = svertex(reo(ireo)%to_vtx)
c dbg
c          print *,'svertex(ivtx),idxsuper_old1,idxsuper_old2:',
c     &             svertex(ivtx),idxsuper_old1,idxsuper_old2
c dbg
          if (svertex(ivtx).ne.idxsuper_old1.and.
     &        svertex(ivtx).ne.idxsuper_old2) cycle
          jvtx = jvtx+1
          iocc_op1op2(1:ngastp,1:2,jvtx) =
     &         occ_vtx(1:ngastp,1:2,nj_res+ivtx)
          irst_op1op2(1:2,1:ngas,1:2,1:2,jvtx) =
     &         rst_vtx(1:2,1:ngas,1:2,1:2,nj_res+ivtx)
c dbg
c          print *,'stored: ',ivtx,' as ',jvtx
c dbg
          mst  = info_vtx(1,nj_res+ivtx)
          gamt = info_vtx(2,nj_res+ivtx)
          exit ! ireo loop: consider ivtx only once
        end do
      end do
      njoined_op1op2 = jvtx

      iocc_op1op2tmp(1:ngastp,1:2,1:njoined_op1op2)
     &     = iocc_op1op2(1:ngastp,1:2,1:njoined_op1op2)
      irst_op1op2tmp(1:2,1:ngas,1:2,1:2,1:njoined_op1op2)
     &     = irst_op1op2(1:2,1:ngas,1:2,1:2,1:njoined_op1op2)


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

        ! mode=1: loop backwards over reorderings (to avoid negative occ.s)
        nreo_op1op2 = 0
        do ireo = 1, nreo
          if (reo(ireo)%idxsuper.ne.idxsuper) cycle
          nreo_op1op2 = nreo_op1op2+1

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
        if (mode.eq.1) then
          streo = nreo
          ndreo = 1
        else
          streo = 1
          ndreo = nreo
        end if
        do ireo = streo, ndreo, -mode !mode=1: backwards loop
          if (reo(ireo)%idxsuper.ne.idxsuper) cycle

          ! shift direction and target depends on mode
          if (mode.eq.-1) then
            ivtx1 = reo(ireo)%from
            ivtx2 = reo(ireo)%to
            occ_target!(1:ngastp,1:2,1:njoined_op1op2)
     &           => iocc_op1op2(1:ngastp,1:2,1:njoined_op1op2)
            rst_target!(1:2,1:ngas,1:2,1:2,1:njoined_op1op2)
     &           =>
     &      irst_op1op2(1:2,1:ngas,1:2,1:2,1:njoined_op1op2)
          else
            ivtx1 = reo(ireo)%to
            ivtx2 = reo(ireo)%from
            occ_target!(1:ngastp,1:2,1:njoined_op1op2)
     &           => iocc_op1op2tmp(1:ngastp,1:2,1:njoined_op1op2)
            rst_target!(1:2,1:ngas,1:2,1:2,1:njoined_op1op2)
     &           =>
     &      irst_op1op2tmp(1:2,1:ngas,1:2,1:2,1:njoined_op1op2)
          end if
          occ_shift => reo(ireo)%occ_shift
          
          ! handle restrictions first ...
          call fit_restr(rst_shift,occ_shift,
     &         rst_target(1,1,1,1,ivtx1),orb_info%ihpvgas,
     &                                   ngas)

          rst_target(1:2,1:ngas,1:2,1:2,ivtx1)
     &         = rst_target(1:2,1:ngas,1:2,1:2,ivtx1)
     &         - rst_shift
          rst_target(1:2,1:ngas,1:2,1:2,ivtx2)
     &         = rst_target(1:2,1:ngas,1:2,1:2,ivtx2)
     &         + rst_shift

          ! ... and adapt now occupations
          occ_target(1:ngastp,1:2,ivtx1) 
     &         = occ_target(1:ngastp,1:2,ivtx1)
     &         - occ_shift
          occ_target(1:ngastp,1:2,ivtx2) 
     &         = occ_target(1:ngastp,1:2,ivtx2)
     &         + occ_shift
          
          
c          if (mode.eq.-1) then
c            ! update op1op2, part I: remove shift occupation
c            iocc_op1op2(1:ngastp,1:2,reo(ireo)%from) 
c     &         = iocc_op1op2(1:ngastp,1:2,reo(ireo)%from)
c     &         - reo(ireo)%occ_shift
c
c            ! update op1op2, part II: add shift occupation
c            iocc_op1op2(1:ngastp,1:2,reo(ireo)%to) 
c     &         = iocc_op1op2(1:ngastp,1:2,reo(ireo)%to)
c     &         + reo(ireo)%occ_shift
c          else
c            ! update op1op2, part I: remove shift occupation
c            iocc_op1op2tmp(1:ngastp,1:2,reo(ireo)%from) 
c     &         = iocc_op1op2tmp(1:ngastp,1:2,reo(ireo)%from)
c     &         + reo(ireo)%occ_shift
c
c            ! update op1op2, part II: add shift occupation
c            iocc_op1op2tmp(1:ngastp,1:2,reo(ireo)%to) 
c     &         = iocc_op1op2tmp(1:ngastp,1:2,reo(ireo)%to)
c     &         - reo(ireo)%occ_shift
c          end if
        end do

        ! store shift occupations in iocc_reo array
        allocate(reo_info%iocc_reo(ngastp,2,nreo_op1op2))
        allocate(reo_info%from_to(2,nreo_op1op2))
        allocate(from_to_vtx(2,nreo_op1op2),
     &           is_op1op2(max(reo_info%nvtx_contr,njoined_op1op2)))
        ! ... and the remainder of the operator during reo:
        allocate(reo_info%iocc_opreo0(ngastp,2,njoined_op1op2))

        reo_info%iocc_opreo0(1:ngastp,1:2,1:njoined_op1op2)
     &       = iocc_op1op2tmp(1:ngastp,1:2,1:njoined_op1op2)

        ireo_op1op2 = 0
        from_to_vtx = 0
        is_op1op2   = 0
        do ireo = 1, nreo
          if (reo(ireo)%idxsuper.ne.idxsuper) cycle
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
          write(luout,*) 'OP1OP2 tmp:'
          call wrt_occ_n(luout,iocc_op1op2tmp,njoined_op1op2)
          write(luout,*) 'OP1OP2:'
          call wrt_occ_n(luout,iocc_op1op2,njoined_op1op2)
          write(luout,*) 'REO:'
          call wrt_occ_n(luout,reo_info%iocc_reo,nreo_op1op2)
          write(luout,*) 'OPREO_0:'
          call wrt_occ_n(luout,reo_info%iocc_opreo0,njoined_op1op2)
          write(luout,*) 'is_op1op2:'
          write(luout,'(1x,20i3)') is_op1op2(1:reo_info%nvtx_contr)
        end if
c patch
c        print *,'redefining is_op1op2:'
        is_op1op2 = 0
        do ivtx = 1, max(reo_info%nvtx_contr,njoined_op1op2)
          if (svertex(ivtx).eq.idxsuper) is_op1op2(ivtx) = 1
        end do
c        write(luout,'(1x,20i3)') is_op1op2(1:reo_info%nvtx_contr)
c        write(luout,'(1x,20i3)') reo_info%nca_vtx(1:reo_info%nvtx_contr)
c        write(luout,*) ' idxsuper = ',idxsuper
c        write(luout,'(1x,20i3)') svertex(1:reo_info%nvtx_contr)
c patch

        reo_info%sign_reo = sign_reo(
     &       iocc_op1op2tmp,reo_info%iocc_opreo0,
     &       njoined_op1op2,reo_info%iocc_reo,
     &       reo_info%from_to,nreo_op1op2,
     &       from_to_vtx,reo_info%nca_vtx,is_op1op2,
     &       max(reo_info%nvtx_contr,njoined_op1op2))

c        ! to be commented out ...
c        call dummy_restr(irst_op1op2,
c     &       iocc_op1op2,njoined_op1op2,orb_info)
c        call dummy_restr(irst_op1op2tmp,
c     &       iocc_op1op2tmp,njoined_op1op2,orb_info)

        ! transform merge-map to condensed representation
        ! length of map: 1 entry for each target vertex
        !  + number of non-zero entries in merge_map_xxxx
        call condense_merge_map(merge_stp1,
     &       merge_map_stp1,njoined_op1op2,njoined_op1op2,.false.)
        call condense_merge_map(merge_stp1inv,
     &       merge_map_stp1,njoined_op1op2,njoined_op1op2,.true.)
        call condense_merge_map(merge_stp2,
     &       merge_map_stp2,njoined_op1op2,njoined_op1op2,.false.)
        call condense_merge_map(merge_stp2inv,
     &       merge_map_stp2,njoined_op1op2,njoined_op1op2,.true.)

        deallocate(merge_map_stp1,merge_map_stp2)

        deallocate(is_op1op2,from_to_vtx)

      end if

      if (mode.eq.2) then

        nreo_op1op2 = reo_info%nreo
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
        allocate(reo_info%map_reo1c(len_map),
     &           reo_info%map_reo1a(len_map),
     &           reo_info%map_reo2c(len_map),
     &           reo_info%map_reo2a(len_map))

        call set_mapping_info(reo_info%map_reo1c,reo_info%map_reo1a,
     &       1,
     &       reo_info%iocc_opreo0,njoined_op1op2,.false.,
     &       reo_info%iocc_reo,nreo_op1op2,.false.,
     &       iocc_op1op2tmp,merge_stp1,njoined_op1op2,hpvxblkseq)
        call set_mapping_info(reo_info%map_reo1c,reo_info%map_reo1a,
     &       2,
     &       reo_info%iocc_reo,nreo_op1op2,.false.,
     &       reo_info%iocc_opreo0,njoined_op1op2,.false.,
     &       iocc_op1op2tmp,merge_stp1inv,njoined_op1op2,hpvxblkseq)
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
     &       irst0(2,ngas,2,2,max(njoined_op1op2,nreo_op1op2)),
     &       irst_reo(2,ngas,2,2,max(njoined_op1op2,nreo_op1op2)))

c        call dummy_restr(irst,
c     &       reo_info%iocc_opreo0,njoined_op1op2,orb_info)

c        call fit_restr2(irst,
c     &       reo_info%iocc_opreo0,njoined_op1op2,
c     &       irst_op1op2,njoined_op1op2,orb_info%ihpvgas,ngas)

c dbg
c        print *,'new call to fit_restr3'
c dbg
        call fit_restr3(irst0,irst_reo,
     &       iocc_op1op2tmp, reo_info%iocc_opreo0,reo_info%iocc_reo,
     &       irst_op1op2tmp, merge_stp1inv,
     &       njoined_op1op2,nreo_op1op2,
     &       orb_info%ihpvgas,ngas,orb_info%nspin)
c dbg
c        print *,'OP1OP2TMP'
c         do ij = 1, njoined_op1op2
c            call wrt_occ_rstr(luout,ij,iocc_op1op2tmp(:,:,ij),
c     &                                 irst_op1op2tmp(:,:,:,:,ij),
c     &                                 ngas,orb_info%nspin)
c          end do
c        print *,'OP0'
c         do ij = 1, njoined_op1op2
c            call wrt_occ_rstr(luout,ij,reo_info%iocc_opreo0(:,:,ij),
c     &                                    irst0(:,:,:,:,ij),
c     &                                    ngas,orb_info%nspin)
c          end do
c        print *,'REO'
c         do ij = 1, nreo_op1op2
c            call wrt_occ_rstr(luout,ij,reo_info%iocc_reo(:,:,ij),
c     &                                    irst_reo(:,:,:,:,ij),
c     &                                    ngas,orb_info%nspin)
c          end do
c dbg

        call get_grph4occ(igrph,
     &       reo_info%iocc_opreo0,irst0,njoined_op1op2,
     &       str_info,orb_info,.true.)
        call condense_occ(reo_info%cinfo_opreo0c,reo_info%cinfo_opreo0a,
     &       reo_info%cinfo_opreo0c(1,3),reo_info%cinfo_opreo0a(1,3),
     &       reo_info%iocc_opreo0,njoined_op1op2,hpvxblkseq)
        call condense_occ(reo_info%cinfo_opreo0c(1,2),
     &                                      reo_info%cinfo_opreo0a(1,2),
     &       reo_info%cinfo_opreo0c(1,3),reo_info%cinfo_opreo0a(1,3),
     &       igrph,njoined_op1op2,hpvxblkseq)

c        call dummy_restr(irst,
c     &       reo_info%iocc_reo,nreo_op1op2,orb_info)

        call get_grph4occ(igrph,reo_info%iocc_reo,irst_reo,nreo_op1op2,
     &       str_info,orb_info,.true.)
        call condense_occ(reo_info%cinfo_reo_c,reo_info%cinfo_reo_a,
     &       reo_info%cinfo_reo_c(1,3),reo_info%cinfo_reo_a(1,3),
     &       reo_info%iocc_reo,nreo_op1op2,hpvxblkseq)
        call condense_occ(reo_info%cinfo_reo_c(1,2),
     &                                       reo_info%cinfo_reo_a(1,2),
     &       reo_info%cinfo_reo_c(1,3),reo_info%cinfo_reo_a(1,3),
     &       igrph,nreo_op1op2,hpvxblkseq)

        deallocate(igrph,irst0,irst_reo)

      end if

      return
      end

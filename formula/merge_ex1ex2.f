*----------------------------------------------------------------------*
      subroutine merge_ex1ex2(iocc_op1op2,njoined_op1op2,
     &                      merge_map,
     &                      ld_map,
     &                      ivtxsuper1, ivtxsuper2, last_cntr,
     &                      iocc_ex1ex2,inum_ori,njoined12,
     &                      arc_list,len_list,
     &                      contr,occ_vtx,njoined_res)
*----------------------------------------------------------------------*
*     merge the raw external occupations giving the occupation of the
*     final result of the binary contraction
*     last_cntr: this is the last contraction so make sure that 
*        the occupation fits to njoined_res
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'ifc_operators.h'

      integer, parameter ::
     &     ntest = 00

      logical, intent(in) ::
     &     last_cntr
      integer, intent(in) ::
     &     njoined12, njoined_res, len_list, ld_map,
     &     ivtxsuper1,ivtxsuper2,
     &     iocc_ex1ex2(ngastp,2,njoined12), inum_ori(2,njoined12),
     &     arc_list(len_list),
     &     occ_vtx(ngastp,2,*)
      type(contraction), intent(in) ::
     &     contr
      integer, intent(out) ::
     &     iocc_op1op2(ngastp,2,njoined12), njoined_op1op2,
     &     merge_map(ld_map,2,njoined12)

      integer ::
     &     ijoin, nvtx, ivtx1, ivtx2, ilist, iarc, idx_merge,
     &     nmvleft, idx1, idx2, nmap1, nmap2, ii, ivtx, idum, mdx
      integer ::
     &     svmap(contr%nvtx), topomap(contr%nvtx,contr%nvtx),
     &     imvleft(contr%nvtx)
      integer, external ::
     &     idx_merge_vtx1vtx2, imltlist
      logical, external ::
     &     merge_check

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'merge_ex1ex2')
        write(luout,*) 'EX1EX2 on entry:'
        call wrt_occ_n(luout,iocc_ex1ex2,njoined12)
c dbg
        call prt_contr3(6,contr,occ_vtx(1,1,1+njoined_res))
c dbg
      end if

      ! set up some help arrays for merge-or-not routine
      if (njoined_res.eq.1) then
        svmap(1:contr%nvtx) = 1
      else
        call svmap4contr(svmap,contr,occ_vtx,njoined_res)
      end if
c dbg
c      print *,'svmap : ',svmap(1:contr%nvtx)
c      print *,'njoined_res: ',njoined_res
c      print *,'njoined12:   ',njoined12
c dbg

      ! easy game then: svmap contains the info needed:
      if (last_cntr) then
        merge_map = 0
        iocc_op1op2(1:ngastp,1:2,1:njoined_res) = 0
        do ivtx1 = 1, contr%nvtx
          idx1 = svmap(ivtx1)
          if (idx1.le.0) cycle
          idx2 = 1
          if (contr%svertex(ivtx1).eq.ivtxsuper2) idx2 = 2
          iocc_op1op2(1:ngastp,1:2,idx1) =
     &         iocc_op1op2(1:ngastp,1:2,idx1) +
     &         iocc_ex1ex2(1:ngastp,1:2,ivtx1)
          ! update map
          mdx = 1
          do while(merge_map(mdx,idx2,idx1).gt.0.and.mdx.le.ld_map)
            mdx = mdx+1
          end do
          if (mdx.gt.ld_map)
     &         call quit(1,'merge_ex1ex2','ld_map too small?')
          merge_map(mdx,idx2,idx1) = inum_ori(2,ivtx1)
        end do

        do ivtx = 1, njoined_res
          ! make sure that the map entries are sorted
          do ii = 1, 2
            nmap1 = ld_map - imltlist(0,merge_map(1,ii,ivtx),ld_map,1)
            if (nmap1.gt.1) call isort(merge_map(1,ii,ivtx),nmap1,+1)
          end do
        end do

        njoined_op1op2 = njoined_res

      else
      
        call topomap4contr(1,-1,-1,topomap,idum,idum,idum,contr,
     &                   occ_vtx(1,1,njoined_res+1))
      
        ! reset map
        merge_map = 0
        ! initial setting (before merge): each vertex contributes to
        ! the vertices as given in inum_ori
        do ijoin = 1, njoined12
          merge_map(1,inum_ori(1,ijoin),ijoin) = inum_ori(2,ijoin)
        end do

        iocc_op1op2 = iocc_ex1ex2

c      idxop1op2 = 0
        nvtx = contr%nvtx
        ! loop over vertex pairs
        ! (reversely, in order to process multi-merges correctly)
        do ivtx1 = nvtx-1, 1, -1
          ! interesting for us?
          if (contr%svertex(ivtx1).ne.ivtxsuper1.and.
     &        contr%svertex(ivtx1).ne.ivtxsuper2) cycle
          do ivtx2 = ivtx1+1, nvtx
            ! interesting for us?
            if (contr%svertex(ivtx2).ne.ivtxsuper1.and.
     &          contr%svertex(ivtx2).ne.ivtxsuper2) cycle

c          ! is there a contraction ?          
c          do ilist = 1, len_list
c            iarc = arc_list(ilist)
c            if (contr%arc(iarc)%link(1).eq.ivtx1.and.
c     &          contr%arc(iarc)%link(2).eq.ivtx2) exit
c            iarc = -1
c          end do
c          if (iarc.le.0) cycle

            ! idx1 on ex1ex2, as corresponding to ivtx1 
            idx1 = imltlist(ivtxsuper1,contr%svertex,ivtx1,1)
     &           + imltlist(ivtxsuper2,contr%svertex,ivtx1,1)
            ! idx2 on ex1ex2, as corresponding to ivtx2
            idx2 = imltlist(ivtxsuper1,contr%svertex,ivtx2,1)
     &           + imltlist(ivtxsuper2,contr%svertex,ivtx2,1)

            ! check whether we may merge the vertices
            idx_merge = idx_merge_vtx1vtx2(
     &           ivtx1,ivtx2,ivtxsuper1,ivtxsuper2,
     &           nmvleft,imvleft,contr%svertex,svmap,topomap,contr%nvtx)
          
            if (idx_merge.gt.0) then
              ! check whether a full merge introduces an unwanted
              ! symmetrization
              if (.not.merge_check(contr,ivtxsuper1,ivtxsuper2,
     &             ivtx1,ivtx2,iocc_op1op2(1,1,idx1),
     &                         iocc_op1op2(1,1,idx2),topomap)) then
                idx_merge = 0
              end if
            end if
          
            ! check for partial symmetrization (additional reforming
            ! step after contraction)

c          if (idx_merge.le.0) cycle

            if (ntest.ge.150)
     &           write(luout,*) 'vertices, indices, merge:',
     &           ivtx1,ivtx2,idx1,idx2,idx_merge.gt.0

            if (idx_merge.gt.0) then
              ! update occupation
              iocc_op1op2(1:ngastp,1:2,idx1) =
     &             iocc_op1op2(1:ngastp,1:2,idx1)
     &            +iocc_op1op2(1:ngastp,1:2,idx2)
              iocc_op1op2(1:ngastp,1:2,idx2) = 0
              ! update map
              do ii = 1, 2
                nmap1 = ld_map-imltlist(0,merge_map(1,ii,idx1),ld_map,1)
                nmap2 = ld_map-imltlist(0,merge_map(1,ii,idx2),ld_map,1)
                if (nmap1+nmap2.gt.ld_map)
     &               call quit(1,'merge_ex1ex2','ld_map too small?')
                merge_map(nmap1+1:nmap1+nmap2,ii,idx1) =
     &               merge_map(1:nmap2,ii,idx2)
                merge_map(1:nmap2,ii,idx2) = 0

              end do

              if (ntest.ge.150) then
                write(luout,*) 'updated OP1OP2:'
                call wrt_occ_n(luout,iocc_op1op2,njoined12)
              end if

            end if          

          end do
        end do

        if (ntest.ge.150) then
          write(luout,*) 'raw OP1OP2:'
          call wrt_occ_n(luout,iocc_op1op2,njoined12)
        end if

        njoined_op1op2 = 0
        do ivtx = 1, njoined12
          ! make sure that the map entries are sorted
          do ii = 1, 2
            nmap1 = ld_map - imltlist(0,merge_map(1,ii,ivtx),ld_map,1)
            if (nmap1.gt.1) call isort(merge_map(1,ii,ivtx),nmap1,+1)
          end do
          if (iocc_zero(iocc_op1op2(1:ngastp,1:2,ivtx))) cycle
          njoined_op1op2 = njoined_op1op2+1
          if (njoined_op1op2.lt.ivtx) then
            iocc_op1op2(1:ngastp,1:2,njoined_op1op2) =
     &           iocc_op1op2(1:ngastp,1:2,ivtx)
            merge_map(1:ld_map,1:2,njoined_op1op2) =
     &           merge_map(1:ld_map,1:2,ivtx)
          end if
        end do
        ! well, we keep at least one:
        njoined_op1op2 = max(1,njoined_op1op2)

      end if

      if (ntest.ge.100) then
        write(luout,*) 'final OP1OP2:'
        call wrt_occ_n(luout,iocc_op1op2,njoined_op1op2)
        write(luout,*) 'merge map:'
        do ivtx = 1, njoined_op1op2
          do ii = 1, 2
            write(luout,'(4x,i4," <-",i4," :",10i4)')
     &           ivtx,ii,merge_map(1:ld_map,ii,ivtx)
          end do
        end do
      end if


      end

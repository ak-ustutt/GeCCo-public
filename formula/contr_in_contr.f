*----------------------------------------------------------------------*
      logical function contr_in_contr(contra,contrb,op_info)
*----------------------------------------------------------------------*
*     check whether contraction A on contra is contained in contraction
*     B on contrb
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'

      integer, parameter ::
     &     ntest = 00

      type(contraction), intent(in) ::
     &     contra, contrb
      type(operator_info), intent(in) ::
     &     op_info

      integer ::
     &     nvtx_a, nvtx_b, narc_a, narc_b, base1, base2, idum

      integer, pointer ::
     &     occ_vtx_a(:,:,:), occ_vtx_b(:,:,:),
     &     topomap_a(:,:), topomap_b(:,:),
     &     vtxinf_a(:), vtxinf_b(:), vtxmap(:)
c dbg
c      integer ::
c     &     ieqvfac
c      integer, pointer ::
c     &     ivtx_reo(:), occ_vtx(:,:,:)
c      logical, pointer ::
c     &     fix_vtx(:)
c      logical ::
c     &     reo
c dbg
      integer, external ::
     &     maxblk_in_contr, ifndmax

      if (ntest.ge.100) then
        write(luout,*) '--------------------------'
        write(luout,*) 'output from contr_in_contr'
        write(luout,*) '--------------------------'
        write(luout,*) 'vertices A/B: ',contra%nvtx,contrb%nvtx
        write(luout,*) 'arcs A/B:     ',contra%narc,contrb%narc
c dbg
        call prt_contr2(6,contra,op_info)
        call prt_contr2(6,contrb,op_info)
c dbg
      end if

      contr_in_contr = .false.

      nvtx_a = contra%nvtx
      nvtx_b = contrb%nvtx
      narc_a = contra%narc
      narc_b = contrb%narc
      
      ! number of vertices and arcs of B must be >= than those of A
      ! else we can directly say good-bye
      if (nvtx_a.gt.nvtx_b.or.
     &    narc_a.gt.narc_b) return

      allocate(occ_vtx_a(ngastp,2,nvtx_a),
     &         occ_vtx_b(ngastp,2,nvtx_b),
     &         topomap_a(nvtx_a,nvtx_a), topomap_b(nvtx_b,nvtx_b),
     &         vtxinf_a(nvtx_a), vtxinf_b(nvtx_b),vtxmap(nvtx_b))

      call occvtx4contr(1,occ_vtx_a,contra,op_info)
      call occvtx4contr(1,occ_vtx_b,contrb,op_info) 

      base1 = ifndmax(occ_vtx_a,1,ngastp*2*nvtx_a,1)
      base1 = max(base1,ifndmax(occ_vtx_b,1,ngastp*2*nvtx_b,1))
      base2 = maxblk_in_contr(contra)
      base2 = max(base2,maxblk_in_contr(contrb))

c dbg
      base1 = base1 + 1
      base2 = base2 + 1
c      print *,'base 1, base 2', base1,base2

c      allocate(ivtx_reo(nvtx_a),fix_vtx(nvtx_a),
c     &     occ_vtx(ngastp,2,nvtx_a))
c      call topo_contr(ieqvfac,reo,ivtx_reo,contra,
c     &     occ_vtx_a,.false.)
c      call canon_contr(contra,reo,ivtx_reo)
c      deallocate(ivtx_reo,fix_vtx,occ_vtx)
c
c      allocate(ivtx_reo(nvtx_b),fix_vtx(nvtx_b),
c     &     occ_vtx(ngastp,2,nvtx_b))
c      call topo_contr(ieqvfac,reo,ivtx_reo,contrb,
c     &     occ_vtx_b,.false.)
c      call canon_contr(contrb,reo,ivtx_reo)
c      deallocate(ivtx_reo,fix_vtx,occ_vtx)
c dbg

      call topomap4contr(2,base1,base2,
     &     topomap_a,vtxinf_a,idum,idum,
     &     contra,occ_vtx_a)
      call topomap4contr(2,base1,base2,
     &     topomap_b,vtxinf_b,idum,idum,
     &     contrb,occ_vtx_b)

c dbg
      if (ntest.ge.100) then
        print *,'A'
        call iwrtma(vtxinf_a,1,nvtx_a,1,nvtx_a)
        call iwrtma(topomap_a,nvtx_a,nvtx_a,nvtx_a,nvtx_a)
        print *,'B'
        call iwrtma(vtxinf_b,1,nvtx_b,1,nvtx_b)
        call iwrtma(topomap_b,nvtx_b,nvtx_b,nvtx_b,nvtx_b)
      end if
c dbg

      call identify_vertices3(vtxmap,contr_in_contr,
     &     vtxinf_a,topomap_a,nvtx_a,
     &     vtxinf_b,topomap_b,nvtx_b)


      if (ntest.ge.100) then
        write(luout,*) 'result: ',contr_in_contr
        write(luout,*) 'vtxmap: ',vtxmap
      end if

      deallocate(occ_vtx_a,occ_vtx_b,vtxinf_a,vtxinf_b,
     &     topomap_a,topomap_b,vtxmap)

      return
      end

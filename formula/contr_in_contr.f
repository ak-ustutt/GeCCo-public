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
     &     ntest = 100

      type(contraction), intent(in) ::
     &     contra, contrb
      type(operator_info), intent(in) ::
     &     op_info

      integer ::
     &     nvtx_a, nvtx_b, narc_a, narc_b,
     &     nj_a, nj_b, lenlist, ivtx, jvtx, nvtx_int, nvtx_new

      integer, pointer ::
c     &     occ_vtx_a(:,:,:), occ_vtx_b(:,:,:),
c     &     topomap_a(:,:), topomap_b(:,:),
c     &     vtxinf_a(:), vtxinf_b(:), vtxmap(:)
cc dbg
cc      integer ::
cc     &     ieqvfac
cc      integer, pointer ::
cc     &     ivtx_reo(:), occ_vtx(:,:,:)
cc      logical, pointer ::
cc     &     fix_vtx(:)
cc      logical ::
cc     &     reo
cc dbg
     &     svertex_a(:), svertex_b(:),
     &     vtxmap(:), list(:), list_reo(:), ireo(:)
      integer(8), pointer ::
     &     xlines_a(:,:), xlines_b(:,:),
     &     topo_a(:,:), topo_b(:,:),
     &     vtx_a(:), vtx_b(:)

      integer, external ::
     &     maxblk_in_contr, ifndmax, njres_contr

      if (ntest.ge.100) then
        write(luout,*) '--------------------------'
        write(luout,*) 'output from contr_in_contr'
        write(luout,*) '--------------------------'
        write(luout,*) 'vertices A/B: ',contra%nvtx,contrb%nvtx
        write(luout,*) 'arcs A/B:     ',contra%narc,contrb%narc
        call prt_contr2(luout,contra,op_info)
        call prt_contr2(luout,contrb,op_info)
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

c      allocate(occ_vtx_a(ngastp,2,nvtx_a),
c     &         occ_vtx_b(ngastp,2,nvtx_b),
c     &         topomap_a(nvtx_a,nvtx_a), topomap_b(nvtx_b,nvtx_b),
c     &         vtxinf_a(nvtx_a), vtxinf_b(nvtx_b),vtxmap(nvtx_b))
c
c      call occvtx4contr(1,occ_vtx_a,contra,op_info)
c      call occvtx4contr(1,occ_vtx_b,contrb,op_info) 
c
c      base1 = ifndmax(occ_vtx_a,1,ngastp*2*nvtx_a,1)
c      base1 = max(base1,ifndmax(occ_vtx_b,1,ngastp*2*nvtx_b,1))
c      base2 = maxblk_in_contr(contra)
c      base2 = max(base2,maxblk_in_contr(contrb))
cc dbg
c      base1 = base1 + 1
c      base2 = base2 + 1
cc      print *,'base 1, base 2', base1,base2
cc      allocate(ivtx_reo(nvtx_a),fix_vtx(nvtx_a),
cc     &     occ_vtx(ngastp,2,nvtx_a))
cc      call topo_contr(ieqvfac,reo,ivtx_reo,contra,
cc     &     occ_vtx_a,.false.)
cc      call canon_contr(contra,reo,ivtx_reo)
cc      deallocate(ivtx_reo,fix_vtx,occ_vtx)
cc
cc      allocate(ivtx_reo(nvtx_b),fix_vtx(nvtx_b),
cc     &     occ_vtx(ngastp,2,nvtx_b))
cc      call topo_contr(ieqvfac,reo,ivtx_reo,contrb,
cc     &     occ_vtx_b,.false.)
cc      call canon_contr(contrb,reo,ivtx_reo)
cc      deallocate(ivtx_reo,fix_vtx,occ_vtx)
cc dbg
c
c      call topomap4contr(2,base1,base2,
c     &     topomap_a,vtxinf_a,idum,idum,
c     &     contra,occ_vtx_a)
c      call topomap4contr(2,base1,base2,
c     &     topomap_b,vtxinf_b,idum,idum,
c     &     contrb,occ_vtx_b)

      nj_a = njres_contr(contra)
      nj_b = njres_contr(contrb)

      allocate(svertex_a(nvtx_a),
     &         svertex_b(nvtx_b),
     &         topo_a(nvtx_a,nvtx_a), topo_b(nvtx_b,nvtx_b),
     &         vtx_a(nvtx_a), vtx_b(nvtx_b),vtxmap(nvtx_b),
     &         xlines_a(nvtx_a,nj_a),xlines_b(nvtx_b,nj_b),
     &         list(nvtx_b*(nvtx_b+1)),ireo(nvtx_b),
     &         list_reo(nvtx_b*(nvtx_b+1)))

      call pack_contr(svertex_a,vtx_a,topo_a,xlines_a,contra,nj_a)
      call pack_contr(svertex_b,vtx_b,topo_b,xlines_b,contrb,nj_b)

      if (ntest.ge.100) then
        print *,'A'
        call prt_contr_p(luout,svertex_a,vtx_a,topo_a,
     &       xlines_a,nvtx_a,nj_a)
        print *,'B'
        call prt_contr_p(luout,svertex_b,vtx_b,topo_b,
     &       xlines_b,nvtx_b,nj_b)
      end if

      call identify_vertices_i8(vtxmap,contr_in_contr,
     &                       vtx_a,topo_a,nvtx_a,
     &                       vtx_b,topo_b,nvtx_b)

      if (ntest.ge.100) then
        write(luout,*) 'result (prel.): ',contr_in_contr
        write(luout,*) 'vtxmap: ',vtxmap
      end if

      if (contr_in_contr.and.nvtx_a.gt.nj_a) then

        ! extra test: can all vertices of a be merged to the
        ! correct number ?

        lenlist = 0
        do ivtx = 1, nvtx_b
          if (vtxmap(ivtx).eq.0) cycle
          do jvtx = ivtx, nvtx_b
            if (vtxmap(jvtx).eq.0) cycle
            lenlist = lenlist+1
            list(1+(lenlist-1)*2) = ivtx
            list(2+(lenlist-1)*2) = jvtx
          end do
        end do
        
        call topo_remove_arcs(topo_b,nvtx_b,list,lenlist)
c dbg
c        print *,'I:'
c        call prt_contr_p(luout,svertex_b,vtx_b,topo_b,
c     &       xlines_b,nvtx_b,nj_b)
c dbg
        
        lenlist = lenlist*2
        call unique_list(list,lenlist)

        call topo_approach_vtxs(ireo,
     &       svertex_b,vtx_b,topo_b,xlines_b,
     &       nvtx_b,nj_b,list,lenlist)

        do ivtx = 1, lenlist
          list_reo(ivtx) = ireo(list(ivtx))
        end do
        call unique_list(list_reo,lenlist)
c dbg
c        print *,'II:'
c        call prt_contr_p(luout,svertex_b,vtx_b,topo_b,
c     &       xlines_b,nvtx_b,nj_b)
c dbg

        call topo_merge_vtxs(ireo,nvtx_new,nvtx_int,
     &                     topo_b,xlines_b,nvtx_b,nj_b,
     &                     list_reo,lenlist)

        contr_in_contr = nvtx_int.le.nj_a

c      call identify_vertices3(vtxmap,contr_in_contr,
c     &     vtxinf_a,topomap_a,nvtx_a,
c     &     vtxinf_b,topomap_b,nvtx_b)

        if (ntest.ge.100)
     &       write(luout,*) '> ',nvtx_int, nj_a

      end if

      if (ntest.ge.100) then
        write(luout,*) 'result: ',contr_in_contr
        write(luout,*) 'vtxmap: ',vtxmap
      end if

      deallocate(xlines_a,xlines_b,topo_a,topo_b,
     &     list,ireo,vtxmap,vtx_a,vtx_b,list_reo)

      return
      end

*----------------------------------------------------------------------*
      subroutine reduce_contr2(sh_sign,iocc_op1op2,njoined_op1op2,
     &     ireo_vtx_no,ireo_vtx_on,ireo0,
     &     ivtx_op1op2,nvtx_red,
     &     mergemap,ld_mmap, 
     &     make_contr_red,contr_red,idxnew_op1op2,
     &     contr,isvtx1,isvtx2,arc_list,nlist,njoined_res)
*----------------------------------------------------------------------*
*     new version of reduce_contr using the topo-representation of
*     contractions for much more straight-forward processing
*
*     given a contraction on contr, and the supervertices to contract
*       (currently redundant additional info: list of arcs involved)
*     figure out the (formal) result of the contraction (including
*     mappings for multi-component operators) and return (if requested)
*     the "reduced" contraction on contr_red
*
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'ifc_operators.h'
      include 'def_contraction.h'
      include 'def_reorder_info.h'
      include 'mdef_operator_info.h'
      include 'def_orbinf.h'
      include 'multd2h.h'

      integer, parameter ::
     &     ntest = 100

      type(contraction), intent(in) ::
     &     contr
      type(contraction), intent(out) ::
     &     contr_red
      logical, intent(in) ::
     &     make_contr_red
      integer, intent(in) ::
     &     njoined_res, idxnew_op1op2, isvtx1, isvtx2, nlist,
     &     arc_list(nlist), ld_mmap
      integer, intent(out) ::
     &     sh_sign,
     &     nvtx_red,
     &     njoined_op1op2,
     &     iocc_op1op2(ngastp,2,*),
     &     mergemap(ld_mmap,2,*),
     &     ireo0(*), ireo_vtx_no(*), ireo_vtx_on(*), ivtx_op1op2(*)

      integer(8), pointer ::
     &     vtx(:), topo(:,:), xlines(:,:), scr(:),
     &     vtx_new(:), topo_new(:,:), xlines_new(:,:), op1op2(:)
      integer, pointer ::
     &     svertex(:), vtx_list(:),vtx_list_reo(:),vtx_list_new(:),
     &     ireo2(:),
     &     svertex_new(:), svertex_reo(:)

      integer ::
     &     nvtx, nvtx_new, nvtx_op1op2, nvtx_cnt, ivtx, idx

      integer, external ::
     &     idxlist

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'reduce_contr2')
        write(luout,*) 'contraction on entry:'
        call prt_contr3(luout,contr,-1)
      end if

      nvtx = contr%nvtx
      allocate(vtx(nvtx), topo(nvtx,nvtx), xlines(nvtx,njoined_res),
     &         scr(nvtx), ireo2(nvtx))
      allocate(vtx_list(max(nvtx**2,nlist*2)),
     &         vtx_list_reo(max(nvtx**2,nlist*2)),
     &         vtx_list_new(max(nvtx**2,nlist*2)) )
      allocate(svertex(nvtx),svertex_reo(nvtx))

      call pack_contr(svertex,vtx,topo,xlines,contr,njoined_res)

      if (ntest.ge.100) then
        write(luout,*) 'contraction in topo form:'
        call prt_contr_p(luout,svertex,vtx,topo,xlines,nvtx,njoined_res)
      end if

      ! get list of contracted vertex pairs
      call set_cnt_vtx_list(vtx_list,contr,arc_list,nlist)

      ! carry out contraction (formally): remove arcs
      call topo_remove_arcs(topo,nvtx,vtx_list,nlist)

      if (ntest.ge.1000) then
        write(luout,*) 'after removing contracted arcs:'
        call prt_contr_p(luout,svertex,vtx,topo,xlines,nvtx,njoined_res)
      end if

      ! augment with joined super vertices, if necessary, and
      ! reduce to unique list of vertices involved in contraction
      nvtx_cnt = 2*nlist
      call make_vtxlist(vtx_list,nvtx_cnt,svertex,nvtx)
      call unique_list(vtx_list,nvtx_cnt) ! contains actual length on exit
      
      if (ntest.ge.1000)
     &     write(luout,*) 'involved vtx: ',vtx_list(1:nvtx_cnt)

      ! move contracted vertices as close together as possible
      ! report reordering in ireo array
      ! idx_new = ireo(idx_old)
      svertex_reo = svertex
      call topo_approach_vtxs(ireo0,sh_sign,
     &     svertex_reo,vtx,topo,xlines,
     &     nvtx,njoined_res,vtx_list,nvtx_cnt)

      do ivtx = 1, nvtx_cnt        
        idx = ireo0(vtx_list(ivtx)) 
        vtx_list_reo(ivtx) = idx
      end do

      if (ntest.ge.1000)
     &     write(luout,*) 'involved vtx(reo): ',vtx_list_reo(1:nvtx_cnt)

      if (ntest.ge.1000) then
        write(luout,*) 'after approaching contracted vertices'
        call prt_contr_p(luout,svertex_reo,
     &       vtx,topo,xlines,nvtx,njoined_res)
      end if

      ! resort list
      call unique_list(vtx_list_reo,nvtx_cnt)
      if (ntest.ge.1000)
     &     write(luout,*) 'updated vtx(reo): ',vtx_list_reo(1:nvtx_cnt)

      ! merge (=symmetrize) contracted vertices, if possible
      call topo_merge_vtxs(ireo2,nvtx_new,nvtx_op1op2,
     &                     topo,xlines,nvtx,njoined_res,
     &                     vtx_list_reo,nvtx_cnt)

      do ivtx = 1, nvtx_cnt
        idx = ireo2(vtx_list_reo(ivtx))
        vtx_list_new(ivtx) = idx
      end do
      call unique_list(vtx_list_new,nvtx_cnt)

      if (ntest.ge.1000) then
        write(luout,*) 'nvtx, nvtx_new, nvtx_cnt, nvtx_op1op2: ',
     &       nvtx, nvtx_new, nvtx_cnt, nvtx_op1op2
        write(luout,*) 'ireo2: ',ireo2(1:nvtx)
        write(luout,*) 'final: ',vtx_list_new(1:nvtx_cnt)
      end if

      if (nvtx_op1op2.ne.nvtx_cnt)
     &     call quit(1,'reduce_contr2','testing?')

      ! store in new topo array
      allocate(topo_new(nvtx_new,nvtx_new),
     &         xlines_new(nvtx_new,njoined_res),
     &         vtx_new(nvtx_new),svertex_new(nvtx_new))

      call topo_reo(svertex_new,vtx_new,topo_new,xlines_new,nvtx_new,
     &              svertex_reo,vtx,    topo,    xlines,    nvtx,
     &              ireo2, njoined_res)

      call topo_rename_vtxs(svertex_new,vtx_new,
     &     nvtx+1,idxnew_op1op2,0,0,
     &     vtx_list_new,nvtx_new,nvtx_op1op2)

      if (ntest.ge.100) then
        write(luout,*) 'reduced term'
        call prt_contr_p(luout,svertex_new,vtx_new,
     &       topo_new,xlines_new,nvtx_new,njoined_res)
      end if

      allocate(op1op2(nvtx_op1op2))
      ! extract occupation of binary contraction result
      call topo_extract_bcres(op1op2,
     &     svertex_new,topo_new,xlines_new,vtx_list_new,
     &     njoined_res,nvtx_new,nvtx_op1op2)

      ! extract merge-map for binary contraction result
      call mergemap_bcres(mergemap,
     &     ld_mmap,
     &     svertex_reo,isvtx1,isvtx2,xlines_new,
     &     ireo2,vtx_list_new,njoined_res,nvtx_new,nvtx,nvtx_op1op2)

c dbg
c      print *,'svertex_reo',svertex_reo
c      print *,'ireo2',ireo2
c      print *,'vtx_new',vtx_new
c dbg

      ! unpack op1op2 to iocc_op1op2
      call unpack_occ(iocc_op1op2,op1op2,nvtx_op1op2)

      njoined_op1op2 = nvtx_op1op2

      ! unpack updated contraction, if requested
      if (make_contr_red) then
        call init_contr(contr_red)
        ! copy header
        contr_red%idx_res  = contr%idx_res
        contr_red%iblk_res = contr%iblk_res
        contr_red%dagger = contr%dagger
        contr_red%fac = contr%fac
        call unpack_contr(contr_red,
     &                  svertex_new,vtx_new,topo_new,xlines_new,
     &                  nvtx_new,njoined_res)

        if (ntest.ge.100) then
          write(luout,*) 'reduced contraction on exit:'
          call prt_contr3(luout,contr_red,-1)
        end if

      else
        if (ntest.ge.100) then
          write(luout,*) 'no reduced contraction requested'
        end if
      end if

      nvtx_red = nvtx_new
      ! new -> old reo; idx_old = ireo_vtx_no(idx_new)
      do ivtx = 1, nvtx
        idx = ireo2(ivtx)
        ireo_vtx_no(idx) = idxlist(ivtx,ireo0,nvtx,1)
      end do
      ! old -> new reo; idx_new = ireo_vtx_no(idx_old)
      do ivtx = 1, nvtx
        ireo_vtx_on(ivtx) = ireo2(ireo0(ivtx))
      end do
      ivtx_op1op2(1:nvtx_op1op2) = vtx_list_new(1:nvtx_op1op2)
      if (ntest.ge.100) then
        write(luout,*) 'nvtx_red = ',nvtx_red
        write(luout,*) 'ireo_vtx_no = ',ireo_vtx_no(1:nvtx_red)
        write(luout,*) 'ireo_vtx_on = ',ireo_vtx_on(1:nvtx)
      end if

      deallocate(vtx, topo, xlines, scr, ireo2,
     &     vtx_list,vtx_list_reo,vtx_list_new, svertex,
     &     topo_new,xlines_new,vtx_new,svertex_new,svertex_reo,
     &     op1op2)

      return
      end

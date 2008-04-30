      subroutine reduce_contr_info(
     &     irestr_vtx_red,info_vtx_red,
     &     irestr_vtx,    info_vtx,
     &     irestr_op1op2, mst_op1op2, gamt_op1op2,
     &     ireo_vtx,ivtx_op1op2,
     &     nvtx,nvtx_red,nvtx_op1op2,nvtx_res,ngas)

      implicit none

      integer, intent(in) ::
     &     nvtx, nvtx_red, nvtx_op1op2, nvtx_res, ngas,
     &     ireo_vtx(nvtx), ivtx_op1op2(nvtx_op1op2),
     &     irestr_vtx(2,ngas,2,2,nvtx+nvtx_res),
     &     irestr_op1op2(2,ngas,2,2,nvtx_op1op2),
     &     info_vtx(2,nvtx+nvtx_res),
     &     mst_op1op2, gamt_op1op2
      integer, intent(out) ::
     &     irestr_vtx_red(2,ngas,2,2,nvtx+nvtx_res),
     &     info_vtx_red(2,nvtx+nvtx_res)

      integer ::
     &     ivtx, ivtx_old, ivtx_new

      irestr_vtx_red(1:2,1:ngas,1:2,1:2,1:nvtx_res) =
     &    irestr_vtx(1:2,1:ngas,1:2,1:2,1:nvtx_res)
      info_vtx_red(1:2,1:nvtx_res) =
     &    info_vtx(1:2,1:nvtx_res)
      do ivtx = 1, nvtx_red
        ivtx_old = ireo_vtx(ivtx)
        irestr_vtx_red(1:2,1:ngas,1:2,1:2,nvtx_res+ivtx) =
     &      irestr_vtx(1:2,1:ngas,1:2,1:2,nvtx_res+ivtx_old)
        info_vtx_red(1:2,nvtx_res+ivtx) =
     &      info_vtx(1:2,nvtx_res+ivtx_old)
      end do

      do ivtx = 1, nvtx_op1op2
        ivtx_new = ivtx_op1op2(ivtx)
        irestr_vtx_red(1:2,1:ngas,1:2,1:2,nvtx_res+ivtx_new) =
     &   irestr_op1op2(1:2,1:ngas,1:2,1:2,ivtx)
        info_vtx_red(1,nvtx_res+ivtx_new) = mst_op1op2
        info_vtx_red(2,nvtx_res+ivtx_new) = gamt_op1op2
      end do

      return
      end

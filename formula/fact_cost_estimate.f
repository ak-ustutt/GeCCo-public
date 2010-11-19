*----------------------------------------------------------------------*
      subroutine fact_cost_estimate(flops,xmemtot,xmemblk,
     &     cnt_info,
     &     str_info,ngas,nsym)
*----------------------------------------------------------------------*
*     fast estimate of contraction expense:
*     use maximum over string lengthes for MS,GAMMA
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      
      real(8), intent(out) ::
     &     flops, xmemtot, xmemblk
      integer, intent(in) ::
     &     ngas, nsym
      type(contraction_info), intent(in) ::
     &     cnt_info
      type(strinf), intent(in) ::
     &     str_info

      real(8) ::
     &     xlen_ex1, xlen_ex2, xlen_cnt, xlen_o12

      real(8), external ::
     &     string_estimate

      ! get estimate of string length for each occupation
      xlen_ex1 = string_estimate(cnt_info%cinfo_ex1c(1:,1),
     &                           cnt_info%cinfo_ex1c(1:,2),
     &                           cnt_info%ncblk_ex1,str_info%g,nsym)
     &         * string_estimate(cnt_info%cinfo_ex1a(1:,1),
     &                           cnt_info%cinfo_ex1a(1:,2),
     &                           cnt_info%nablk_ex1,str_info%g,nsym)

      xlen_ex2 = string_estimate(cnt_info%cinfo_ex2c(1:,1),
     &                           cnt_info%cinfo_ex2c(1:,2),
     &                           cnt_info%ncblk_ex2,str_info%g,nsym)
     &         * string_estimate(cnt_info%cinfo_ex2a(1:,1),
     &                           cnt_info%cinfo_ex2a(1:,2),
     &                           cnt_info%nablk_ex2,str_info%g,nsym)

      xlen_cnt = string_estimate(cnt_info%cinfo_cntc(1:,1),
     &                           cnt_info%cinfo_cntc(1:,2),
     &                           cnt_info%ncblk_cnt,str_info%g,nsym)
     &         * string_estimate(cnt_info%cinfo_cnta(1:,1),
     &                           cnt_info%cinfo_cnta(1:,2),
     &                           cnt_info%nablk_cnt,str_info%g,nsym)

      xlen_o12 = string_estimate(cnt_info%cinfo_op1op2c(1:,1),
     &                           cnt_info%cinfo_op1op2c(1:,2),
     &                           cnt_info%ncblk_op1op2,str_info%g,nsym)
     &         * string_estimate(cnt_info%cinfo_op1op2a(1:,1),
     &                           cnt_info%cinfo_op1op2a(1:,2),
     &                           cnt_info%nablk_op1op2,str_info%g,nsym)

      ! use these for FLOP and MEM estimates
      flops = xlen_ex1*xlen_ex2*xlen_cnt/(dble(nsym)**2)
      xmemtot = xlen_o12
      xmemblk = xlen_o12/dble(nsym)

      return
      end

      real(8) function string_estimate(cnt_occ,cnt_grph,
     &                                 nblk,grph,nsym)

      implicit none

      include 'def_graph.h'

      integer, intent(in) ::
     &     nblk, nsym, cnt_occ(nblk), cnt_grph(nblk)
      type(graph) ::
     &     grph(*)
      
      integer, pointer ::
     &     lenstr(:,:)
      integer ::
     &     iblk, ims, isym, nms, lenmax
      real(8) ::
     &     xlen

      xlen = 1d0
      do iblk = 1, nblk
        nms = cnt_occ(iblk)+1
        lenmax = 0
        lenstr => grph(cnt_grph(iblk))%lenstr_gm
        do ims = 1, nms
          do isym = 1, nsym
cmh            lenmax = max(lenstr(isym,ims),lenmax)
            lenmax = lenmax + lenstr(isym,ims)
          end do
        end do
        xlen = xlen*dble(lenmax)
      end do

      string_estimate = xlen
      
      return
      end

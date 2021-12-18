*----------------------------------------------------------------------*
      subroutine reduce_vtx_info(irestr_vtx,info_vtx,
     &                           contr,occ_vtx,iarc_red,
     &                           irestr_res,orb_info)
*----------------------------------------------------------------------*
*     analogon to reduce_contr() for updating info on intermediates
*     call either before reduce_contr(), or with a copy of the old
*     contr and occ_vtx
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'def_orbinf.h'
      include 'ifc_operators.h'
      include 'multd2h.h'

      integer, parameter ::
     &     ntest = 00
      
      type(contraction), intent(in) ::
     &     contr
      type(orbinf), intent(in) ::
     &     orb_info
      integer, intent(in) ::
     &     iarc_red, 
     &     occ_vtx(ngastp,2,*),
     &     irestr_res(2,orb_info%ngas,2,2)
      integer, intent(inout) ::
     &     irestr_vtx(2,orb_info%ngas,2,2,*),
     &     info_vtx(2,*)

      integer ::
     &     ivtx1, ivtx2, ivtx, nvtx, ngas, njoined_res,
     &     occ_new(ngastp,2), arc_list(contr%narc), len_list,
     &     iarc_prm

      call quit(1,'reduce_vtx_info','should be obsolete now')

      njoined_res = 1
      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'reduce_vtx_info')
        write(lulog,*) 'on entry:'
        write(lulog,*) 'info_vtx:'
        write(lulog,'(3x,2(i2,x,i2,2x))')
     &       info_vtx(1:2,1:contr%nvtx+njoined_res)
        write(lulog,*) 'irestr_vtx'
        do ivtx = 1, contr%nvtx+njoined_res
          call wrt_rstr(lulog,irestr_vtx(1,1,1,1,ivtx),orb_info%ngas)
        end do
      end if


      ngas = orb_info%ngas
      nvtx = contr%nvtx

      iarc_prm =  iarc_red      !arc_list(ilist)

      ivtx1 = contr%arc(iarc_prm)%link(1)
      ivtx2 = contr%arc(iarc_prm)%link(2)

      info_vtx(1,ivtx1+1) = info_vtx(1,ivtx1+1)+info_vtx(1,ivtx2+1)
      info_vtx(2,ivtx1+1) = multd2h(info_vtx(2,ivtx1+1),
     &                                info_vtx(2,ivtx2+1))
      
      do ivtx = ivtx2, nvtx-1
        info_vtx(1:2,ivtx+1) = info_vtx(1:2,ivtx+2)
      end do

      occ_new = -contr%arc(iarc_prm)%occ_cnt
      occ_new = occ_new + iocc_dagger(occ_new)
      occ_new = occ_new +
     &     occ_vtx(1:ngastp,1:2,ivtx1+1) + occ_vtx(1:ngastp,1:2,ivtx2+1)

      ! still a preliminary thing:
      call fit_restr(irestr_vtx(1,1,1,1,ivtx1+1),occ_new,
     &     irestr_res,orb_info%ihpvgas,orb_info%ngas)

      do ivtx = ivtx2, nvtx-1
        irestr_vtx(1:2,1:ngas,1:2,1:2,ivtx+1) =
     &       irestr_vtx(1:2,1:ngas,1:2,1:2,ivtx+2)
      end do

      if (ntest.ge.100) then
        write(lulog,*) 'on exit:'
        write(lulog,*) 'info_vtx:'
        write(lulog,'(3x,2(i2,x,i2,2x))')
     &       info_vtx(1:2,1:contr%nvtx+njoined_res)
        write(lulog,*) 'irestr_vtx'
        do ivtx = 1, contr%nvtx+1
          call wrt_rstr(lulog,irestr_vtx(1,1,1,1,ivtx),orb_info%ngas)
        end do
      end if

      return
      end

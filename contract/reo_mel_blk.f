*----------------------------------------------------------------------*
      subroutine reo_mel_blk(meinp,meout,iblkinp,iblkout,
     &     str_info,strmap_info,orb_info,ifrom,ito,idxinp)
*----------------------------------------------------------------------*
*     reorder vertex ifrom to vertex ito and write result the list
*     of meout
*     
*     matthias, dec 2009
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_orbinf.h'
      include 'mdef_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'def_reorder_info.h'
      include 'ifc_memman.h'

      integer, parameter ::
     &     ntest = 00

      type(orbinf), intent(in) ::
     &     orb_info
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf), intent(in) ::
     &     strmap_info
      type(me_list), intent(in), target ::
     &     meinp, meout
      integer, intent(in) ::
     &     iblkinp, iblkout, ifrom, ito, idxinp

      integer ::
     &     njinp, njout, idoff_in, idoff_out, idum, imode,
     &     ij, ijinp, ijout, ngas
      real(8) ::
     &     xret_dum
      logical ::
     &     problem

      type(operator), pointer ::
     &     opinp, opout
      type(filinf), pointer ::
     &     ffin, ffout

      type(reorder_info) ::
     &     reo_info

      integer, pointer ::
     &     iocc_inp(:,:,:), iocc_out(:,:,:),
     &     irst_inp(:,:,:,:,:)

      integer, allocatable ::
     &     iocc(:,:,:), iocc_tmp(:,:,:),
     &     irst(:,:,:,:,:), irst_tmp(:,:,:,:,:),
     &     merge_stp1(:), merge_stp2(:), merge_stp1inv(:),
     &     merge_stp2inv(:), svertex(:),occ_vtx(:,:,:),
     &     info_vtx(:,:), irst_vtx(:,:,:,:,:)

      logical, external ::
     &     iocc_zero, iocc_equal

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'reo_mel_blk')
      end if

      opinp => meinp%op
      opout => meout%op

      njinp = opinp%njoined
      njout = opout%njoined

      ngas = orb_info%ngas

      ! translate records to offset in file:
      ffin => meinp%fhand
      idoff_in = ffin%length_of_record*(ffin%current_record-1)
      ffout => meout%fhand
      idoff_out = ffout%length_of_record*(ffout%current_record-1)

      iocc_inp => opinp%ihpvca_occ(1:ngastp,1:2,
     &           (iblkinp-1)*njinp+1:iblkinp*njinp)
      iocc_out => opout%ihpvca_occ(1:ngastp,1:2,
     &           (iblkout-1)*njout+1:iblkout*njout)
      irst_inp => opinp%igasca_restr(1:2,1:ngas,1:2,1:2,1,
     &           (iblkinp-1)*njinp+1:iblkinp*njinp)

      call init_reo_info(reo_info)

      ! set reorder information: shift entire vertex
      allocate(reo_info%reo(2*njinp),reo_info%nca_vtx(njinp))
      reo_info%nvtx_contr = njinp
      call set_nca_vtx(reo_info%nca_vtx,iocc_inp,njinp)

      reo_info%nreo = 1
      reo_info%reo(1)%idxsuper = 1
      reo_info%reo(1)%idxop_ori = idxinp
      reo_info%reo(1)%iblkop_ori = 1
      reo_info%reo(1)%is_bc_result = .true.
      reo_info%reo(1)%reo_before = .false.
      reo_info%reo(1)%shift_i0 = .false.
      reo_info%reo(1)%to = ito
      reo_info%reo(1)%from = ifrom
      reo_info%reo(1)%to_vtx   = ito
      reo_info%reo(1)%from_vtx = ifrom
      reo_info%reo(1)%occ_shift = iocc_inp(1:ngastp,1:2,ifrom)

      ! most arrays here are just dummies
      allocate(iocc(ngastp,2,njinp), iocc_tmp(ngastp,2,njinp),
     &     irst(2,ngas,2,2,njinp),
     &     irst_tmp(2,ngas,2,2,njinp),
     &     merge_stp1(2*njinp*njinp), merge_stp2(2*njinp*njinp),
     &     merge_stp1inv(2*njinp*njinp), merge_stp2inv(2*njinp*njinp),
     &     svertex(njinp),occ_vtx(ngastp,2,2*njinp),
     &     info_vtx(2,2*njinp),irst_vtx(2,ngas,2,2,2*njinp))

      iocc = iocc_inp
      iocc(1:ngastp,1:2,ito) = iocc(1:ngastp,1:2,ito)
     &       + iocc(1:ngastp,1:2,ifrom)
      iocc(1:ngastp,1:2,ifrom) = 0
      irst = irst_inp
      irst(1:2,1:ngas,1:2,1:2,ito) = irst(1:2,1:ngas,1:2,1:2,ito)
     &       + irst(1:2,1:ngas,1:2,1:2,ifrom)
      irst(1:2,1:ngas,1:2,1:2,ifrom) = 0

      ! now there should be a one-to-one correspondence to the output op.
      problem = .false.
      do ij = 1, njout
        if (.not.iocc_zero(iocc_out(1:ngastp,1:2,ij))) then
          ijout = ij
          exit
        end if
      end do
      do ijinp = 1, njinp
        if (iocc_zero(iocc(1:ngastp,1:2,ijinp))) cycle
        problem = problem.or.
     &            .not.iocc_equal(iocc_out(1:ngastp,1:2,ijout),
     &              .false.,iocc(1:ngastp,1:2,ijinp),.false.)
        do ij = ijout+1, njout
          if (.not.iocc_zero(iocc_out(1:ngastp,1:2,ij))) then
            ijout = ij
            exit
          end if
        end do
      end do
      do ij = ijout+1, njout
        problem = problem.or..not.iocc_zero(iocc_out(1:ngastp,1:2,ij))
      end do
      if (problem) call quit(1,'reo_mel_blk','ops do not match?')

      svertex = 1
      occ_vtx(1:ngastp,1:2,1:njinp) = iocc_inp
      occ_vtx(1:ngastp,1:2,njinp+1:2*njinp) = iocc
      irst_vtx(1:2,1:ngas,1:2,1:2,1:njinp) = irst_inp
      irst_vtx(1:2,1:ngas,1:2,1:2,njinp+1:2*njinp) = irst

      do imode = 1, 2
        call get_reo_info2(imode,1,
     &         iocc,iocc_tmp,
     &         irst,irst_tmp,
     &         njinp,idum,idum,
     &         merge_stp1,merge_stp1inv,merge_stp2,merge_stp2inv,
     &         occ_vtx,irst_vtx,svertex,info_vtx,njinp,njinp,
     &         reo_info,reo_info%nreo,str_info,orb_info)
      end do
c dbg
c      if (trim(meout%label).ne.'ME_Dproj'
c     &    .and.reo_info%sign_reo.ne.1d0) then
      if (reo_info%sign_reo.ne.1d0) then
        write(luout,*) 'setting sign_reo = +1 for block ',iblkout
        reo_info%sign_reo = 1d0
      end if
c dbgend

      call reo_op_wmaps_c(
     &     .false.,xret_dum,0,
     &     meinp,meout,
     &     .false., .false.,
     &     iblkinp,iblkout,
     &     idoff_in,idoff_out,
     &     reo_info,
     &     str_info,strmap_info,orb_info)

      call dealloc_reo_info(reo_info)

      deallocate(iocc,iocc_tmp,
     &     irst,irst_tmp,merge_stp1,merge_stp2,
     &     merge_stp1inv,merge_stp2inv,
     &     svertex,occ_vtx,info_vtx,irst_vtx)

      return
      end

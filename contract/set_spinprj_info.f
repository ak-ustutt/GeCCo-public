      subroutine set_spinprj_info(offsets,maps_c,maps_a,
     &           ldim_op_c,ldim_op_a,gsign,
     &           nsplc,mel,iblk,ioff0,igama,ngam,graphs,
     &           msdis_c,gamdis_c,hpvx_csub,occ_csub,graph_csub,ncblk,
     &           msdis_a,gamdis_a,hpvx_asub,occ_asub,graph_asub,nablk,
     &           s2,ms2,ncoup,coeff)

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_graph.h'
      include 'hpvxseq.h'

      integer, parameter ::
     &    ntest = 00

      integer, intent(in) ::
     &    nsplc, iblk, ioff0, igama, ngam, ncblk, nablk,
     &    s2, ms2, ncoup
      integer, intent(in) ::
     &    msdis_c(ncblk), gamdis_c(ncblk), hpvx_csub(ncblk), 
     &    occ_csub(ncblk), graph_csub(ncblk),
     &    msdis_a(nablk), gamdis_a(nablk), hpvx_asub(nablk), 
     &    occ_asub(nablk), graph_asub(nablk)
      type(me_list), intent(inout) ::
     &     mel
      type(graph), intent(inout) ::
     &    graphs(*)
      integer, intent(out) ::
     &    gsign(nsplc), offsets(nsplc), 
     &    maps_c(nsplc,ncblk), maps_a(nsplc,nablk),
     &    ldim_op_c(ncblk,nsplc), ldim_op_a(nablk,nsplc)
      real(8), intent(out) ::
     &    coeff(nsplc,ncoup)

      integer ::
     &    isplc, idxmsa, idxdis, occ, msc, msa
      integer ::
     &    msdiscmp_c(ncblk,nsplc), msdiscmp_a(nablk,nsplc),
     &    len_str(ncblk+nablk), idxmsdis_c(ncblk), idxmsdis_a(nablk)

      integer, external ::
     &    std_spsign_msdis, ielsum, idx_msgmdst2

      occ = ielsum(occ_asub,nablk)
      if (occ.ne.ielsum(occ_csub,ncblk)) 
     &     call quit(1,'set_spinprj_info','occ(A)/=occ(C) not allowed!')

      ! for the cases below, the spin distributions of the
      ! first sceleton element differ a bit, cf. also spin_prj_blk
      if (occ.eq.1) then
        call set_case((/+1/),(/+1/),1)
      else if (occ.eq.2.and.ncblk.eq.1.and.nablk.eq.1) then
        call set_case((/-1,+1/),(/-1,+1/),2)
      else if (occ.eq.2.and.ncblk.eq.2.and.nablk.eq.1) then
        call set_case((/+1,-1/),(/-1,+1/),2)
      else if (occ.eq.2.and.ncblk.eq.1.and.nablk.eq.2) then
        call set_case((/-1,+1/),(/+1,-1/),2)
      else if (occ.eq.2.and.ncblk.eq.2.and.nablk.eq.2) then
        call set_case((/+1,-1/),(/+1,-1/),2)
      else if (occ.eq.3.and.ncblk.eq.1.and.nablk.eq.1) then
        call set_case((/-1,+1,+1/),(/-1,+1,+1/),3)
      else if (occ.eq.3.and.ncblk.eq.2.and.nablk.eq.2.and.
     &         occ_csub(1).eq.2.and.occ_asub(1).eq.1) then
        call set_case((/-1,+1,+1/),(/+1,-1,+1/),3)
      else if (occ.eq.3.and.ncblk.eq.2.and.nablk.eq.2.and.
     &         occ_csub(1).eq.1.and.occ_asub(1).eq.2) then
        call set_case((/+1,-1,+1/),(/-1,+1,+1/),3)
      else
        write(luout,*) 'occ, ncblk, nablk = ',occ, ncblk, nablk
        write(luout,*) 'occ_csub: ',occ_csub(1:ncblk)
        write(luout,*) 'occ_asub: ',occ_asub(1:nablk)
        call quit(1,'set_spinprj_info','case not covered')
      end if

      do isplc = 1, nsplc

        if (ntest.ge.100) then
          write(luout,*) 'generated spin-case and maps: case = ',isplc
          write(luout,*) 'msdiscmp_c = ',msdiscmp_c(1:ncblk,isplc)
          write(luout,*) 'msdiscmp_a = ',msdiscmp_a(1:nablk,isplc)
          write(luout,*) 'maps_c = ',maps_c(isplc,1:ncblk)
          write(luout,*) 'maps_a = ',maps_a(isplc,1:nablk)
          write(luout,'(a,10f10.6)') 'coupling coefficients: ',
     &                               coeff(isplc,1:ncoup)
        end if

        gsign(isplc)=
     &         std_spsign_msdis(msdiscmp_c(1,isplc),occ_csub,ncblk)
        gsign(isplc)=gsign(isplc)*
     &         std_spsign_msdis(msdiscmp_a(1,isplc),occ_asub,nablk)

        call ms2idxms(idxmsdis_c,msdiscmp_c(1,isplc),occ_csub,ncblk)
        call ms2idxms(idxmsdis_a,msdiscmp_a(1,isplc),occ_asub,nablk)

        msa = ielsum(msdiscmp_a(1,isplc),nablk)
        idxmsa = (occ-msa)/2+1

        idxdis = 1
        if (mel%off_op_gmox(iblk)%ndis(igama,idxmsa).gt.1)
     &           idxdis =
     &               idx_msgmdst2(.false., !don't give error
     &                iblk,idxmsa,igama,
     &                occ_csub,idxmsdis_c,gamdis_c,ncblk,
     &                occ_asub,idxmsdis_a,gamdis_a,nablk,
     &                .false.,-1,-1,mel,ngam)
        if (idxdis.lt.0) then
          offsets(isplc) = 0
        else
        offsets(isplc) = mel%off_op_gmox(iblk)%
     &             d_gam_ms(idxdis,igama,idxmsa) - ioff0
        end if

        call set_len_str(len_str,ncblk,nablk,
     &                   graphs,
     &                   graph_csub,idxmsdis_c,gamdis_c,hpvx_csub,
     &                   graph_asub,idxmsdis_a,gamdis_a,hpvx_asub,
     &                   hpvxseq,.false.)

        call set_op_ldim_c(ldim_op_c(1,isplc),ldim_op_a(1,isplc),
     &           hpvx_csub,hpvx_asub,
     &           len_str,ncblk,nablk,.false.)

        if (ntest.ge.100) then
          write(luout,*) 'offset = ',offsets(isplc)
          write(luout,*) 'ldim_op_c',ldim_op_c(1:ncblk,isplc)
          write(luout,*) 'ldim_op_a',ldim_op_a(1:nablk,isplc)
          write(luout,*) 'gsign = ',gsign(isplc)
        end if

      end do

      contains

      subroutine set_case(spin_ref_c,spin_ref_a,npart)

      integer, intent(in) ::
     &     npart, spin_ref_c(npart), spin_ref_a(npart)

      integer ::
     &     idx_spc, len, ioff, idx,
     &     mind(npart), maxd(npart),
     &     spin_dis_a(npart), spin_dis_c(npart)

      logical, external ::
     &     next_dist2

      mind(1:npart) = -1 
      maxd(1:npart) = +1

      idx_spc = 0
      spin_dis_a(1:npart) = maxd(1:npart)
      ! enumerate spin cases
      aloop: do
        msa = sum(spin_dis_a(1:npart))
        spin_dis_c(1:npart) = maxd(1:npart)
        cloop: do
          msc = sum(spin_dis_c(1:npart))
          if (msc.eq.msa) then
            idx_spc = idx_spc+1
            ! collect ms according to subblock structure
            ioff = 0
            do idx = 1, ncblk
              len = occ_csub(idx)
              ! set msdiscmp_c
              msdiscmp_c(idx,idx_spc) = 
     &             sum(spin_dis_c(ioff+1:ioff+len))
              ! set maps_c
              maps_c(idx_spc,idx) = 
     &             map_type(spin_dis_c(ioff+1:ioff+len),
     &             spin_ref_c(ioff+1:ioff+len),len)
              ioff = ioff+len
            end do
            ioff = 0
            do idx = 1, nablk
              len = occ_asub(idx)
              ! set msdiscmp_a
              msdiscmp_a(idx,idx_spc) = 
     &             sum(spin_dis_a(ioff+1:ioff+len))
              ! set maps_a
              maps_a(idx_spc,idx) = 
     &             map_type(spin_dis_a(ioff+1:ioff+len),
     &             spin_ref_a(ioff+1:ioff+len),len)
              ioff = ioff+len
            end do

            ! Get genealogical coupling coefficients for this string
            call get_genealogical_coeff(s2,ms2,npart,npart,ncoup,
     &               coeff(idx_spc,1:ncoup),spin_dis_c,spin_dis_a)

          end if
          if (.not.next_dist2(spin_dis_c,npart,mind,maxd,-2)) exit cloop
        end do cloop
        if (.not.next_dist2(spin_dis_a,npart,mind,maxd,-2)) exit aloop
      end do aloop
          
      end subroutine

      integer function map_type(distgt,disref,npart)

      integer, intent(in) ::
     &     npart, distgt(npart), disref(npart)
      integer ::
     &     res

      res=-1
      select case(npart)
      case(1)
        if (distgt(1).eq.disref(1)) res = 0
        if (distgt(1).ne.disref(1)) res = 1
      case(2)
        if (disref(1).eq.-1.and.disref(2).eq.+1) then
          if (distgt(1).eq.+1) then
            if (distgt(2).eq.+1) res = 1
            if (distgt(2).eq.-1) res = 2
          else if (distgt(1).eq.-1) then
            if (distgt(2).eq.+1) res = 0
            if (distgt(2).eq.-1) res = 3
          end if
        end if
      case(3)
        ! an idea is welcome ...
      end select
c dbg
c      print *,'tgt = ',distgt(1:npart)
c      print *,'ref = ',disref(1:npart)
c      print *,' -> res = ',res
c dbg

      if (res.lt.0) call quit(1,'set_spinprj_info','error in map_type')

      map_type = res

      return

      end function

      end subroutine

*----------------------------------------------------------------------*
      logical function nondia_blk(mapca,diag_idx,diag_ca,
     &                            iocc,nj,nc,na,diag_type)
*----------------------------------------------------------------------*
*     determines whether an operator block is diagonal.
*     if diagonal, a map for corresponding indices and
*     a list containing indices of the c/a (sub)strings which
*     characterizes the diagonal elements is returned.
*
*     matthias, feb 2010
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'hpvxseq.h'

      integer, intent(out) ::
     &     mapca(nc), diag_idx(nc), diag_ca(nc)
      integer, intent(in) ::
     &     iocc(ngastp,2,nj), nj, nc, na, diag_type

      integer ::
     &     occ_csub(nc), occ_asub(na), hpvx_csub(nc), hpvx_asub(na),
     &     mapac(na), idx


      ! set HPVX and OCC info
      call condense_occ(occ_csub, occ_asub,
     &                  hpvx_csub,hpvx_asub,
     &                  iocc,nj,hpvxblkseq)

      ! get transposition maps
      call set_dis_tra_map(mapca,mapac,
     &          hpvx_csub,hpvx_asub,nc,na)
      nondia_blk = .false.
      ! only diagonal blocks
      if (nc.eq.na) then
        do idx = 1, nc
          if (occ_csub(idx).ne.occ_asub(mapca(idx)).or.
     &        hpvx_csub(idx).ne.hpvx_asub(mapca(idx))) then
            nondia_blk = .true.
            exit
          end if
        end do
        if (.not.nondia_blk) then
          if (diag_type.ne.1) call quit(1,'nondia_blk',
     &          'only diag_type=1 available so far')
          ! get index tuples which constitute diagonal
          call get_diag_tuple1(diag_idx,diag_ca,
     &                         iocc,nj,hpvxblkseq)
        end if
      else
        nondia_blk = .true.
      end if

      return
      end
*----------------------------------------------------------------------*

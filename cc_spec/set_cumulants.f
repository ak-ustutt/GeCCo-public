*----------------------------------------------------------------------*
      subroutine set_cumulants(form_cum,
     &     title,label,nlabels,mode,level,op_info)
*----------------------------------------------------------------------*
*     defines cumulants in terms of reduced densities (mode='CUMULANT')
*     or reduced densities in terms of cumulants (default)
*
*     matthias feb 2010
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      include 'def_formula_item.h'
      include 'def_formula.h'
      include 'ifc_input.h'

      type(formula), intent(inout), target ::
     &     form_cum
      integer, intent(in) ::
     &     nlabels, level
      character(*), intent(in) ::
     &     label(nlabels), title, mode

      type(operator_info), intent(inout) ::
     &     op_info

      ! local variables
      character ::
     &     name*(form_maxlen_label*2)

      type(formula_item), target ::
     &     flist
      type(formula_item), pointer ::
     &     flist_pnt

      integer ::
     &     ilabel, idx1, idxd, idxscal,
     &     idx, nn, ii, jj, iblk, mini, maxc, off, max_n

      logical ::
     &     cum

      real(8) ::
     &     fac

      integer, allocatable ::
     &     idx_op(:), iblk_min(:), iblk_max(:), avoid(:,:),
     &     idx_op_vtx(:), dist(:), imnmx(:,:)
 
      integer, external::
     &     idx_oplist2, ifac, imltlist
      logical, external ::
     &     next_dist

      ! for timings:
      real(8) ::
     &     cpu, wall, sys, cpu0, wall0, sys0

      off = 0 ! set to 1 if first occ. cls. of idx1, idxd is scalar

      if (ntest.eq.100) then
        call write_title(luout,wst_dbg_subr,
     &       'output from set_cumulants')
        write(luout,*) ' mode = ',trim(mode)
      end if

      cum = mode(1:8).eq.'CUMULANT'

      call atim_csw(cpu0,sys0,wall0)

      ! transform labels into indices in this way
      do ilabel = 1, nlabels
        idx = idx_oplist2(label(ilabel),op_info)
        if (idx.le.0)
     &       call quit(1,'set_cumulants',
     &       'label not on list: '//trim(label(ilabel)))
        if (ilabel.eq.1) idxscal = idx
        if (ilabel.eq.2) idx1 = idx
        if (ilabel.eq.3) idxd = idx
      end do

      maxc = op_info%op_arr(idxd)%op%n_occ_cls-off
      ! exclude formal classes
      do idx = op_info%op_arr(idxd)%op%n_occ_cls, 1, -1
        if (.not.op_info%op_arr(idxd)%op%formal_blk(idx)) exit
        maxc = maxc - 1
      end do

      ! initialize formula
      call init_formula(flist)
      flist_pnt => flist
      ! put [INIT] at the beginning
      call new_formula_item(flist_pnt,command_set_target_init,idxscal)
      flist_pnt => flist_pnt%next

      fac = 1d0

      ! scalar element if requested
      if (off.eq.1) then
        call expand_op_product2(flist_pnt,idxscal,
     &       fac,3,2,
     &       (/idxd,idx1,idxd/),
     &       (/2,3,2/),
     &       (/1,1,1/),(/1,1,1/),
     &       0,0,
     &       0,0,
     &       0,0,
     &       .false.,op_info)
        do while(flist_pnt%command.ne.command_end_of_formula)
          flist_pnt => flist_pnt%next
        end do
      end if

      if (level.gt.0) then
        max_n = min(level,op_info%op_arr(idx1)%op%n_occ_cls-off)
      else
        max_n = op_info%op_arr(idx1)%op%n_occ_cls-off
      end if
      do nn = 1, max_n

        if (cum) then
          fac = dble(((-1)**(nn-1))*ifac(nn-1))
        end if        

        allocate(idx_op(2*nn+1),iblk_min(2*nn+1),iblk_max(2*nn+1),
     &           avoid(1:2,nn**2),idx_op_vtx(2*nn+1),dist(nn),
     &           imnmx(2,nn))
        do ii = 1, 2*nn+1
          idx_op(ii) = nn+2-abs(nn+1-ii)
          if (ii.eq.nn+1) then
            idx_op_vtx(ii) = idx1
          else
            idx_op_vtx(ii) = idxd
          end if
        end do
        ! no contractions between densities
        idx = 0
        do ii = 1, nn
          do jj = nn+2,2*nn+1
            idx = idx + 1
            avoid(1,idx) = ii
            avoid(2,idx) = jj
          end do
        end do

        dist = 1
        dist(1) = 0
        imnmx(1,1:nn) = 1
        imnmx(2,1:nn) = maxc

        do while(next_dist(dist,nn,imnmx,1))
          iblk = sum(dist)
          if (iblk.gt.op_info%op_arr(idx1)%op%n_occ_cls) then
            mini = 0
            do ii = nn-1,1,-1
              if (mini.eq.0.and.dist(ii).eq.maxc) mini=dist(ii+1)
              if (mini.gt.0) imnmx(1,ii) = mini+1
            end do
            cycle
          end if
          do ii = 1, 2*nn+1
            if (ii.eq.nn+1) then
              iblk_min(ii) = iblk+off
              iblk_max(ii) = iblk+off
            else
              iblk_min(ii) = dist(nn+1-abs(nn+1-ii))+off
              iblk_max(ii) = dist(nn+1-abs(nn+1-ii))
            end if
          end do
c dbg testing prefactor correction (need change in topo_contr)
          fac = 1d0
          if (cum) then
            fac = dble(((-1)**(nn-1))*ifac(nn-1))
          end if
          do ii = 1, maxc
            fac = fac*ifac(imltlist(ii,dist,nn,1))
          end do
c dbgend
c dbg
c          print *,'nn, iblk: ',nn, iblk
c          print *,'idx_op_vtx: ',idx_op_vtx
c          print *,'idx_op    : ',idx_op
c          print *,'iblk_min  : ',iblk_min
c          print *,'iblk_max  : ',iblk_max
c          print *,'avoid: ',avoid
c          print *,'fac: ',fac
c dbgend  

          ! expand <0|D D D ... 1 ... D D D|0>
          call expand_op_product2(flist_pnt,idxscal,
     &         fac,2*nn+1,nn+1,
     &         idx_op_vtx,
     &         idx_op,
     &         iblk_min,iblk_max,
     &         0,0,
     &         avoid,nn**2,
     &         0,0,
     &         .false.,op_info)

          do while(flist_pnt%command.ne.command_end_of_formula)
            flist_pnt => flist_pnt%next
          end do
          mini = 0
          do ii = nn-1,1,-1
            if (mini.eq.0.and.dist(ii).eq.imnmx(2,1)) mini=dist(ii+1)
            if (mini.gt.0) imnmx(1,ii) = mini+1
          end do
        end do
        deallocate(idx_op_vtx,idx_op,iblk_min,iblk_max,avoid,
     &             dist,imnmx)
      end do

      if (ntest.ge.100) then
        call write_title(luout,wst_title,'Final formula')
        call print_form_list(luout,flist,op_info)
      end if

      ! assign comment
      form_cum%comment = trim(title)
      ! write to disc
      write(name,'(a,".fml")') trim(form_cum%label)
      call file_init(form_cum%fhand,name,ftyp_sq_unf,0)
      call write_form_list(form_cum%fhand,flist,
     &     form_cum%comment)

      call dealloc_formula_list(flist)

      call atim_csw(cpu,sys,wall)
      call prtim(luout,'set_cumulants',cpu-cpu0,sys-sys0,wall-wall0)

      return
      end

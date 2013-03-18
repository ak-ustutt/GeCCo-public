*----------------------------------------------------------------------*
      subroutine factorize_new(fl_fact,fl_raw,
     &     op_info,str_info,orb_info,label)
*----------------------------------------------------------------------*
*     convert formula on linked list fl_raw to sequence of operations
*     (i.pt. binary contractions) on fl_fact
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'mdef_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'

      integer, parameter ::
     &     ntest = 00

      type(operator_info), intent(inout) ::
     &     op_info
      type(strinf), intent(inout) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info
      type(formula_item), intent(in), target ::
     &     fl_raw
      type(formula_item), intent(inout), target ::
     &     fl_fact
      character*(*), intent(in) ::
     &     label

      integer, parameter ::
     &     max_stat = 1000000
      type(formula_item), pointer ::
     &     fl_ptr, fl_fact_ptr
      integer ::
     &     iterm, nterms, istat, binmx(ngastp), binmxtmp(ngastp),
     &     ibin1, ibin2, ibin3, ibin4, iitem
      integer, pointer ::
     &     iscale_stat(:,:,:), ireo_t(:), ireo_m(:), ireo_s(:),
     &     binning(:,:,:,:), tmp(:,:,:,:)
      real(8), pointer ::
     &     time_stat(:), mem_stat(:), scale_stat(:)
      real(8) ::
     &     cpu0, sys0, wall0, cpu, sys, wall, xsum
      type(filinf) ::
     &     ffstat

      logical, external ::
     &     check_contr4zeroop

      real(8), external ::
     &     scale_rank

      call atim_csw(cpu0,sys0,wall0)

      call write_title(luout,wst_subsection,
     &     'Formula factorization')

      allocate(iscale_stat(ngastp,2,max_stat),
     &         time_stat(max_stat),mem_stat(max_stat),
     &         scale_stat(max_stat))

      iscale_stat = 0
      time_stat   = 0d0
      mem_stat    = 0d0

      if (iprlvl.ge.3) then
        binmx = (/4,4,0,0/) ! just initial, could also set all to 0
        allocate(binning(binmx(1)+1,binmx(2)+1,binmx(3)+1,binmx(4)+1))
        binning = 0
      end if

      fl_ptr => fl_raw
      fl_fact_ptr => fl_fact
      iterm = 0
      iitem = 0
      do while(fl_ptr%command.ne.command_end_of_formula)
        if (fl_ptr%command.eq.command_add_contribution) then
          iterm = iterm + 1
          if (iprlvl.ge.10)
     &         write(luout,*) 'factorizing term # ',iterm
          if (lustat.gt.0)
     &         write(lustat,*) 'factorizing term # ',iterm
          if (iterm.eq.max_stat+1)
     &         call warn('factorize_new','max_stat exceeded')
          istat = min(iterm,max_stat)

          if (.not.check_contr4zeroop(fl_ptr%contr,op_info)) then 
           call form_fact_new(fl_fact_ptr,fl_ptr%contr,
     &       op_info,str_info,orb_info,
     &       iscale_stat(1,1,istat),time_stat(istat),mem_stat(istat),
     &       iitem)
           scale_stat(istat) = scale_rank(iscale_stat(1,1,istat))

           if (iprlvl.ge.3) then
            ! binning
            if (any(iscale_stat(1:4,1,istat)-binmx(1:4).gt.0)) then
             ! resize binning matrix
             do ibin1 = 1, 4
              binmxtmp(ibin1) = 
     &              max(binmx(ibin1),iscale_stat(ibin1,1,istat))
             end do
             allocate(tmp(binmxtmp(1)+1,binmxtmp(2)+1,binmxtmp(3)+1,
     &                    binmxtmp(4)+1))
             tmp = 0
             tmp(1:binmx(1)+1,1:binmx(2)+1,1:binmx(3)+1,1:binmx(4)+1) =
     &          binning(1:binmx(1)+1,1:binmx(2)+1,1:binmx(3)+1,
     &                  1:binmx(4)+1)
             deallocate(binning)
             binmx = binmxtmp
             binning => tmp
             tmp => null()
            end if
            binning(iscale_stat(1,1,istat)+1,iscale_stat(2,1,istat)+1,
     &             iscale_stat(3,1,istat)+1,iscale_stat(4,1,istat)+1) =
     &        binning(iscale_stat(1,1,istat)+1,iscale_stat(2,1,istat)+1,
     &             iscale_stat(3,1,istat)+1,iscale_stat(4,1,istat)+1)+1

           end if
          end if

          ! advance fl_fact_ptr
          do while(fl_fact_ptr%command.ne.command_end_of_formula)
            fl_fact_ptr => fl_fact_ptr%next
          end do
        end if
        if (fl_ptr%command.eq.command_set_target_init.or.
     &      fl_ptr%command.eq.command_set_target_update.or.
     &      fl_ptr%command.eq.command_symmetrise) then
          call new_formula_item(fl_fact_ptr,
     &                          fl_ptr%command,fl_ptr%target)
          if (lustat.gt.0)
     &       call print_form_item(lustat,iitem,fl_fact_ptr,op_info)
          fl_fact_ptr => fl_fact_ptr%next
        end if
        if (.not.associated(fl_ptr%next))
     &       call quit(1,'factorize_new','buggy formula list')
        fl_ptr => fl_ptr%next
      end do
      if (lustat.gt.0)
     &       call print_form_item(lustat,iitem,fl_fact_ptr,op_info)

      nterms = min(iterm,max_stat)
      allocate(ireo_t(nterms),ireo_m(nterms),ireo_s(nterms))

      do iterm = 1, nterms
        ireo_t(iterm) = iterm
      end do
      do iterm = 1, nterms
        ireo_m(iterm) = iterm
      end do
      do iterm = 1, nterms
        ireo_s(iterm) = iterm
      end do

      call idxsortx(time_stat,ireo_t,nterms,-1)
      call idxsortx(mem_stat, ireo_m,nterms,-1)
      call idxsortx(scale_stat,ireo_s,nterms,-1)

      xsum = 0d0
      do iterm = 1, nterms
        xsum = xsum + time_stat(iterm)
      end do

      call write_title(luout,wst_subsection,
     &     'Summary')
      
      write(luout,'(x,"Most expensive contractions: ")') 
      do iterm = 1, min(5,nterms)
        write(luout,'(x," term #",i5,'//
     &            '" - H^",i2," P^",i2," V^",i2," X^",i2'//
     &            '" - flops: ",e10.3,"(",f6.1"%)")')
     &       ireo_t(iterm),iscale_stat(1:4,1,ireo_t(iterm)),
     &       time_stat(iterm),time_stat(iterm)/xsum*100d0
      end do
      write(luout,'(x,"Formally most expensive contractions: ")') 
      do iterm = 1, min(5,nterms)
        write(luout,'(x," term #",i5,'//
     &            '" - H^",i2," P^",i2," V^",i2," X^",i2)')
     &       ireo_s(iterm),iscale_stat(1:4,1,ireo_s(iterm))
      end do
      write(luout,'(x,"Largest intermediates occur in: ")') 
      do iterm = 1, min(5,nterms)
        write(luout,'(x," term #",i5,'//
     &            '" - H^",i2," P^",i2," V^",i2," X^",i2'//
     &            '" - Mb:    ",e10.3)')
     &       ireo_m(iterm),iscale_stat(1:4,2,ireo_m(iterm)),
     &       mem_stat(iterm)/(128d0*1024d0)
      end do

      call file_init(ffstat,trim(label)//'.statistics',ftyp_sq_frm,0)
      call file_open(ffstat)

      write(ffstat%unit,'(x,"Computational cost of contractions: ")') 
      do iterm = 1, nterms
        write(ffstat%unit,'(x," term #",i5,'//
     &            '" - H^",i2," P^",i2," V^",i2," X^",i2'//
     &            '" - flops: ",e10.3,"(",f6.1"%)")')
     &       ireo_t(iterm),iscale_stat(1:4,1,ireo_t(iterm)),
     &       time_stat(iterm),time_stat(iterm)/xsum*100d0
      end do
      write(ffstat%unit,'(x,"Max. size of intermediates: ")') 
      do iterm = 1, nterms
        write(ffstat%unit,'(x," term #",i5,'//
     &            '" - H^",i2," P^",i2," V^",i2," X^",i2'//
     &            '" - Mb:    ",e10.3)')
     &       ireo_m(iterm),iscale_stat(1:4,2,ireo_m(iterm)),
     &       mem_stat(iterm)/(128d0*1024d0)
      end do

      call file_close_keep(ffstat)

      if (iprlvl.ge.3) then
        ! write binning statistics
        write(luout,'(x,55("-"))')
        write(luout,'(x,a)') 'Numbers of terms per formal scaling'
        write(luout,'(x,55("-"))')
        do ibin4 = 1, binmx(4)+1
          do ibin2 = 1, binmx(2)+1
            do ibin1 = 1, binmx(1)+1
              do ibin3 = 1, binmx(3)+1
                if (binning(ibin1,ibin2,ibin3,ibin4).eq.0) cycle
                write(luout,'(x,"H^",i2," P^",i2," V^",i2," X^",i2,'//
     &              '" - number of terms: ",i16)')
     &              ibin1-1,ibin2-1,ibin3-1,ibin4-1,
     &              binning(ibin1,ibin2,ibin3,ibin4)
              end do
            end do
          end do
        end do
        write(luout,'(x,55("-"))')
        deallocate(binning)
      end if

      deallocate(ireo_t,ireo_m,ireo_s,
     &     mem_stat,time_stat,scale_stat,iscale_stat)

      call atim_csw(cpu,sys,wall)
      if (iprlvl.ge.3) 
     &     call prtim(luout,'factorization',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      return
      end


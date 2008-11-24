*----------------------------------------------------------------------*
      subroutine factorize(form_head,op_info,str_info,orb_info,label)
*----------------------------------------------------------------------*
*     read formula from file ffform to linked list 
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
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info
      type(formula_item), intent(in), target ::
     &     form_head
      character*(*), intent(in) ::
     &     label

      integer, parameter ::
     &     max_stat = 1000
      type(formula_item), pointer ::
     &     form_ptr
      integer ::
     &     iterm, nterms
      integer, pointer ::
     &     iscale_stat(:,:,:), ireo_t(:), ireo_m(:), ireo_s(:)
      real(8), pointer ::
     &     time_stat(:), mem_stat(:), scale_stat(:)
      real(8) ::
     &     cpu0, sys0, wall0, cpu, sys, wall, xsum
      type(filinf) ::
     &     ffstat

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

      form_ptr => form_head
      iterm = 0
      do while(form_ptr%command.ne.command_end_of_formula)
        if (form_ptr%command.eq.command_add_contribution) then
          iterm = iterm + 1
          if (iprlvl.ge.10)
     &         write(luout,*) 'factorizing term # ',iterm
          if (iterm.le.max_stat) then
            call form_fact2(form_ptr%contr,
     &         op_info,str_info,orb_info,
     &         iscale_stat(1,1,iterm),time_stat(iterm),mem_stat(iterm))
            scale_stat(iterm) = scale_rank(iscale_stat(1,1,iterm))
          else
            if (iterm.eq.max_stat+1)
     &           call warn('factorize','max_stat exceeded')
            call form_fact2(form_ptr%contr,
     &         op_info,str_info,orb_info,
     &         iscale_stat(1,1,max_stat),
     &           time_stat(max_stat),mem_stat(max_stat))
            scale_stat(max_stat) = scale_rank(iscale_stat(1,1,iterm))
          end if
        end if
        if (.not.associated(form_ptr%next))
     &       call quit(1,'form_opt','buggy formula list')
        form_ptr => form_ptr%next
      end do

      nterms = iterm
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

      deallocate(ireo_t,ireo_m,ireo_s,
     &     mem_stat,time_stat,scale_stat,iscale_stat)

      call atim_csw(cpu,sys,wall)
      if (iprlvl.ge.5) 
     &     call prtim(luout,'factorization',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      return
      end

      real(8) function scale_rank(iscale)

      implicit none

      integer, intent(in) ::
     &     iscale(4)

      scale_rank = (sum(iscale(1:4))       )*100d0*100d0*100d0
     &           + (    iscale(2)+iscale(4))*100d0*100d0
     &           + (    iscale(4)          )*100d0
     &           + (    iscale(1)          )

      end

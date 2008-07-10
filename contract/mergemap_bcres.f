      subroutine mergemap_bcres(mergemap,
     &     ld_map,
     &     svertex,isvtx1,isvtx2,xlines,
     &     ireo,vtx_list_final,nj,nvtx_new,nvtx,nvtx_final)

      implicit none

      include 'stdunit.h'

      integer, intent(in) ::
     &     ld_map, nvtx, nvtx_final, nj, nvtx_new,
     &     svertex(nvtx), isvtx1, isvtx2,
     &     ireo(nvtx),vtx_list_final(nvtx_final)
      integer(8), intent(in) ::
     &     xlines(nvtx_final,nj)
      integer, intent(out) ::
     &     mergemap(ld_map,2,nvtx_final)

      integer ::
     &     ivtx1_ini, ivtx2_ini, ivtx_ini, ivtx_final,
     &     ivtx, jvtx, idx, jdx, mdx
      logical ::
     &     error1, error2

c dbg
c      print *,'nvtx,nvtx_final,nj: ',nvtx,nvtx_final,nj
c      print *,'reo:            ',ireo
c      print *,'vtx_list_final: ',vtx_list_final
c dbg

      mergemap = 0

      error1 = .false.
      error2 = .false.
      ivtx1_ini = 0
      ivtx2_ini = 0
      ! loop over all vertices
      do ivtx = 1, nvtx
        ! does vertex originate from either operator
        ! involved in contraction?
        if (svertex(ivtx).eq.isvtx1) then
          idx = 1
          ivtx1_ini = ivtx1_ini+1
          ivtx_ini = ivtx1_ini
        else if (svertex(ivtx).eq.isvtx2) then
          idx = 2
          ivtx2_ini = ivtx2_ini+1
          ivtx_ini = ivtx2_ini
        else
          cycle
        end if
        ! the vertex to which it contributes finally is
        ! given by the reordering array: 
c dbg
c        print *,'ivtx = ',ivtx
c        print *,'was operator: ',idx
c        print *,'current ivtx/ivtx1/ivtx2_ini:',
c     &                   ivtx_ini,ivtx1_ini,ivtx2_ini 
c dbg
        ivtx_final = ireo(ivtx)
c dbg
c        print *,'ivtx_final = ',ivtx_final
c dbg
        ! last contraction?
        if (nj.eq.nvtx_new) then
          ! look at xlines to find result vertices to which the 
          ! present vertex contributes
          do jdx = 1, nj
            if (xlines(ivtx_final,jdx).eq.0.and.nj.gt.1) cycle
            ! store in mergemap
            mdx = 1
            do while(mdx.le.ld_map.and.mergemap(mdx,idx,jdx).gt.0)
              mdx = mdx+1
            end do
            error2 = mdx.gt.ld_map
            if (.not.error2) mergemap(mdx,idx,jdx) = ivtx_ini    
          end do
        else 
          ! translate into vertex number of intermediate
          jdx = -1
          do jvtx = 1, nvtx_final
            if (vtx_list_final(jvtx).eq.ivtx_final) then
              jdx = jvtx
              exit
            end if
          end do
          error1 = jdx.eq.-1
          if (error1) exit
          ! store in mergemap
          mdx = 1
          do while(mdx.le.ld_map.and.mergemap(mdx,idx,jdx).gt.0)
            mdx = mdx+1
          end do
          error2 = mdx.gt.ld_map
          if (.not.error2) mergemap(mdx,idx,jdx) = ivtx_ini
        end if
      end do

      if (error1) then
        write(luout,*) 'svertex: ',svertex
        write(luout,*) 'isvtx1,isvtx2: ',isvtx1,isvtx2
        write(luout,*) 'reo: ',ireo
        write(luout,*) 'vtx_list_final: ',vtx_list_final(1:nvtx_final)
        write(luout,*) 'not found: ivtx_final = ',ivtx_final
        call quit(1,'mergemap_bcres','something is buggy!')
      end if
      if (error2) then
        write(luout,*) 'ld_map = ',ld_map
        call quit(1,'mergemap_bcres','ld_map too small?')
      end if
c dbg
c      write(luout,*) 'ld_map: ',ld_map
c      write(luout,*) 'mergemap: '
c      do idx = 1, nvtx_final
c        do jdx = 1, 2
c          write(luout,'(x,i4,i4,":",10i4)') 
c     &               idx,jdx,mergemap(1:ld_map,jdx,idx)
c        end do
c      end do
c dbg

      return
      end

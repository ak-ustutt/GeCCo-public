      subroutine mergemap_bcres(mergemap,
     &     ld_map,
     &     svertex,isvtx1,isvtx2,
     &     ireo,vtx_list_final,nvtx,nvtx_final)

      implicit none

      include 'stdunit.h'

      integer, intent(in) ::
     &     ld_map, nvtx, nvtx_final,
     &     svertex(nvtx), isvtx1, isvtx2,
     &     ireo(nvtx),vtx_list_final(nvtx_final)
      integer, intent(out) ::
     &     mergemap(ld_map,2,nvtx_final)

      integer ::
     &     ivtx1_ini, ivtx2_ini, ivtx_ini, ivtx_final,
     &     ivtx, jvtx, idx, jdx, mdx
      logical ::
     &     error1, error2

      mergemap = 0

      error1 = .false.
      error2 = .false.
      ivtx1_ini = 0
      ivtx2_ini = 0
      do ivtx = 1, nvtx
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
        ivtx_final = ireo(ivtx)
        jdx = -1
        do jvtx = 1, nvtx_final
          if (vtx_list_final(jvtx).eq.ivtx_final) then
            jdx = jvtx
            exit
          end if
        end do
        error1 = jdx.eq.-1
        if (error1) exit
        mdx = 1
        do while(mdx.le.ld_map.and.mergemap(mdx,idx,jdx).gt.0)
          mdx = mdx+1
        end do
        error2 = mdx.gt.ld_map
c dbg
c        print '(a,5i5)','ivtx, idx, jdx, mdx, ivtx_ini: ',
c     &           ivtx, idx, jdx, mdx, ivtx_ini
c dbg
        mergemap(mdx,idx,jdx) = ivtx_ini
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
c      write(luout,*) 'mergemap: ',mergemap
c dbg

      return
      end

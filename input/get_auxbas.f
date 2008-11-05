      subroutine get_auxbas(auxcmo,unit)

      implicit none

      integer, intent(in) ::
     &     unit
      real(8), intent(out) ::
     &     auxcmo(*)

      integer ::
     &     nsym, aobas_i, loop, i, j, k, joffs, linind, nbas, idum

      rewind(unit)
      read(unit,'(i5)') nsym

      joffs = 0
      do i = 1, nsym
        read(unit,*)idum,aobas_i,linind,nbas
        do j = 1, linind
          read(unit,*)
          read(unit,*)
          do k = 1, nbas/4
            read(unit,'(4e30.20)') auxcmo(joffs+1:joffs+4)
            joffs=joffs+4
          end do
          k=mod(nbas,4)
          if (k.gt.0) then
            if (k.eq.1) then
              read(unit,'(e30.20)')auxcmo(joffs+1)
            else if (k.eq.2) then
              read(unit,'(2e30.20)')auxcmo(joffs+1:joffs+2)
            else
              read(unit,'(3e30.20)')auxcmo(joffs+1:joffs+3)
            endif              
            joffs=joffs+k
          endif  
        enddo
          
      enddo  
        
      return
      end

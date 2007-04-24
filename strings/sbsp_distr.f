*----------------------------------------------------------------------*
      subroutine sbsp_distr(ipass,
     &     idspc_dis,nsbsp_dis,
     &     nel,nsbsp,mnmxsbsp)
*----------------------------------------------------------------------*
*     set up subspace and IRREP distributions
*     the subspace distribution is equivalent to an occupation
*
*     idspc_dis(1:nel,nsbsp_dis) ---
*       subspace distr:  subspace of components
*
*     ipass = 1   get dimensions
*     ipass = 2   set arrays
*
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'

      integer, intent(in) ::
     &     ipass,
     &     nel, nsbsp,          ! number of el., subspaces
     &     mnmxsbsp(2,nsbsp)    ! min and max number of particles in subsp
                                ! assumed to be adapted to restrictions by
                                ! number of orbitals etc.pp.

      integer, intent(inout) ::
     &     nsbsp_dis

      integer, intent(out) ::
     &     idspc_dis(nel,*)
      
      logical ::
     &     lnext
      integer ::
     &     icurdis(nel), imindis(nel), imaxdis(nel)

      logical, external ::
     &     allow_sbsp_dis, next_rvlex

*----------------------------------------------------------------------*

      nsbsp_dis = 0
      ! start with (1,1,.....,1)
      icurdis(1:nel) = 1
      imindis(1:nel) = 1
      imaxdis(1:nel) = nsbsp
      lnext = .true.
      ! loop over subspace distributions
      dis_loop: do
        ! current distribution conforming with restrictions?
        if (lnext.and.allow_sbsp_dis(icurdis,nel,nsbsp,mnmxsbsp)) then
              ! yes, so we have a new distribution
          nsbsp_dis = nsbsp_dis+1
          if (ipass.eq.2)
     &         idspc_dis(1:nel,nsbsp_dis) = icurdis(1:nel)
        end if
        ! get next distribution
        lnext = next_rvlex(nel,icurdis,imindis,imaxdis) 

        ! last distribution --> exit
        if (.not.lnext) exit dis_loop

      ! end of loop over subspace distributions
      end do dis_loop

      if (nsbsp_dis.eq.0) then
        write(luout,*)
     &       'ERROR in sbsp_distr: cannot create any distribution'
        write(luout,*) ' nel   = ',nel
        write(luout,*) ' nsbsp = ',nsbsp
        write(luout,*) ' mnmxsbsp = ',mnmxsbsp(1:2,nsbsp)
        stop 'ERROR in sbsp_distr'
      end if

      return
      end

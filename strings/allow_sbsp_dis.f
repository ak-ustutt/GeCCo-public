*----------------------------------------------------------------------*
      logical function allow_sbsp_dis(idistr,nel,nsbsp,mnmxsbsp)
*----------------------------------------------------------------------*
*     check whether distribution of subspaces conforms with minimum
*     and maximum accumulated occupations given on mnmxsbsp
*
*       mnmxsbsp(1,ispc) :   minimum number of orbitals to be
*                            occupied after spaces 1-ispc are occ.
*       mnmxsbsp(2,ispc) :   maximum number of orbitals to be
*                            occupied after spaces 1-ispc are occ.
*----------------------------------------------------------------------*

      implicit none

      integer, intent(in) ::
     &     nel, nsbsp,
     &     idistr(nel), mnmxsbsp(2,nsbsp)

      logical ::
     &     ok
      integer ::
     &     iocc(nsbsp), iel, ispc, iocct


      ! set up occupation
      iocc(1:nsbsp) = 0
      do iel = 1, nel
        iocc(idistr(iel)) = iocc(idistr(iel)) + 1
      end do

      ! get accumulated occupation after space ispc and 
      ! compare to restrictions
      iocct = 0
      ok = .true.
      do ispc = 1, nsbsp
        iocct = iocct + iocc(ispc)
        ok = ok.and.iocct.ge.mnmxsbsp(1,ispc)
     &         .and.iocct.le.mnmxsbsp(2,ispc)
      end do

      allow_sbsp_dis = ok

      return
      end

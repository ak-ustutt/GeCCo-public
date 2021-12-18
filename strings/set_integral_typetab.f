*----------------------------------------------------------------------*
      subroutine set_integral_typetab(ntypes,nelmax,typetab)
*----------------------------------------------------------------------*
*     hard-wired for nelmax = 4
*     ntypes = 3 :  12-symmetry + CA-symmetry + individual CA-symmetry
*                 (e.g. 1/r12 and f(r12) integrals 
*                 -> 3 distinct types per index quadruple
*     ntypes = 6 :  no individual CA-symmetry (e.g. [T12,r12] integrals)
*                 -> 6 distinct types per index quadruple
*     ntypes = 12:  no CA-symmetry at all (e.g. K-transformed r12 integrals)
*                 -> 12 distinct types per index quadruple
*     typetab contains for all the 24 permuations the corresponding
*     integral type and a sign
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'stdunit.h'
      include 'def_orbinf.h'
      include 'multd2h.h'

      integer, parameter ::
     &     p12(4)   = (/2,1,3,4/),
     &     p34(4)   = (/1,2,4,3/),
     &     p1234(4) = (/2,1,4,3/),
     &     p1324(4) = (/3,4,1,2/)

      integer, intent(in) ::
     &     ntypes, nelmax
      integer, intent(out) ::
     &     typetab(24)
      
      integer ::
     &     rank(4), perm(4,24), perm0(4), perm1(4), perm2(4),
     &     idx, jdx, itype, fac
      integer, external ::
     &     rank_ivec

      if (nelmax.ne.4)
     &     call quit(1,'set_integral_typetab',
     &     'hard-wired for nelmax==4')

      ! set up all permutations
      perm(1:4,1) = (/1,2,3,4/)
      do idx = 2, 24
        perm(1:4,idx) = perm(1:4,idx-1)
        call next_perm(perm(1,idx),4)
      end do

      typetab = 0

      itype = 1
      do
        ! look for next free permutaion
        jdx = -1
        do idx = 1, 24
          if (typetab(idx).eq.0) then
            jdx = idx
            exit
          end if
        end do
        if (jdx.eq.-1) exit

        if (ntypes.eq.3) then
        ! 1
        typetab(jdx) = itype
        perm0(1:4) = perm(1:4,jdx)
        if (ntest.ge.100)
     &       write(lulog,*) '1: ',perm0
        
        ! 2
        call perm_mult(perm1,perm0,p12,4)
        jdx = rank_ivec(rank,perm1,4)+1
        typetab(jdx) = itype
        if (ntest.ge.100)
     &       write(lulog,*) '2: ',perm1

        ! 3
        call perm_mult(perm1,perm0,p34,4)
        jdx = rank_ivec(rank,perm1,4)+1
        typetab(jdx) = itype
        if (ntest.ge.100)
     &       write(lulog,*) '3: ',perm1

        ! 4
        call perm_mult(perm1,perm0,p1324,4)
        jdx = rank_ivec(rank,perm1,4)+1
        typetab(jdx) = itype
        if (ntest.ge.100)
     &       write(lulog,*) '4: ',perm1

        ! 5
        call perm_mult(perm1,perm0,p12,4)
        call perm_mult(perm2,perm1,p34,4)
        jdx = rank_ivec(rank,perm2,4)+1
        typetab(jdx) = itype
        if (ntest.ge.100)
     &       write(lulog,*) '5: ',perm2

        ! 6
        call perm_mult(perm1,perm0,p12,4)
        call perm_mult(perm2,perm1,p1324,4)
        jdx = rank_ivec(rank,perm2,4)+1
        typetab(jdx) = itype
        if (ntest.ge.100)
     &       write(lulog,*) '6: ',perm2

        ! 7
        call perm_mult(perm1,perm0,p34,4)
        call perm_mult(perm2,perm1,p1324,4)
        jdx = rank_ivec(rank,perm2,4)+1
        typetab(jdx) = itype
        if (ntest.ge.100)
     &       write(lulog,*) '7: ',perm2

        ! 8
        call perm_mult(perm1,perm0,p12,4)
        call perm_mult(perm2,perm1,p34,4)
        call perm_mult(perm1,perm2,p1324,4)
        jdx = rank_ivec(rank,perm1,4)+1
        typetab(jdx) = itype
        if (ntest.ge.100)
     &       write(lulog,*) '8: ',perm1

        else if (ntypes.eq.6) then

        ! 1          
        perm0(1:4) = perm(1:4,jdx)
        if ((perm0(1)-1)*4+perm0(3).gt.perm0(2)*4+perm0(4)) fac = -1
        fac = 1
        typetab(jdx) = fac*itype
        if (ntest.ge.100)
     &       write(lulog,*) '1: ',perm0,fac,jdx
        
        ! 2
        call perm_mult(perm1,perm0,p1234,4)
        jdx = rank_ivec(rank,perm1,4)+1
        fac = 1
        if (rank(1)*4+rank(3).gt.rank(2)*4+rank(4)) fac = -1
        typetab(jdx) = fac*itype
        if (ntest.ge.100)
     &       write(lulog,*) '2: ',perm1,fac,jdx

        ! 3
        call perm_mult(perm1,perm0,p1324,4)
        jdx = rank_ivec(rank,perm1,4)+1
        fac = 1
        if (rank(1)*4+rank(3).gt.rank(2)*4+rank(4)) fac = -1
        typetab(jdx) = fac*itype
        if (ntest.ge.100)
     &       write(lulog,*) '3: ',perm1,fac,jdx

        ! 4
        call perm_mult(perm1,perm0,p1234,4)
        call perm_mult(perm2,perm1,p1324,4)
        jdx = rank_ivec(rank,perm2,4)+1
        fac = 1
        if (rank(1)*4+rank(3).gt.rank(2)*4+rank(4)) fac = -1
        typetab(jdx) = fac*itype
        if (ntest.ge.100)
     &       write(lulog,*) '4: ',perm2,fac,jdx

        else if (ntypes.eq.12) then

        ! 1          
        perm0(1:4) = perm(1:4,jdx)
c        if ((perm0(1)-1)*4+perm0(3).gt.perm0(2)*4+perm0(4)) fac = -1
        fac = 1
        typetab(jdx) = fac*itype
        if (ntest.ge.100)
     &       write(lulog,*) '1: ',perm0,fac,jdx
        
        ! 2
        call perm_mult(perm1,perm0,p1324,4)
        jdx = rank_ivec(rank,perm1,4)+1
        fac = 1
c        if (rank(1)*4+rank(3).gt.rank(2)*4+rank(4)) fac = -1
        typetab(jdx) = fac*itype
        if (ntest.ge.100)
     &       write(lulog,*) '2: ',perm1,fac,jdx

        end if

        itype = itype+1

        if (itype.gt.12) exit
      end do

      if (ntest.ge.100) then
        write(lulog,*) 'typetab: '
        write(lulog,'(5x,5i4,x,5i4)') typetab
      end if

      
      return
      end

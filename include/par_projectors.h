* new definition of projectors: count H and V as occupied (O)
* needed for MR generalization
      integer, parameter ::
     &     nproj = 11,
*     blocks of projector
     &     nblk(nproj) = (/1,9,2,8,2,4,1,6,3,2,1/),
*     offset of projector
     &     ioff(nproj) = (/0,1,10,12,20,22,26,27,33,36,38/)!,
c PGF90 only accepts this
*     the projectors:
      integer ::
     &     occ_prj(ngastp,39)
      data occ_prj /
*     ansatz 1, Q
     &     0,0,0,2, !  1
*     ansatz 1, P
     &     2,0,0,0, !  2
     &     1,0,1,0, !  3 V
     &     0,0,2,0, !  4 V
     &     1,1,0,0, !  5
     &     0,1,1,0, !  6 V
     &     1,0,0,1, !  7
     &     0,0,1,1, !  8 V
     &     0,2,0,0, !  9
     &     0,1,0,1, ! 10
*     ansatz 2, Q
     &     0,1,0,1, ! 11
     &     0,0,0,2, ! 12
*     ansatz 2, P
     &     2,0,0,0, ! 13
     &     1,0,1,0, ! 14 V
     &     0,0,2,0, ! 15 V
     &     1,1,0,0, ! 16
     &     0,1,1,0, ! 17 V
     &     1,0,0,1, ! 18
     &     0,0,1,1, ! 19 V
     &     0,2,0,0, ! 20
*     O1+O2-2O1O2-V1O2-V2O1  Z for ansatz 1 
     &     1,0,0,1, ! 21 
     &     0,0,1,1, ! 22 V
*     O1+O2-2O1O2            Z for ansatz 2
     &     1,1,0,0, ! 23
     &     0,1,1,0, ! 24 V
     &     1,0,0,1, ! 25 
     &     0,0,1,1, ! 26 V
*     P1P2                   
     &     0,2,0,0, ! 27
*     O1O2+O1P2+P1O2+P1P2
     &     2,0,0,0, ! 28
     &     1,0,1,0, ! 29 V
     &     0,0,2,0, ! 30 V
     &     1,1,0,0, ! 31
     &     0,1,1,0, ! 32 V
     &     0,2,0,0, ! 33
*     O1X2+X1O2+P1X2+X2P1  for ansatz 1
     &     1,0,0,1, ! 34
     &     0,0,1,1, ! 35 V
     &     0,1,0,1, ! 36
*     O1X2+X1O2            for ansatz 3
     &     1,0,0,1, ! 37
     &     0,0,1,1, ! 38 V
*     V1X2+X1V2            for ansatz 3
     &     0,1,0,1/ ! 39    

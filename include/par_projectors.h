      integer, parameter ::
     &     nproj = 11,
*     blocks of projector
     &     nblk(nproj) = (/1,5,2,4,1,2,1,3,2,1,1/),
*     offset of projector
     &     ioff(nproj) = (/0,1,6,8,12,13,15,16,19,21,22/)!,
c*     the projectors:
c     &     occ_prj(ngastp,19) = (/
c*     ansatz 1, Q
c     &     (/0,0,0,2/),
c*     ansatz 1, P
c     &     (/2,0,0,0/),
c     &     (/1,1,0,0/),
c     &     (/1,0,0,1/),
c     &     (/0,2,0,0/),
c     &     (/0,1,0,1/),
c*     ansatz 2, Q
c     &     (/0,1,0,1/),
c     &     (/0,0,0,2/),
c*     ansatz 2, P
c     &     (/2,0,0,0/),
c     &     (/1,1,0,0/),
c     &     (/1,0,0,1/),
c     &     (/0,2,0,0/),
c*     O1+O2-2O1O2-V1O2-V2O1  Z for ansatz 1 
c     &     (/1,0,0,1/),
c*     O1+O2-2O1O2            Z for ansatz 2
c     &     (/1,1,0,0/),
c     &     (/1,0,0,1/),
c*     P1P2
c     &     (/0,2,0,0/)/)
c PGF90 only accepts this:
*     the projectors:
      integer ::
     &     occ_prj(ngastp,23)
      data occ_prj /
*     ansatz 1, Q
     &     0,0,0,2,
*     ansatz 1, P
     &     2,0,0,0,
     &     1,1,0,0,
     &     1,0,0,1,
     &     0,2,0,0,
     &     0,1,0,1,
*     ansatz 2, Q
     &     0,1,0,1,
     &     0,0,0,2,
*     ansatz 2, P
     &     2,0,0,0,
     &     1,1,0,0,
     &     1,0,0,1,
     &     0,2,0,0,
*     O1+O2-2O1O2-V1O2-V2O1  Z for ansatz 1 
     &     1,0,0,1,
*     O1+O2-2O1O2            Z for ansatz 2
     &     1,1,0,0,
     &     1,0,0,1,
*     P1P2                   
     &     0,2,0,0,
*     O1O2+O1P2+P1O2+P1P2
     &     2,0,0,0,
     &     1,1,0,0,
     &     0,2,0,0,
*     O1X2+X1O2+P1X2+X2P1  for ansatz 1
     &     1,0,0,1,
     &     0,1,0,1,
*     O1X2+X1O2            for ansatz 3
     &     1,0,0,1,
*     V1X2+X1V2            for ansatz 3
     &     0,1,0,1/     

      integer, parameter ::
     &     nproj = 6,
*     blocks of projector
     &     nblk(nproj) = (/1,5,2,4,1,2/),
*     offset of projector
     &     ioff(nproj) = (/0,1,6,8,12,13/),
*     the projectors:
     &     occ_prj(ngastp,15) = (/
*     ansatz 1, Q
     &     (/0,0,0,2/),
*     ansatz 1, P
     &     (/2,0,0,0/),
     &     (/1,1,0,0/),
     &     (/1,0,0,1/),
     &     (/0,2,0,0/),
     &     (/0,1,0,1/),
*     ansatz 2, Q
     &     (/0,1,0,1/),
     &     (/0,0,0,2/),
*     ansatz 2, P
     &     (/2,0,0,0/),
     &     (/1,1,0,0/),
     &     (/1,0,0,1/),
     &     (/0,2,0,0/),
*     O1+O2-2O1O2-V1O2-V2O1  Z for ansatz 1 
     &     (/1,0,0,1/),
*     O1+O2-2O1O2            Z for ansatz 2
     &     (/1,1,0,0/),
     &     (/1,0,0,1/)
     &     /)

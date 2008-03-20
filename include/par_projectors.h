      integer, parameter ::
     &     nproj = 4,
*     blocks of projector
     &     nblk(nproj) = (/1,5,2,4/),
*     offset of projector
     &     ioff(nproj) = (/0,1,6,8/),
*     the projectors:
     &     occ_prj(ngastp,12) = (/
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
     &     (/0,2,0,0/)
     &     /)

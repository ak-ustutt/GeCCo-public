      
      type vector_t
      integer::
     &     len=0
      type(me_list_array),pointer::
     &     me_lists(:) => null()
      end type

      type vector_space_t
      integer::
     &     maxvec,              !maximum dimension of the subspace
     &     nvec,             !current dimension of the subspace
     &     nlists,              !number of lists a vector consists of
     &     ivec                 ! currently active vector == current_record for all associated files (hopefully)
      type(me_list_array),dimension(:),pointer::
     &     me_lists=>null()            ! saves information about the shapeof the vectors (do not have to be associated to the files)
      type(file_array),dimension(:),pointer::
     &     vectors=>null()              ! Note: a vector may consist of multiple lists. the vector number is the requested record.
                                        !the index in vectors is only the actual list requested as part of this vector
      end type

!       { (c)_1  (c)_2  ...}       <- vectors(1)
!     U={ (t1)_1 (t1)_2 ...}       <- vectors(2)
!       { (t2)_1 (t2)_2 ...}       <- vectors(3)




      type davidson_subspace_t
      integer::
     &     nmaxsub,             !maximum number of vectors 
     &     ncursub,             ! current number of vectors
     &     icursub,              ! last updated vectors (if new vectors are added after nmaxsub is reached, vectors are overridden)
     &     lcursub       !
      real(8),dimension(:),pointer::
     &     vMv_mat=>null()           ! matrix of all vMv products it is a maxsub x maxsub matrix  with ncursub xncursub  elements !=0 
! (vMv_mat)_{i,j} = (v)_i*(Mv)_j

      type(vector_space_t)::
     &     vspace,        ! vector space
     &     Mvspace        ! matrix-vector product space

      end type

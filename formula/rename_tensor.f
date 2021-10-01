*----------------------------------------------------------------------*
      pure function rename_tensor(string, rank, nops, itf_names)
*----------------------------------------------------------------------*
!     Rename tensor according to taste
!     This should be expanded to rename all tensors so they correspond
!     with the ITF algo file
*----------------------------------------------------------------------*

      implicit none
      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'

      character(len=MAXLEN_BC_LABEL), intent(in) ::
     &    string
      integer, intent(in) ::
     &    rank,          ! Rank of tensor
     &    nops(ngastp)   ! Number of creation/annhilation ops per excitation space
      type(tensor_names), intent(in) ::
     &     itf_names            ! contains renaming information
      character(len=MAXLEN_BC_LABEL) ::
     &     rename_tensor

      logical ::
     &     found
      integer ::
     &     ii

      found = .false.
      do ii = 1, itf_names%nrename
        if (trim(string)==trim(itf_names%gc_itf_rename(1,ii))) then
          found = .true.
          select case(itf_names%rename_type(ii))
          case(RENAME_BASIC)
            rename_tensor = trim(itf_names%gc_itf_rename(2,ii))
          case(RENAME_HAMILTONIAN)
            if (rank==2) then
              rename_tensor='f'
            else
              if (nops(2)==3) then
                rename_tensor='J' ! Currently use J:eeec for 3 external integrals
              else
                rename_tensor='K'
              end if
            end if
          case(RENAME_ADD_RANK)
            if (rank>18) then
              !write(lulog,*) 'rank = ',rank
              rename_tensor='ERROR_unexp_lg_rank'
            end if
            write(rename_tensor,'(a,i1)')
     &           trim(itf_names%gc_itf_rename(2,ii)),rank/2
          case default
            !write(lulog,*) 'case ',itf_names%rename_type
            rename_tensor='ERROR_unknown_case'
          end select
          exit
        end if
      end do

      if (.not.found) rename_tensor = trim(string)
      
      if (rename_tensor(1:1)=='_') rename_tensor(1:1)=''

      end function

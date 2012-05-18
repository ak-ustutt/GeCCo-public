      subroutine parse_hpvx_descr(occ_list,nlist,string,maxlist)

      implicit none

      integer, intent(in) ::
     &     maxlist
      character(len=*), intent(in) ::
     &     string
      integer, intent(out) ::
     &     nlist, occ_list(4,maxlist)
      
      logical ::
     &     brackets, allowed(4)
      integer ::
     &     len_string, ipos, nset, iset, hpvx, nlist_uniq
      

      len_string = len_trim(string)

      nlist = 0
      ipos = 1
      brackets = .false.
      nset = 0
      allowed(1:4) = .false.
      ! loop over string
      do while(ipos.le.len_string)
        ! single character or '[' ?
        select case (string(ipos:ipos))
          case ('H','h') 
            allowed(1) = .true.
            nset = nset + 1
          case ('P','p') 
            allowed(2) = .true.
            nset = nset + 1
          case ('V','v') 
            allowed(3) = .true.
            nset = nset + 1
          case ('X','x') 
            allowed(4) = .true.
            nset = nset + 1
          case ('[') 
            brackets = .true.
            nset = 0
          case (']')
            if (.not.brackets)
     &           call quit(1,'parse_hpvx_descr',
     &           'unexpected "]" in "'//trim(string)//'"')
            brackets = .false.
          case default
            call quit(1,'parse_hpvx_descr',
     &           'syntax error in "'//trim(string)//'"')
        end select
        if (.not.brackets) then
          if (nlist.eq.0) then
            nlist = 1
            occ_list(1:4,1) = 0
          end if
          if (nset*nlist.gt.maxlist)
     &         call quit(1,'parse_hpvx_descr',
     &           'maxlist is too small!')
          do iset = 2, nset
            occ_list(1:4,(iset-1)*nlist+1:iset*nlist)
     &           = occ_list(1:4,1:nlist)
          end do
          iset = 0
          do hpvx = 1, 4
            if (.not.allowed(hpvx)) cycle
            iset = iset + 1
            occ_list(hpvx,(iset-1)*nlist+1:iset*nlist) =
     &           occ_list(hpvx,(iset-1)*nlist+1:iset*nlist) + 1
          end do
          if (iset.ne.nset)
     &         call quit(1,'parse_hpvx_descr','inconsistency!')
          nlist = nlist*nset
          nset = 0
          allowed(1:4) = .false.
        end if
        ipos = ipos + 1
      end do

      if (brackets)
     &     call quit(1,'parse_hpvx_descr',
     &           'missing "]" in "'//trim(string)//'"')

      ! at least: empty occupation
      if (nlist.eq.0) then
        occ_list(1:4,1) = 0
        nlist = 1
      end if

      ! remove duplicate entries
      call unique_ntupel(occ_list,4,nlist,nlist_uniq)
      nlist = nlist_uniq

      end

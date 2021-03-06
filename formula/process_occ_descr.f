*--------------------------------------------------------------------------*
!>    parse a general occupation descriptor of the syntax
!>
!>    <OP_DESCR>|<OP_DESCR>|... 
!>    <OP_DESCR>: (<HPVX_DESCR>,<HPVX_DSCR>);(<HPVX_DESCR>,<HPVX_DESCR>)...
!>    <HPVX_DSCR>: [HPVX]* 
!>    separators: ' ,;|'
!>    quotes: '"'
!>    @param[out] occ_list 
!>    @param[out] nlist used length of occ_list
!>    @param[in] dstring Descriptor string
!>    @param[in] njoined number of vertices
!>    @maxlist[in] allocated length of occ_list in operator descriptors
*--------------------------------------------------------------------------*
      subroutine process_occ_descr(occ_list,nlist,
     &                             dstring,njoined,maxlist)
*--------------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'def_word_list.h'

      integer, parameter ::
     &     ntest = 00
      character(len=18) ::
     &     i_am = 'process_occ_descr'
      

      integer, intent(in) ::
     &     maxlist, njoined
      character(len=*), intent(in) ::
     &     dstring
      integer, intent(out) ::
     &     nlist, occ_list(4,2,njoined,maxlist)

      type(word_list) ::
     &     wlist
      type(word_list_entry) ::
     &     wl_pnt

      character(len=maxlen_word) ::
     &     cur_entry
      character ::
     &     sep
      logical ::
     &     next
      integer ::
     &     icnt, ilist, ijoin, kjoin, nlist_uniq
      integer ::
     &     nl_c(njoined), nl_a(njoined), 
     &     il_a(njoined), il_c(njoined), n_one(njoined),
     &     occ_list_c(4,maxlist,njoined),
     &     occ_list_a(4,maxlist,njoined)

      logical, external ::
     &     advance_word_list_entry, next_dist2

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,i_am)
        write(lulog,*) 'string: "',trim(dstring),'"'
        write(lulog,*) ' njoined, maxlist: ',njoined,maxlist
      end if

      call init_word_list(wlist)
      call lex_line(wlist,dstring,0,' ,;|',' ','"')

      if (ntest.ge.100) then
        write(lulog,*) 'lex''ed string:'
        call print_word_list(lulog,wlist)
      end if

      nl_c = 1
      nl_a = 1
      occ_list_c = 0
      occ_list_a = 0

      nlist = 0

      ! move internal pointer in wlist to beginning
      call reset_word_list_pointer(wlist)

      icnt = 0
      ijoin = 0
      main_lp: do 
        call get_word_list_entry(cur_entry,sep,wlist)
        next = advance_word_list_entry(wlist,' ')
        icnt = icnt+1

        if (mod(icnt,2).eq.1) ijoin = ijoin+1

        if (ntest.ge.100) then
          write(lulog,*) 'ijoin,icnt: ',ijoin,icnt
          write(lulog,*) 'cur_entry: "',trim(cur_entry),'" sep=',sep
        end if

        if (mod(icnt,2).eq.1) then
          call parse_hpvx_descr(occ_list_c(1:,1:,ijoin),
     &                          nl_c(ijoin),cur_entry,maxlist)
        else
          call parse_hpvx_descr(occ_list_a(1:,1:,ijoin),
     &                          nl_a(ijoin),cur_entry,maxlist)
        end if

        if (ijoin.eq.njoined.and.mod(icnt,2).eq.0) then

          if (ntest.ge.100) then
            write(lulog,*) 'C: ',nl_c
            write(lulog,*) 'A: ',nl_a
          end if

          n_one = 1

          il_c = 1
          il_c_lp: do 
            il_a = 1
            il_a_lp: do
              nlist = nlist + 1
              if (nlist.gt.maxlist) exit main_lp
              do kjoin = 1, njoined

           occ_list(1:4,1,kjoin,nlist)=occ_list_c(1:4,il_c(kjoin),kjoin)
           occ_list(1:4,2,kjoin,nlist)=occ_list_a(1:4,il_a(kjoin),kjoin)

              end do
              if (.not.next_dist2(il_a,njoined,n_one,nl_a,1)) 
     &          exit il_a_lp
            end do il_a_lp
            if (.not.next_dist2(il_c,njoined,n_one,nl_c,1)) 
     &        exit il_c_lp
          end do il_c_lp

          ! missing here:
          ! merge this list with previous entries
          call unique_ntupel(occ_list,4*2*njoined,nlist,nlist_uniq)
          nlist = nlist_uniq

          ! reset for next round
          nl_c = 1
          nl_a = 1
          occ_list_c = 0
          occ_list_a = 0
        end if

        if (.not.next) then
          if (ijoin.eq.njoined) then
            exit
          else
            call quit(0,i_am,
     &    'string ended before definition was complete: "'
     &    //trim(dstring)//'"')
          end if 
        else
          if (ijoin.eq.njoined.and.mod(icnt,2).eq.0) ijoin = 0
        end if
      end do main_lp

      if (ntest.ge.100)
     &     write(lulog,*) 'nlist = ',nlist

      if (nlist.gt.maxlist)
     &     call quit(1,i_am,'maxlen too small')

      if (ntest.ge.100) then
        write(lulog,*) 'final occupation list'
        do ilist = 1, nlist
          call wrt_occ_n(lulog,occ_list(1:,1:,1:,ilist),njoined)
        end do
      end if
 
      return

      end subroutine


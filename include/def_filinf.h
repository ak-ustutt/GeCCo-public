*------------------------------------------------------------------------*
*	info about file
*------------------------------------------------------------------------*
      integer, parameter ::
     &     ftyp_da_unf = 1,
     &     ftyp_sq_unf = 2,
     &     ftyp_sq_frm = 3
*------------------------------------------------------------------------*
      integer, parameter ::
     &     maxfilnam = 512 
      type filinf
         character*(maxfilnam) ::
     &     name
         integer ::
     &     unit 
         integer ::
     &     type, reclen
         logical ::
     &     buffered
         integer, pointer ::
     &     idxrec(:)
         integer, pointer ::
     &     incore(:)
         real(8), pointer ::
     &     buffer(:)
      end type filinf
*------------------------------------------------------------------------*
*	name, unit: self-explaining, unit == -1 means not assigned
*       type: 1 direct access
*             2 sequential unformatted
*             3 sequential formatted
*       reclen: record length for direct access
*       buffered: set to .false. by default, if set to .true., the
*             behaviour is "user-defined", i.e. some of the below
*             array should be set to the appropriate values
*	idxrec: an index array for direct access files ("user defined")
*       incore: an array for indicating that certain parts of the
*               file are incore ("user defined", e.g. blocks of operator
*               or records of the file ...)
*       buffer: a buffer for incore parts of the file contents
*------------------------------------------------------------------------*
*	a note on buffering for "operators":
*	incore(iocc_cls) holds info whether the occupation iocc_cls
*	is held incore (in the buffer) (if > 0)
*	as the length increases with iocc_cls, we currently will hold 
*       the first few occupations (say 1 to 4 for the Hamilton) incore, 
*       so the addressing array on the operator structure can be used 
*       for "buffer" as well
*------------------------------------------------------------------------*

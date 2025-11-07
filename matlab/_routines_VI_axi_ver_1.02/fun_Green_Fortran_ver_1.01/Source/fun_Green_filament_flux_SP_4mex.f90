#include "fintrf.h"      

subroutine mexFunction(nlhs, plhs, nrhs, prhs)

implicit none

! mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs

! stuff for mexing
      mwPointer mxCreateDoubleMatrix
      mwPointer mxGetPr
	  mwPointer npt_source_pr, r_source_pr, z_source_pr, I_source_pr
      mwPointer npt_point_pr, r_point_pr, z_point_pr
      mwPointer OPT_PARALLEL_pr
      mwPointer n_thread_pr
      mwPointer vec_flux_pr

      mwPointer m, n
      mwSize mo,no,size
      integer flag,i,j
      mwPointer mxGetM, mxGetN
      integer mxIsNumeric 


! fortran subroutine arguments source
    real*8, allocatable, dimension(:) :: r_source, z_source, I_source
    real*8 npt_sourcer
    integer npt_source

! fortran subroutine arguments target points
    real*8, allocatable, dimension(:) :: r_point, z_point
    real*8 npt_pointr
    integer npt_point

! fortran subroutine arguments OPT_PARALLEL
    real*8 OPT_PARALLELr
    integer OPT_PARALLEL

! fortran subroutine arguments n_thread
    real*8 n_threadr
    integer n_thread

! fortran output
    real*8, allocatable, dimension(:) :: vec_flux


!                    i1         i2       i3       i4       i5        i6      i7      i8           i9       o1
!call fun_Green_flux(npt_source,r_source,z_source,I_source,npt_point,r_point,z_point,OPT_PARALLEL,n_thread,vec_flux)


! check for proper number of arguments.
      if (nrhs .ne. 9) then
        call mexErrMsgIdAndTxt ('MATLAB:fun_Green_flux:nInput', &
                                '8 input arguments required.')
      elseif (nlhs .ne. 1) then
        call mexErrMsgIdAndTxt ('MATLAB:fun_mutua:nOutput', &
                                '1 output argument required.')
      endif

!    Check to see inputs are numeric.
      do i=1,9
        if (mxIsNumeric(prhs(i)) .ne. 1) then
            call mexErrMsgIdAndTxt ('MATLAB:fun_Green_flux:nInput', &
                                    'Inputs must be numeric')
        endif
	  enddo



!     Check that input #1 is scalar and fetch it
      m = mxGetM(prhs(1))
      n = mxGetN(prhs(1))
      if(m .ne. 1 .or. n .ne. 1) then
        call mexErrMsgIdAndTxt ('MATLAB:fun_Green_flux:nInput', &
                                'Input 1 must be scalar')
      endif	  
      size = m*n
      npt_source_pr = mxGetPr(prhs(1))
      call mxCopyPtrToReal8(npt_source_pr, npt_sourcer, size)
      npt_source=int(npt_sourcer)

!     Check that input #2 is a npt_sourcex1 array and fetch it
      m = mxGetM(prhs(2))
      n = mxGetN(prhs(2))
      if(m .ne. npt_source .or. n .ne. 1) then
        call mexErrMsgIdAndTxt ('MATLAB:fun_Green_flux:nInput', &
                                'Input 2 must be npt_sourcex1')
      endif
      size = m*n
      r_source_pr = mxGetPr(prhs(2))
      allocate(r_source(npt_source))
      call mxCopyPtrToReal8(r_source_pr, r_source, size)

!     Check that input #3 is a npt_sourcex1 array and fetch it
      m = mxGetM(prhs(3))
      n = mxGetN(prhs(3))
      if(m .ne. npt_source .or. n .ne. 1) then
        call mexErrMsgIdAndTxt ('MATLAB:fun_Green_flux:nInput', &
                                'Input 3 must be npt_sourcex1')
      endif
      size = m*n
      z_source_pr = mxGetPr(prhs(3))
      allocate(z_source(npt_source))
      call mxCopyPtrToReal8(z_source_pr, z_source, size)

!     Check that input #4 is a npt_sourcex1 array and fetch it
      m = mxGetM(prhs(4))
      n = mxGetN(prhs(4))
      if(m .ne. npt_source .or. n .ne. 1) then
        call mexErrMsgIdAndTxt ('MATLAB:fun_Green_flux:nInput', &
                                'Input 4 must be npt_sourcex1')
      endif
      size = m*n
      I_source_pr = mxGetPr(prhs(4))
      allocate(I_source(npt_source))
      call mxCopyPtrToReal8(I_source_pr, I_source, size)

!     Check that input #5 is scalar and fetch it
      m = mxGetM(prhs(5))
      n = mxGetN(prhs(5))
      if(m .ne. 1 .or. n .ne. 1) then
        call mexErrMsgIdAndTxt ('MATLAB:fun_Green_flux:nInput', &
                                'Input 5 must be scalar')
      endif	  
      size = m*n
      npt_point_pr = mxGetPr(prhs(5))
      call mxCopyPtrToReal8(npt_point_pr, npt_pointr, size)
      npt_point=int(npt_pointr)

!     Check that input #6 is a npt_pointx1 array and fetch it
      m = mxGetM(prhs(6))
      n = mxGetN(prhs(6))
      if(m .ne. npt_point .or. n .ne. 1) then
        call mexErrMsgIdAndTxt ('MATLAB:fun_Green_flux:nInput', &
                                'Input 6 must be npt_pointx1')
      endif
      size = m*n
      r_point_pr = mxGetPr(prhs(6))
      allocate(r_point(npt_point))
      call mxCopyPtrToReal8(r_point_pr, r_point, size)

!     Check that input #7 is a npt_pointx1 array and fetch it
      m = mxGetM(prhs(7))
      n = mxGetN(prhs(7))
      if(m .ne. npt_point .or. n .ne. 1) then
        call mexErrMsgIdAndTxt ('MATLAB:fun_Green_flux:nInput', &
                                'Input 7 must be npt_pointx1')
      endif
      size = m*n
      z_point_pr = mxGetPr(prhs(7))
      allocate(z_point(npt_point))
      call mxCopyPtrToReal8(z_point_pr, z_point, size)

!     Check that input #8 is scalar and fetch it
      m = mxGetM(prhs(8))
      n = mxGetN(prhs(8))
      if(m .ne. 1 .or. n .ne. 1) then
        call mexErrMsgIdAndTxt ('MATLAB:fun_Green_flux:nInput', &
                                'Input 8 must be scalar')
      endif	  
      size = m*n
      OPT_PARALLEL_pr = mxGetPr(prhs(8))
      call mxCopyPtrToReal8(OPT_PARALLEL_pr, OPT_PARALLELr, size)
      OPT_PARALLEL=int(OPT_PARALLELr)

!     Check that input #9 is scalar and fetch it
      m = mxGetM(prhs(9))
      n = mxGetN(prhs(9))
      if(m .ne. 1 .or. n .ne. 1) then
        call mexErrMsgIdAndTxt ('MATLAB:fun_Green_flux:nInput', &
                                'Input 9 must be scalar')
      endif	  
      size = m*n
      n_thread_pr = mxGetPr(prhs(9))
      call mxCopyPtrToReal8(n_thread_pr, n_threadr, size)
      n_thread=int(n_threadr)




! call the computational subroutine
        allocate(vec_flux(npt_point))
        call fun_Green_filament_flux_SP_4mex(npt_source,r_source,z_source,I_source, &
                            npt_point,r_point,z_point,OPT_PARALLEL,n_thread,vec_flux)

! Create an array for the #1 return argument
      mo=npt_point
      no=1
      plhs(1) = mxCreateDoubleMatrix(mo, no, 0)
! Load the output into a MATLAB array.
      vec_flux_pr = mxGetPr(plhs(1))
      size=mo*no
      call mxCopyReal8ToPtr(vec_flux, vec_flux_pr, size)

!      deallocate(r1,z1,t1,t1r,r2,z2,t2,t2r)
	  deallocate(vec_flux)

    return
    end
!------------------------------------------------------------------------------------------------------------------


!------------------------------------------------------------------------------------------------------------------
    subroutine fun_Green_filament_flux_SP_4mex(npt_source,r_source,z_source,I_source, &
                              npt_point,r_point,z_point,OPT_PARALLEL,n_thread,vec_flux)

    use module_Green_axi_MB


!-
    implicit none

    integer npt_source, npt_point, OPT_PARALLEL, n_thread

    real(kind=8) r_source(npt_source), z_source(npt_source), I_source(npt_source)
    real(kind=8) r_point(npt_point), z_point(npt_point), vec_flux(npt_source)

!-
    if(OPT_PARALLEL == 0) then
        call sub_Green_flux_serial(npt_source,r_source,z_source,I_source,npt_point,r_point,z_point,vec_flux)

    else if(OPT_PARALLEL == 1) then
        call sub_Green_flux_parallel(npt_source,r_source,z_source,I_source,npt_point,r_point,z_point,n_thread,vec_flux)

    endif


    end subroutine fun_Green_filament_flux_SP_4mex
!------------------------------------------------------------------------------------------------------------------











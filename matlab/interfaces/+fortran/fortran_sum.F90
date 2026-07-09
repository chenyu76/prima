#include "fintrf.h"
!
!     MEX wrapper for Fortran-native SUM.  Branches on ndims (1..7)
!     because Fortran SUM requires compile-time rank: SUM(rank-N, dim)
!     returns rank-(N-1), and gfortran rejects reshape with a shape
!     array whose size is allocatable.  We use array constructors like
!     [dims(1),dims(2)] whose element count is known at compile time.
!     TRANSFER recovers mwSize* dims from real*8 bit-patterns.
!     cnt (mwSize) avoids integer*4/mwSize type mismatch on 64-bit.
!
       subroutine mexFunction(nlhs, plhs, nrhs, prhs)
       implicit none
       mwPointer plhs(*), prhs(*)
       integer nlhs, nrhs
       mwPointer mxGetPr, mxCreateDoubleMatrix
       mwPointer mxGetM, mxGetN
       mwSize m, n, total, cnt
       integer*4 :: dim, ndims, d
       real*8, allocatable :: x(:), outFlat(:)
       real*8 :: dim_arr(1), tmp(1)

       mwPointer, external :: mxGetDimensions
       integer*4, external :: mxGetNumberOfDimensions
       mwPointer pDims
       real*8, allocatable :: dimsReal(:)
       integer*8 :: i8tmp
       integer*4, allocatable :: dims(:)
       mwSize :: nBefore, nAfter, outTotal
       mwPointer px, py

       real*8, allocatable :: x2(:,:), x3(:,:,:)
       real*8, allocatable :: x4(:,:,:,:), x5(:,:,:,:,:)
       real*8, allocatable :: x6(:,:,:,:,:,:)
       real*8, allocatable :: x7(:,:,:,:,:,:,:)
       if (nrhs .lt. 1 .or. nrhs .gt. 2) then
           call mexErrMsgIdAndTxt('fortran:sum:nrhs','Bad nrhs.')
       end if

       m = mxGetM(prhs(1))
       n = mxGetN(prhs(1))
       total = m * n

       if (nrhs .ge. 2) then
           px = mxGetPr(prhs(2))
           cnt = 1
           call mxCopyPtrToReal8(px, dim_arr, cnt)
           dim = int(dim_arr(1), 4)
       else
           dim = 0
       end if

       allocate(x(total))
       px = mxGetPr(prhs(1))
       call mxCopyPtrToReal8(px, x, total)

       if (dim .le. 0) then
           plhs(1) = mxCreateDoubleMatrix(1, 1, 0)
           py = mxGetPr(plhs(1))
           tmp(1) = sum(x)
           cnt = 1
           call mxCopyReal8ToPtr(tmp, py, cnt)
           deallocate(x)
           return
       end if

       ndims = mxGetNumberOfDimensions(prhs(1))
       pDims = mxGetDimensions(prhs(1))
       allocate(dimsReal(ndims))
       cnt = ndims
       call mxCopyPtrToReal8(pDims, dimsReal, cnt)
       allocate(dims(ndims))
       i8tmp = 0
       do d = 1, ndims
           dims(d) = int(transfer(dimsReal(d), i8tmp), 4)
       end do
       deallocate(dimsReal)

       if (dim .gt. ndims) then
           plhs(1) = mxCreateDoubleMatrix(m, n, 0)
           py = mxGetPr(plhs(1))
           call mxCopyReal8ToPtr(x, py, total)
           deallocate(x, dims)
           return
       end if

       nBefore = 1
       do d = 1, dim - 1
           nBefore = nBefore * dims(d)
       end do
       nAfter = 1
       do d = dim + 1, ndims
           nAfter = nAfter * dims(d)
       end do
       outTotal = nBefore * nAfter

       plhs(1) = mxCreateDoubleMatrix(nBefore, nAfter, 0)
       allocate(outFlat(outTotal))

       if (ndims .eq. 1) then
           outFlat(1) = sum(reshape(x,[dims(1)]), 1)
       else if (ndims .eq. 2) then
           allocate(x2(dims(1),dims(2)))
           x2 = reshape(x,[dims(1),dims(2)])
           outFlat = reshape(sum(x2,dim),[outTotal])
           deallocate(x2)
       else if (ndims .eq. 3) then
           allocate(x3(dims(1),dims(2),dims(3)))
           x3 = reshape(x,[dims(1),dims(2),dims(3)])
           outFlat = reshape(sum(x3,dim),[outTotal])
           deallocate(x3)
       else if (ndims .eq. 4) then
           allocate(x4(dims(1),dims(2),dims(3),dims(4)))
           x4 = reshape(x,[dims(1),dims(2),dims(3),dims(4)])
           outFlat = reshape(sum(x4,dim),[outTotal])
           deallocate(x4)
       else if (ndims .eq. 5) then
           allocate(x5(dims(1),dims(2),dims(3),dims(4),dims(5)))
           x5 = reshape(x,[dims(1),dims(2),dims(3),dims(4),dims(5)])
           outFlat = reshape(sum(x5,dim),[outTotal])
           deallocate(x5)
       else if (ndims .eq. 6) then
           allocate(x6(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)))
           x6 = reshape(x,[dims(1),dims(2),dims(3),dims(4),dims(5),dims(6)])
           outFlat = reshape(sum(x6,dim),[outTotal])
           deallocate(x6)
       else if (ndims .eq. 7) then
           allocate(x7(dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7)))
           x7 = reshape(x,[dims(1),dims(2),dims(3),dims(4),dims(5),dims(6),dims(7)])
           outFlat = reshape(sum(x7,dim),[outTotal])
           deallocate(x7)
       else
           call mexErrMsgIdAndTxt('fortran:sum:ndims','>7 dims.')
       end if

       py = mxGetPr(plhs(1))
       call mxCopyReal8ToPtr(outFlat, py, outTotal)
       deallocate(x, outFlat, dims)
       end subroutine mexFunction

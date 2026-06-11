#include "fintrf.h"

      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
      implicit none
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs
      mwPointer mxGetPr, mxCreateDoubleMatrix
      mwPointer mxGetM, mxGetN
      mwSize m, n, total
      real*8, allocatable :: x(:), y(:)
      mwPointer px, py

      if (nrhs .ne. 1) then
          call mexErrMsgIdAndTxt('fortran:sin:nrhs', 'One input required.')
      end if

      m = mxGetM(prhs(1))
      n = mxGetN(prhs(1))
      total = m * n

      plhs(1) = mxCreateDoubleMatrix(m, n, 0)

      allocate(x(total), y(total))

      px = mxGetPr(prhs(1))
      call mxCopyPtrToReal8(px, x, total)

      y = sin(x)

      py = mxGetPr(plhs(1))
      call mxCopyReal8ToPtr(y, py, total)

      deallocate(x, y)
      end subroutine mexFunction

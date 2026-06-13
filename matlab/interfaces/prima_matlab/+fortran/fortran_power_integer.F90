#include "fintrf.h"

      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
      implicit none
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs
      mwPointer mxGetPr, mxCreateDoubleMatrix
      mwPointer mxGetM, mxGetN
      mwSize m, n, total
      real*8, allocatable :: x(:), y(:)
      real*8 :: exp_arr(1)
      integer*4 :: exp_int
      mwPointer px, py, pexp

      if (nrhs .ne. 2) then
          call mexErrMsgIdAndTxt('fortran:power_integer:nrhs', 'Two inputs required.')
      end if

      m = mxGetM(prhs(1))
      n = mxGetN(prhs(1))
      total = m * n

      plhs(1) = mxCreateDoubleMatrix(m, n, 0)

      pexp = mxGetPr(prhs(2))
      call mxCopyPtrToReal8(pexp, exp_arr, 1)
      exp_int = int(exp_arr(1), 4)

      allocate(x(total), y(total))

      px = mxGetPr(prhs(1))
      call mxCopyPtrToReal8(px, x, total)

      y = x**exp_int

      py = mxGetPr(plhs(1))
      call mxCopyReal8ToPtr(y, py, total)

      deallocate(x, y)
      end subroutine mexFunction

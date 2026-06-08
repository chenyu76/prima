classdef linalg_mod
    %--------------------------------------------------------------------------------------------------%
    % This module provides some basic linear algebra procedures.
    %
    % The procedures are NOT intended to be optimized but to be sufficient for my projects. The projects
    % are mainly the development and maintenance of derivative-free optimization software, where the
    % major expense comes from the function evaluations, NOT the numerical linear algebraic computations,
    % and the sizes of matrices/vectors involved are relatively SMALL, the order being at most 10^3.
    %
    % If your needs are of a different nature, you may still use these procedures to prototype your
    % ideas, but keep in mind that the implementations here are mostly STRAIGHTFORWARD and NAIVE, and
    % some algorithms selected here may be suboptimal for your problems.
    %
    % If it is needed to enhance the performance of these procedures, one can optimize their
    % implementations according to the resources (hardware, e.g., C/GPU, cache, and libraries, e.g.,
    % BLAS, LAPACK) available and the sizes of the matrices/vectors concerned.
    %
    % In case you need similar procedures in MATLAB/Python/Julia/R, note the following.
    % 1. Most of the procedures here are intrinsic to the languages or available in standard libraries.
    % If available, they should NOT be implemented from scratch like we do here.
    % 2. For the procedures that are not available, it may be better to code them inline instead of as
    % external functions, because the code is usually short using matrix/vector operations, and because
    % the overhead of function calling can be high in these languages.
    % 3. In Fortran, we implement the procedures as subroutines/functions here for several reasons.
    % 3.1.) Most of the procedures are not intrinsically available in Fortran.
    % 3.2.) When using these procedures for the modernization of Powell's derivative-free software, we
    % want to start with an implementation that is verifiably faithful to Powell's original code. To
    % achieve such faithfulness, it is not always possible to use the intrinsic matrix/vector procedures
    % in Fortran, the most noticeable examples being DOT_PRODUCT (v.s. INPROD) and MATMUL (v.s. MATPROD).
    % Powell implemented all matrix/vector operations by loops, which may not be the case for intrinsic
    % procedures such as MATMUL and DOT_PRODUCT. Different implementations lead to slightly different
    % results due to rounding, and hence the verification of faithfulness will fail.
    % 3.3.) As of 20220507, with some compilers, the performance of Fortran's intrinsic matrix/vector
    % procedures may not be as good as naive loops, let alone highly optimized libraries like BLAS.
    % Concentrating all the linear algebra procedures at one place as we do here, it will be relatively
    % easy to optimize them when necessary.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net).
    %
    % Started: July 2020
    %
    % Last Modified: Fri 13 Feb 2026 05:11:41 PM CET
    %--------------------------------------------------------------------------------------------------%

    methods
        function obj = linalg_mod()
            % Mathematically, INPROD = DOT_PRODUCT, MATPROD = MATMUL

        end
        function varargout = matprod(obj, varargin)
            if numel(varargin) == 2 && isvector(varargin{1}) && (~isvector(varargin{2}) && ~isscalar(varargin{2}))
                [varargout{1:nargout}] = obj.matprod12(varargin{:});
            elseif numel(varargin) == 2 && (~isvector(varargin{1}) && ~isscalar(varargin{1})) && isvector(varargin{2})
                [varargout{1:nargout}] = obj.matprod21(varargin{:});
            else
                [varargout{1:nargout}] = obj.matprod22(varargin{:});
            end
        end
        function varargout = r1update(obj, varargin)
            if numel(varargin) == 3
                [varargout{1:nargout}] = obj.r1_sym(varargin{:});
            else
                [varargout{1:nargout}] = obj.r1(varargin{:});
            end
        end
        function varargout = r2update(obj, varargin)
            if numel(varargin) == 4
                [varargout{1:nargout}] = obj.r2_sym(varargin{:});
            else
                [varargout{1:nargout}] = obj.r2(varargin{:});
            end
        end
        function varargout = eye(obj, varargin)
            if numel(varargin) == 1
                [varargout{1:nargout}] = obj.eye1(varargin{:});
            else
                [varargout{1:nargout}] = obj.eye2(varargin{:});
            end
        end
        function varargout = project(obj, varargin)
            if numel(varargin) == 2 && isvector(varargin{2})
                [varargout{1:nargout}] = obj.project1(varargin{:});
            else
                [varargout{1:nargout}] = obj.project2(varargin{:});
            end
        end
        function varargout = lsqr(obj, varargin)
            if numel(varargin) == 3 && isvector(varargin{1}) && (~isvector(varargin{2}) && ~isscalar(varargin{2}))
                [varargout{1:nargout}] = obj.lsqr_Rfull(varargin{:});
            else
                [varargout{1:nargout}] = obj.lsqr_Rdiag(varargin{:});
            end
        end
        function varargout = isminor(obj, varargin)
            if numel(varargin) == 2 && isscalar(varargin{1}) && isscalar(varargin{2})
                [varargout{1:nargout}] = obj.isminor0(varargin{:});
            else
                [varargout{1:nargout}] = obj.isminor1(varargin{:});
            end
        end
        function varargout = sort(obj, varargin)
            if numel(varargin) >= 1 && numel(varargin) <= 2 && isvector(varargin{1}) && (numel(varargin) < 2 || (ischar(varargin{2})))
                [varargout{1:nargout}] = obj.sort_i1(varargin{:});
            else
                [varargout{1:nargout}] = obj.sort_i2(varargin{:});
            end
        end
        function varargout = minimum(obj, varargin)
            if numel(varargin) == 1 && isvector(varargin{1})
                [varargout{1:nargout}] = obj.minimum1(varargin{:});
            else
                [varargout{1:nargout}] = obj.minimum2(varargin{:});
            end
        end
        function varargout = maximum(obj, varargin)
            if numel(varargin) == 1 && isvector(varargin{1})
                [varargout{1:nargout}] = obj.maximum1(varargin{:});
            else
                [varargout{1:nargout}] = obj.maximum2(varargin{:});
            end
        end
        function varargout = norm(obj, varargin)
            if numel(varargin) >= 1 && numel(varargin) <= 2 && isvector(varargin{1}) && (numel(varargin) < 2 || (isfloat(varargin{2})))
                [varargout{1:nargout}] = obj.p_norm(varargin{:});
            elseif numel(varargin) == 2 && isvector(varargin{1}) && ischar(varargin{2})
                [varargout{1:nargout}] = obj.named_norm_vec(varargin{:});
            else
                [varargout{1:nargout}] = obj.named_norm_mat(varargin{:});
            end
        end
        function varargout = linspace(obj, varargin)
            if numel(varargin) == 3 && isfloat(varargin{1}) && isfloat(varargin{2})
                [varargout{1:nargout}] = obj.linspace_r(varargin{:});
            else
                [varargout{1:nargout}] = obj.linspace_i(varargin{:});
            end
        end
        function varargout = hessenberg(obj, varargin)
            if numel(varargin) == 3 && isvector(varargin{2}) && isvector(varargin{3})
                [varargout{1:nargout}] = obj.hessenberg_hhd_trid(varargin{:});
            else
                [varargout{1:nargout}] = obj.hessenberg_full(varargin{:});
            end
        end
        function varargout = eigmin(obj, varargin)
            if numel(varargin) >= 2 && numel(varargin) <= 3
                [varargout{1:nargout}] = obj.eigmin_sym_trid(varargin{:});
            end
        end
        function varargout = int(obj, varargin)
            if numel(varargin) == 1
                [varargout{1:nargout}] = obj.logical_to_int(varargin{:});
            end
        end
        function A = r1_sym(obj, A, alpha, x)
            %--------------------------------------------------------------------------------------------------%
            % R1_SYM sets
            % A = A + ALPHA*( X*X^T ),
            % where A is an NxN matrix, ALPHA is a scalar, and X is an N-dimensional vector.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();

            % Inputs


            % In-outputs
            % A(SIZE(X), SIZE(X))
            % Local variables
            srname = "R1_SYM";
            n = NaN; j = NaN;

            % Sizes
            n = fix(numel(x));

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(size(A, 1) == n && size(A, 2) == n, "SIZE(A) == [N, N]", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            % Only update the LOWER TRIANGULAR part of A.
            for j = 1:n
                A(j:n, j) = A(j:n, j) + alpha * x(j:n) * x(j);
            end
            A = obj.symmetrize(A); % Copy A(LOWER_TRI) to A(UPPER_TRI).

            % For some reason, A + alpha*outprod(x,x), A + (outprod(alpha*x, x) + outprod(x, alpha*x))/2,
            % A + symmetrize(x, alpha*x), or A + sign(alpha) * outprod(sqrt(|alpha|) * x, sqrt(|alpha|) * x)
            % does not work as well as the above lines in NEWUOA, where SYMMETRIZE should copy A(LOWER_TRI)
            % to A(UPPER_TRI) rather than set A = (A'+A)/2. When X is rather small or large, calculating
            % OUTPROD(X, X) can be a bad idea, even though it guarantees symmetry in finite-precision arithmetic.

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(size(A, 1) == n && obj.issymmetric(A), "A is N-by-N and symmetric", srname);
            end
        end
        function A = r1(obj, A, alpha, x, y)
            %--------------------------------------------------------------------------------------------------%
            % R1 sets
            % A = A + ALPHA*( X*Y^T ),
            % where A is an MxN matrix, ALPHA is a real scalar, X is an M-dimensional vector, and Y is an
            % N-dimensional vector.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();

            % Inputs



            % In-outputs
            % A(SIZE(X), SIZE(Y))
            % Local variables
            srname = "R1";

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(size(A, 1) == numel(x) && size(A, 2) == numel(y), "SIZE(A) == [SIZE(X), SIZE(Y)]", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            % N.B.: The use of OUTPROD is expensive memory-wise, but it is not our concern in this implementation.
            A(:, :) = A + obj.outprod(alpha * x, y);
            %A = A + alpha * outprod(x, y)

            %====================%
            %  Calculation ends  %
            %====================%
        end
        function A = r2_sym(obj, A, alpha, x, y)
            %--------------------------------------------------------------------------------------------------%
            % R2_SYM sets
            % A = A + ALPHA*( X*Y^T + Y*X^T ),
            % where A is an NxN matrix, X and Y are N-dimensional vectors, and alpha is a scalar.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();

            % Inputs



            % In-outputs
            % A(SIZE(X), SIZE(X))
            % Local variables
            srname = "R2_SYM";
            n = NaN; j = NaN;

            % Sizes
            n = fix(numel(x));

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(y) == n, "SIZE(Y) == N", srname);
                debug_obj.assert(size(A, 1) == n && size(A, 2) == n, "SIZE(A) == [N, N]", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            for j = 1:n
                A(j:n, j) = A(j:n, j) + alpha * x(j:n) * y(j) + alpha * y(j:n) * x(j);
            end
            A = obj.symmetrize(A); % Copy A(LOWER_TRI) to A(UPPER_TRI).

            % For some reason, A = A + ALPHA * (OUTPROD(X, Y) + OUTPROD(Y, X)) does not work as well as the
            % above lines for NEWUOA, where SYMMETRIZE should copy A(LOWER_TRI) to A(UPPER_TRI), although
            % ALPHA*( X*Y^T + Y*X^T) is guaranteed symmetric even in floating-point arithmetic.

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(obj.issymmetric(A), "A is symmetric", srname);
            end
        end
        function A = r2(obj, A, alpha, x, y, beta, u, v)
            %--------------------------------------------------------------------------------------------------%
            % R2 sets
            % A = A + ( ALPHA*( X*Y^T ) + BETA*( U*V^T ) ),
            % where A is an MxN matrix, ALPHA and BETA are real scalars, X and U are M-dimensional vectors,
            % Y and V are N-dimensional vectors.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();

            % Inputs




            % U(SIZE(X))
            % V(SIZE(Y))
            % In-outputs
            % A(SIZE(X), SIZE(Y))
            % Local variables
            srname = "R2";

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(u) == numel(x), "SIZE(U) == SIZE(X)", srname);
                debug_obj.assert(numel(v) == numel(y), "SIZE(V) == SIZE(Y)", srname);
                debug_obj.assert(size(A, 1) == numel(x) && size(A, 2) == numel(y), "SIZE(A) == [SIZE(X), SIZE(Y)]", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            % N.B.: The use of OUTPROD is expensive memory-wise, but it is not our concern in this implementation.
            A(:, :) = A + obj.outprod(alpha * x, y) + obj.outprod(beta * u, v);
            %A = A + (alpha * outprod(x, y) + beta * outprod(u, v))

            %====================%
            %  Calculation ends  %
            %====================%
        end
        function z = matprod12(obj, x, y)
            %--------------------------------------------------------------------------------------------------%
            % This procedure calculates the matrix product of X and Y, where X is an M-dimensional vector
            % considered as a row, and Y is an M-by-N matrix.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();

            % Inputs


            % Outputs
            z = NaN(size(y, 2), 1);
            % Local variables
            srname = "MATPROD12";
            j = NaN;

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(x) == size(y, 1), "SIZE(X) == SIZE(Y, 1)", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            for j = 1:size(y, 2)
                % When interfaced with MATLAB, the following seems more efficient than a loop, which is strange
                % since inprod itself is implemented by a loop. This may depend on the machine (e.g., cache
                % size), compiler, compiling options, and MATLAB version.
                z(j) = obj.inprod(x, y(:, j));
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(z) == size(y, 2), "SIZE(Z) == SIZE(Y, 2)", srname);
            end
        end
        function z = matprod21(~, x, y)
            %--------------------------------------------------------------------------------------------------%
            % This procedure calculates the matrix product of X and Y, where X is an M-by-N matrix, and Y is an
            % M-dimensional vector considered as a column.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();

            % Inputs


            % Outputs
            z = NaN(size(x, 1), 1);
            % Local variables
            srname = "MATPROD21";
            j = NaN;

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(size(x, 2) == numel(y), "SIZE(X, 2) == SIZE(Y)", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            z = repmat(consts_obj.ZERO, size(z));
            for j = 1:size(x, 2)
                z(:) = z + x(:, j) * y(j);
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(z) == size(x, 1), "SIZE(Z) == SIZE(X, 1)", srname);
            end
        end
        function z = matprod22(~, x, y)
            %--------------------------------------------------------------------------------------------------%
            % This procedure calculates the matrix product of X and Y, where X is an M-by-P matrix, and Y is a
            % P-by-N matrix.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();

            % Inputs


            % Outputs
            z = NaN(size(x, 1), size(y, 2));
            % Local variables
            srname = "MATPROD22";
            i = NaN; j = NaN;

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(size(x, 2) == size(y, 1), "SIZE(X, 2) == SIZE(Y, 1)", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            z = repmat(consts_obj.ZERO, size(z));
            for j = 1:size(y, 2)
                for i = 1:size(x, 2)
                    z(:, j) = z(:, j) + x(:, i) * y(i, j);
                end
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(size(z, 1) == size(x, 1) && size(z, 2) == size(y, 2), "[SIZE(Z) == SIZE(X, 1), SIZE(Y, 2)]", srname);
            end
        end
        function z = inprod(~, x, y)
            %--------------------------------------------------------------------------------------------------%
            % INPROD calculates the inner product of X and Y, i.e., Z = X^T*Y, regarding X and Y as columns.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();

            % Inputs


            % Outputs
            z = NaN;
            % Local variables
            srname = "INPROD";
            i = NaN;

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(x) == numel(y), "SIZE(X) == SIZE(Y)", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            z = consts_obj.ZERO;
            for i = 1:fix(numel(x))
                z = z + x(i) * y(i);
            end

            %====================%
            %  Calculation ends  %
            %====================%
        end
        function z = outprod(~, x, y)
            %--------------------------------------------------------------------------------------------------%
            % OUTPROD calculates the outer product of X and Y, i.e., Z = X*Y^T, regarding X and Y as columns.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();

            % Inputs


            % Outputs
            z = NaN(numel(x), numel(y));
            % Local variables
            srname = "OUTPROD";
            i = NaN;

            %====================%
            % Calculation starts %
            %====================%

            for i = 1:fix(numel(y))
                z(:, i) = x * y(i);
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(size(z, 1) == numel(x) && size(z, 2) == numel(y), "SIZE(Z) == [SIZE(X), SIZE(Y)]", srname);
            end
        end
        function x = eye1(~, n)
            %--------------------------------------------------------------------------------------------------%
            % EYE1 is the univariate case of EYE, a function similar to the MATLAB function with the same name.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();

            % Inputs

            % Outputs
            x = NaN(max(n, 0));
            % Local variables
            srname = "EYE1";
            i = NaN;

            %====================%
            % Calculation starts %
            %====================%

            if size(x, 1) * size(x, 2) > 0
                x = repmat(consts_obj.ZERO, size(x));
                for i = 1:fix(min(size(x, 1), size(x, 2)))
                    x(i, i) = consts_obj.ONE;
                end
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(size(x, 1) == max(n, 0) && size(x, 2) == max(n, 0), "SIZE(X) == [N, N]", srname);
            end
        end
        function x = eye2(~, m, n)
            %--------------------------------------------------------------------------------------------------%
            % EYE2 is the bivariate case of EYE, a function similar to the MATLAB function with the same name.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();

            % Inputs


            % Outputs
            x = NaN(max(m, 0), max(n, 0));
            % Local variables
            srname = "EYE2";
            i = NaN;

            %====================%
            % Calculation starts %
            %====================%

            if size(x, 1) * size(x, 2) > 0
                x = repmat(consts_obj.ZERO, size(x));
                for i = 1:fix(min(size(x, 1), size(x, 2)))
                    x(i, i) = consts_obj.ONE;
                end
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(size(x, 1) == max(m, 0) && size(x, 2) == max(n, 0), "SIZE(X) == [M, N]", srname);
            end
        end
        function x = solve(obj, A, b)
            %--------------------------------------------------------------------------------------------------%
            % This function solves the linear system A*X = B. We assume that A is a square matrix that is small
            % and invertible, and B is a vector of length SIZE(A, 1). The implementation is NAIVE.
            % TODO: Better to implement it into several subfunctions: triu, tril, and general square.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();

            % Inputs


            % Outputs
            x = NaN(size(A, 2), 1);
            % Local variables
            srname = "SOLVE";
            P = NaN(size(A, 1), 1);
            i = NaN;
            n = NaN;
            Q = NaN(size(A, 1));
            R = NaN(size(A, 1), size(A, 2));
            tol = NaN;

            % Sizes
            n = size(A, 1);

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(size(A, 1) == size(A, 2), "A is square", srname);
                debug_obj.assert(numel(b) == size(A, 1), "SIZE(B) == SIZE(A, 1)", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            if n <= 0
                % Of course, N < 0 should never happen.
                return
            end

            % Zaikun 20220527: With the following code, Huawei Bisheng flang 2.1.0, Arm Fortran Compiler 23.1,
            % and AOCC 5.1 flang, which raise a false positive error about out-bound subscripts when invoked
            % with the -Mbounds flag. See https://github.com/flang-compiler/flang/issues/1238
            if obj.istril(A)
                for i = 1:n
                    x(i) = (b(i) - obj.inprod(A(i, 1:i - 1), x(1:i - 1))) / A(i, i); % INPROD = 0 if I == 1.
                end
            elseif obj.istriu(A)
                % This case is invoked in LINCOA.
                for i = n:-1:1
                    x(i) = (b(i) - obj.inprod(A(i, i + 1:n), x(i + 1:n))) / A(i, i); % INPROD = 0 if I == N.
                end
            else
                % This is NOT a good algorithm for linear systems, but since the QR subroutine is available ...
                [Q, R, P] = obj.qr(A);
                x(:) = obj.matprod12(b, Q);
                for i = n:-1:1
                    x(i) = (x(i) - obj.inprod(R(i, i + 1:n), x(i + 1:n))) / R(i, i); % INPROD = 0 if I == N.
                end
                x(P) = x; % Handle the permutation.
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(x) == size(A, 2), "SIZE(X) == SIZE(A, 2)", srname);
                if infnan_obj.is_finite(sum(abs(A), 'all') + sum(abs(b), 'all'))
                    tol = max(consts_obj.TEN ^ max(-8, -consts_obj.MAXPOW10), min(0.1, consts_obj.TEN ^ min(8, consts_obj.MAXPOW10) * consts_obj.EPS * double(n + 1)));
                    debug_obj.assert(obj.p_norm(obj.matprod21(A, x) - b) <= tol * max([consts_obj.ONE, obj.p_norm(b), obj.p_norm(x)], [], 'all'), "A*X == B", srname);
                end
            end
        end
        function B = inv(obj, A)
            %--------------------------------------------------------------------------------------------------%
            % This function calculates the inverse of a matrix A, which is ASSUMED TO BE SMALL AND INVERTIBLE.
            % The function is implemented NAIVELY. It is NOT coded for general purposes but only for the usage
            % in this project. Indeed, only the lower triangular case is used.
            % TODO: extend this function to calculate the pseudo inverse of any matrix of full rank. Better to
            % implement it into several subfunctions: triu with M >= N, tril with M <= N; general with M >= N,
            % general with M <= N, etc.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();

            % Inputs

            % Outputs
            B = NaN(size(A, 1));
            % Local variables
            srname = "INV";
            P = NaN(size(A, 1), 1);
            InvP = NaN(size(A, 1), 1);
            i = NaN;
            n = NaN;
            Q = NaN(size(A, 1));
            R = NaN(size(A, 1));
            tol = NaN;

            % Sizes
            n = size(A, 1);

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(size(A, 1) == size(A, 2), "A is square", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            if n <= 0
                % Of course, N < 0 should never happen.
                return
            end

            if obj.istril(A)
                % This case is invoked in COBYLA.
                R(:, :) = A'; % Take transpose to work on columns.
                B = repmat(consts_obj.ZERO, size(B));
                for i = 1:n
                    B(i, i) = consts_obj.ONE / R(i, i);
                    B(1:i - 1, i) = -obj.matprod21(B(1:i - 1, 1:i - 1), R(1:i - 1, i) ./ R(i, i));
                end
                B(:, :) = B';
            elseif obj.istriu(A)
                B = repmat(consts_obj.ZERO, size(B));
                for i = 1:n
                    B(i, i) = consts_obj.ONE / A(i, i);
                    B(1:i - 1, i) = -obj.matprod21(B(1:i - 1, 1:i - 1), A(1:i - 1, i) ./ A(i, i));
                end
            else
                % This is NOT the best algorithm for the inverse, but since the QR subroutine is available ...
                [Q, R, P] = obj.qr(A);
                R(:, :) = R'; % Take transpose to work on columns.
                B = repmat(consts_obj.ZERO, size(B));
                for i = n:-1:1
                    B(:, i) = (Q(:, i) - obj.matprod21(B(:, i + 1:n), R(i + 1:n, i))) ./ R(i, i);
                end
                InvP(P) = obj.linspace_i(1, n, n); % The inverse permutation
                B(:, :) = B(:, InvP)';
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(size(B, 1) == n && size(B, 2) == n, "SIZE(B) == [N, N]", srname);
                debug_obj.assert(obj.istril(B) || ~obj.istril(A), "If A is lower triangular, then so is B", srname);
                debug_obj.assert(obj.istriu(B) || ~obj.istriu(A), "If A is upper triangular, then so is B", srname);
                tol = max(consts_obj.TEN ^ max(-8, -consts_obj.MAXPOW10), min(0.1, consts_obj.TEN ^ min(10, consts_obj.MAXPOW10) * consts_obj.EPS * double(n + 1)));
                debug_obj.assert(obj.isinv(A, B, 'tol', tol), "B = A^{-1}", srname);
            end
        end
        function is_inv = isinv(obj, A, B, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This procedure tests whether A = B^{-1} up to the tolerance TOL.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();

            % Inputs



            % Outputs
            is_inv = false;
            % Local variables
            srname = "ISINV";
            tol_loc = NaN;
            n = NaN;

            % Sizes
            n = size(A, 1);

            % Preconditions
            ipObj = inputParser();
            addParameter(ipObj, 'tol', NaN);
            parse(ipObj, varargin{:});
            tol = ipObj.Results.tol;
            if consts_obj.DEBUGGING
                debug_obj.assert(size(A, 1) == size(A, 2), "A is square", srname);
                debug_obj.assert(size(B, 1) == size(B, 2), "B is square", srname);
                debug_obj.assert(size(A, 1) == size(B, 1), "SIZE(A) == SIZE(B)", srname);
                if ~ismember('tol', ipObj.UsingDefaults)
                    debug_obj.assert(tol >= 0, "TOL >= 0", srname);
                end
            end

            %====================%
            % Calculation starts %
            %====================%

            if ismember('tol', ipObj.UsingDefaults)
                tol_loc = min(1.0e-3, 100.0 * consts_obj.EPS * double(max(size(A, 1), size(A, 2))));
            else
                tol_loc = tol;
            end
            tol_loc = max([tol_loc, tol_loc * max(abs(A), [], 'all'), tol_loc * max(abs(B), [], 'all')], [], 'all');
            is_inv = all(abs(obj.matprod22(A, B) - obj.eye1(n)) <= tol_loc, 'all') || all(abs(obj.matprod22(B, A) - obj.eye1(n)) <= tol_loc, 'all');

            %====================%
            %  Calculation ends  %
            %====================%
        end
        function [Q, R, P] = qr(obj, A, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine calculates the QR factorization of A, possibly with column pivoting, so that
            % A = Q*R (if no pivoting) or A(:, P) = Q*R (if pivoting), where the columns of Q are orthonormal,
            % and R is upper triangular.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();

            % Inputs

            % Outputs



            % Local variables
            srname = "QR";
            pivot = false;
            i = NaN;
            j = NaN;
            k = NaN;
            m = NaN;
            n = NaN;
            G = NaN(2);
            Q_loc = NaN(size(A, 1));
            T = NaN(size(A, 2), size(A, 1));
            tol = NaN;

            ipObj = inputParser();
            addParameter(ipObj, 'Q', NaN);
            addParameter(ipObj, 'R', NaN);
            addParameter(ipObj, 'P', NaN);
            parse(ipObj, varargin{:});
            Q = ipObj.Results.Q;
            R = ipObj.Results.R;
            P = ipObj.Results.P;
            if ~(nargout >= 1 || nargout >= 2 || nargout >= 2)
                return
            end

            % Sizes
            m = size(A, 1);
            n = size(A, 2);

            % Preconditions
            if consts_obj.DEBUGGING
                if nargout >= 1
                    debug_obj.assert(size(Q, 1) == m && (size(Q, 2) == m || size(Q, 2) == min(m, n)), "SIZE(Q) == [M, N] .or. SIZE(Q) == [M, MIN(M, N)]", srname);
                end
                if nargout >= 2
                    debug_obj.assert((size(R, 1) == m || size(R, 1) == min(m, n)) && size(R, 2) == n, "SIZE(R) == [M, N] .or. SIZE(R) == [MIN(M, N), N]", srname);
                end
                if nargout >= 1 && nargout >= 2
                    debug_obj.assert(size(Q, 2) == size(R, 1), "SIZE(Q, 2) == SIZE(R, 1)", srname);
                end
                if nargout >= 3
                    debug_obj.assert(numel(P) == n, "SIZE(P) == N", srname);
                end
            end

            %====================%
            % Calculation starts %
            %====================%

            pivot = (nargout >= 3);
            Q_loc(:, :) = obj.eye1(m);
            T(:, :) = A'; % T is the transpose of R. We consider T in order to work on columns.
            if pivot
                P(:) = obj.linspace_i(1, n, n);
            end

            for j = 1:n
                if pivot
                    k = fix(fortran.maxloc(sum(T(j:n, j:m) .^ 2, 2), 'dim', 1));
                    if k > 1 && k <= n - j + 1
                        k = k + j - 1;
                        P([j, k]) = P([k, j]);
                        T([j, k], :) = T([k, j], :);
                    end
                end
                for i = m:-1:j + 1
                    G(:, :) = obj.planerot(T(j, [j, i]))';
                    T(j, [j, i]) = [obj.hypotenuse(T(j, j), T(j, i)), consts_obj.ZERO]; %T(j, [j, i]) = [sqrt(T(j, j)**2 + T(j, i)**2), ZERO]
                    T(j + 1:n, [j, i]) = obj.matprod22(T(j + 1:n, [j, i]), G);
                    Q_loc(:, [j, i]) = obj.matprod22(Q_loc(:, [j, i]), G);
                end
            end

            if nargout >= 1
                Q(:, :) = Q_loc(:, 1:size(Q, 2));
            end
            if nargout >= 2
                R(:, :) = T(:, 1:size(R, 1))';
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                tol = max(consts_obj.TEN ^ max(-10, -consts_obj.MAXPOW10), min(0.1, consts_obj.TEN ^ min(4, consts_obj.MAXPOW10) * consts_obj.EPS * double(max(m, n) + 1)));
                debug_obj.assert(obj.isorth(Q_loc, 'tol', tol), "The columns of Q are orthonormal", srname);
                debug_obj.assert(obj.istril(T, 'tol', tol), "R is upper triangular", srname);
                if pivot
                    debug_obj.assert(all(abs(obj.matprod22(Q_loc, T') - A(:, P)) <= max(tol, tol * max(abs(A), [], 'all')), 'all'), "A(:, P) == Q*R", srname);
                    for j = 1:min(m, n) - 1
                        % The following test cannot be passed on ill-conditioned problems.
                        %call assert(abs(T(j, j)) + max(tol, tol * abs(T(j, j))) >= &
                        % & abs(T(j + 1, j + 1)), '|R(J, J)| >= |R(J + 1, J + 1)|', srname)
                        debug_obj.assert(all(T(j, j) ^ 2 + max(tol, tol * T(j, j) ^ 2) >= sum(T(j + 1:n, j:min(m, n)) .^ 2, 2), 'all'), "R(J, J)^2 >= SUM(R(J : MIN(M, N), J + 1 : N).^2", srname);
                    end
                else
                    debug_obj.assert(all(abs(obj.matprod22(Q_loc, T') - A) <= max(tol, tol * max(abs(A), [], 'all')), 'all'), "A == Q*R", srname);
                end
            end
        end
        function x = lsqr_Rdiag(obj, A, b, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This function solves the linear least squares problem min ||A*x - b||_2 by the QR factorization.
            % This function is used in COBYLA, where,
            % 1. Q is supplied externally (called Z);
            % 2. RDIAG (the diagonal of R) is supplied externally (called ZDOTA);
            % 3. A HAS FULL COLUMN RANK;
            % 4. It seems that b (CGRAD and DNEW) is in the column space of A (not sure yet).
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();

            % Inputs
            % A(M, N)
            % B(M)
            % Q(M, :), SIZE(Q, 2) = M or MIN(M, N)
            % Rdiag(MIN(M, N))
            % Outputs
            x = NaN(size(A, 2), 1);
            % Local variables
            srname = "LSQR_RDIAG";
            pivot = false;
            i = NaN;
            j = NaN;
            m = NaN;
            n = NaN;
            P = NaN(size(A, 2), 1);
            rank = NaN;
            Q_loc = NaN(size(A, 1), min(size(A, 1), size(A, 2)));
            Rdiag_loc = NaN(min(size(A, 1), size(A, 2)), 1);
            tol = NaN;
            y = NaN(numel(b), 1);
            yq = NaN;
            yqa = NaN;

            % Sizes
            m = size(A, 1);
            n = size(A, 2);

            % Preconditions
            ipObj = inputParser();
            addParameter(ipObj, 'Q', NaN);
            addParameter(ipObj, 'Rdiag', NaN);
            parse(ipObj, varargin{:});
            Q = ipObj.Results.Q;
            Rdiag = ipObj.Results.Rdiag;
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(b) == m, "SIZE(B) == M", srname);
                if ~ismember('Q', ipObj.UsingDefaults)
                    debug_obj.assert(size(Q, 1) == m && (size(Q, 2) == m || size(Q, 2) == min(m, n)), "SIZE(Q) == [M, N] .or. SIZE(Q) == [M, MIN(M, N)]", srname);
                    tol = max(consts_obj.TEN ^ max(-10, -consts_obj.MAXPOW10), min(0.1, consts_obj.TEN ^ min(6, consts_obj.MAXPOW10) * consts_obj.EPS * double(max(m, n) + 1)));
                    debug_obj.assert(obj.isorth(Q, 'tol', tol), "The columns of Q are orthogonal", srname);
                end
                if ~ismember('Rdiag', ipObj.UsingDefaults)
                    debug_obj.assert(numel(Rdiag) == min(m, n), "SIZE(R) == MIN(M, N)", srname);
                    debug_obj.assert(~ismember('Q', ipObj.UsingDefaults), "Rdiag is present only if Q is present", srname);
                end
            end

            %====================%
            % Calculation starts %
            %====================%

            if n <= 0
                % Of course, N < 0 should never happen.
                return
            end

            if ismember('Q', ipObj.UsingDefaults)
                [Q_loc, ~, P] = obj.qr(A);
                Rdiag_loc(:) = reshape(cell2mat(arrayfun(@(i) obj.inprod(Q_loc(:, i), A(:, P(i))), (1:min(m, n)), "UniformOutput", false)), [], 1);
                %%MATLAB: Rdiag_loc = sum(Q_loc(:, 1:min(m,n)) .* A(:, P(1:min(m,n))), 1); % Row vector
                rank = max([0; reshape(obj.trueloc(abs(Rdiag_loc) > 0), [], 1)], [], 'all');
                pivot = true;
            else
                Q_loc(:, :) = Q(:, 1:size(Q_loc, 2));
                if ismember('Rdiag', ipObj.UsingDefaults)
                    Rdiag_loc(:) = reshape(cell2mat(arrayfun(@(i) obj.inprod(Q_loc(:, i), A(:, i)), (1:min(m, n)), "UniformOutput", false)), [], 1);
                    %%MATLAB: Rdiag_loc = sum(Q_loc(:, 1:min(m,n)) .* A(:, 1:min(m,n)), 1); % Row vector
                else
                    Rdiag_loc(:) = Rdiag;
                end
                rank = min(m, n);
                pivot = false;
            end

            x = repmat(consts_obj.ZERO, size(x));
            y(:) = b; % Local copy of B; B is INTENT(IN) and should not be modified.

            for i = rank:-1:1
                if pivot
                    j = P(i);
                else
                    j = i;
                end
                % The following IF comes from Powell. It forces X(J) = 0 if deviations from this value can be
                % attributed to computer rounding errors. This is a favorable choice in the context of COBYLA.
                yq = obj.inprod(y, Q_loc(:, i));
                yqa = obj.inprod(abs(y), abs(Q_loc(:, i)));
                if obj.isminor0(yq, yqa)
                    x(j) = consts_obj.ZERO;
                else
                    x(j) = yq / Rdiag_loc(i);
                    y(:) = y - x(j) * A(:, j);
                end
            end

            %====================%
            %  Calculation ends  %
            %====================%

            %% Postconditions
            %if (DEBUGGING) then
            % % The following test cannot be passed.
            % %call assert(norm(matprod(b - matprod(A, x), A)) <= max(tol, tol * norm(matprod(b, A))), &
            % % & 'A*X is the projection of B to the column space of A', srname)
            %end if
        end
        function x = lsqr_Rfull(obj, b, Q, R)
            %--------------------------------------------------------------------------------------------------%
            % This function solves the linear least squares problem min ||A*x - b||_2 by the QR factorization.
            % This function is used in LINCOA, where,
            % 1. The economy-size QR factorization is supplied externally (Q is called QFAC and R is called RFAC);
            % 2. R is non-singular.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();

            % Inputs
            % B(M)
            % Q(M, N)
            % R(N, N)
            % Outputs
            x = NaN(size(R, 2), 1);
            % Local variables
            srname = "LSQR_RFULL";
            i = NaN;
            j = NaN;
            m = NaN;
            n = NaN;
            tol = NaN;

            % Sizes
            m = size(Q, 1);
            n = size(R, 2);

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(m >= n && n >= 0, "M >= N >= 0", srname);
                debug_obj.assert(numel(b) == m, "SIZE(B) == M", srname);
                debug_obj.assert(size(Q, 1) == m && size(Q, 2) == n, "SIZE(Q) == [M, N]", srname);
                tol = max(consts_obj.TEN ^ max(-10, -consts_obj.MAXPOW10), min(0.1, consts_obj.TEN ^ min(6, consts_obj.MAXPOW10) * consts_obj.EPS * double(m + 1)));
                debug_obj.assert(obj.isorth(Q, 'tol', tol), "The columns of Q are orthogonal", srname);
                debug_obj.assert(size(R, 1) == n && size(R, 2) == n, "SIZE(R) == [N, N]", srname);
                debug_obj.assert(obj.istriu(R), "R is upper triangular", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            if n <= 0
                % Of course, N < 0 should never happen.
                return
            end

            x(:) = obj.matprod12(b, Q);
            for i = n:-1:1
                for j = i + 1:n
                    x(i) = x(i) - R(i, j) * x(j);
                end
                x(i) = x(i) / R(i, i);
            end
            %--------------------------------------------------------------------------------------------------%
            % The following is equivalent to the above, yet the above version works slightly better in LINCOA.
            % %do i = n, 1_IK, -1_IK
            % %    x(i) = (inprod(Q(:, i), b) - inprod(R(i, i + 1:n), x(i + 1:n))) / R(i, i)
            % %end do
            %--------------------------------------------------------------------------------------------------%

            %====================%
            %  Calculation ends  %
            %====================%
        end
        function D = diag(~, A, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This function takes the K-th diagonal of the matrix A, K = 0 (default) corresponding to the main
            % diagonal, K > 0 above the main diagonal, and K < 0 below the main diagonal. When |K| exceeds the
            % number of rows or columns in A, the function returns an empty rank-1 array.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            memory_obj = memory_mod();

            % Inputs


            % Outputs
            D = NaN(1);
            % Local variables
            srname = "DIAG";
            dlen = NaN;
            i = NaN;
            k_loc = NaN;

            %====================%
            % Calculation starts %
            %====================%

            ipObj = inputParser();
            addParameter(ipObj, 'k', 0);
            parse(ipObj, varargin{:});
            k_loc = ipObj.Results.k;


            % DLEN is the length of D. We allow |K| to exceed the number of rows/columns in A.
            dlen = max(0, fix(min(size(A, 1), size(A, 2)) - abs(k_loc)));
            D = memory_obj.alloc_rvector_sp(D, dlen);
            if k_loc >= 0
                D = reshape(cell2mat(arrayfun(@(i) A(i, i + k_loc), (1:dlen), "UniformOutput", false)), [], 1);
            else
                D = reshape(cell2mat(arrayfun(@(i) A(i - k_loc, i), (1:dlen), "UniformOutput", false)), [], 1);
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(D) == dlen, "SIZE(D) == DLEN", srname);
            end
        end
        function is_banded = isbanded(~, A, lwidth, uwidth, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This function tests whether the matrix A banded within the bandwidth specified by LWIDTH and
            % UWIDTH up to the tolerance TOL.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();

            % Inputs




            % Outputs
            is_banded = false;
            % Local variables
            srname = "ISBANDED";
            i = NaN;
            m = NaN;
            n = NaN;
            tol_loc = NaN;

            % Preconditions
            ipObj = inputParser();
            addParameter(ipObj, 'tol', NaN);
            parse(ipObj, varargin{:});
            tol = ipObj.Results.tol;
            if consts_obj.DEBUGGING
                debug_obj.assert(lwidth >= 0 && uwidth >= 0, "LWIDTH >= 0 .and. UWIDTH >= 0", srname);
                if ~ismember('tol', ipObj.UsingDefaults)
                    debug_obj.assert(tol >= 0, "TOL >= 0", srname);
                end
            end

            %====================%
            % Calculation starts %
            %====================%

            tol_loc = consts_obj.ZERO;
            if ~ismember('tol', ipObj.UsingDefaults)
                tol_loc = max(tol, tol * max(abs(A), [], 'all'));
            end
            if infnan_obj.is_nan_sp(tol_loc)
                tol_loc = consts_obj.ZERO;
            end

            m = size(A, 1);
            n = size(A, 2);

            is_banded = true;
            for i = 1:n
                is_banded = (all(abs(A(i + lwidth + 1:m, i)) <= tol_loc, 'all') && all(abs(A(1:i - uwidth - 1, i)) <= tol_loc, 'all'));
                if ~is_banded
                    break
                end
            end

            %====================%
            %  Calculation ends  %
            %====================%
        end
        function is_tril = istril(obj, A, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This function tests whether the matrix A is lower triangular up to the tolerance TOL.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();

            % Inputs


            % Outputs
            is_tril = false;
            % Local variables
            srname = "ISTRIL";
            width = NaN;
            tol_loc = NaN;

            % Preconditions
            ipObj = inputParser();
            addParameter(ipObj, 'tol', NaN);
            parse(ipObj, varargin{:});
            tol = ipObj.Results.tol;
            if consts_obj.DEBUGGING
                if ~ismember('tol', ipObj.UsingDefaults)
                    debug_obj.assert(tol >= 0, "TOL >= 0", srname);
                end
            end

            %====================%
            % Calculation starts %
            %====================%

            if ismember('tol', ipObj.UsingDefaults)
                tol_loc = consts_obj.ZERO;
            else
                tol_loc = tol;
            end
            width = fix(max(0, size(A, 1) - 1));
            is_tril = obj.isbanded(A, width, 0, 'tol', tol_loc);

            %====================%
            %  Calculation ends  %
            %====================%
        end
        function is_triu = istriu(obj, A, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This function tests whether the matrix A is upper triangular up to the tolerance TOL.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();

            % Inputs


            % Outputs
            is_triu = false;
            % Local variables
            srname = "ISTRIU";
            width = NaN;
            tol_loc = NaN;

            % Preconditions
            ipObj = inputParser();
            addParameter(ipObj, 'tol', NaN);
            parse(ipObj, varargin{:});
            tol = ipObj.Results.tol;
            if consts_obj.DEBUGGING
                if ~ismember('tol', ipObj.UsingDefaults)
                    debug_obj.assert(tol >= 0, "TOL >= 0", srname);
                end
            end

            %====================%
            % Calculation starts %
            %====================%

            if ismember('tol', ipObj.UsingDefaults)
                tol_loc = consts_obj.ZERO;
            else
                tol_loc = tol;
            end
            width = fix(max(0, size(A, 2) - 1));
            is_triu = obj.isbanded(A, 0, width, 'tol', tol_loc);

            %====================%
            %  Calculation ends  %
            %====================%
        end
        function is_orth = isorth(obj, A, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This function tests whether the matrix A has orthonormal columns up to the tolerance TOL.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();

            % Inputs


            % Outputs
            is_orth = false;
            % Local variables
            srname = "ISORTH";
            n = NaN;
            tol_loc = NaN;

            % Preconditions
            ipObj = inputParser();
            addParameter(ipObj, 'tol', NaN);
            parse(ipObj, varargin{:});
            tol = ipObj.Results.tol;
            if consts_obj.DEBUGGING
                if ~ismember('tol', ipObj.UsingDefaults)
                    debug_obj.assert(tol >= 0, "TOL >= 0", srname);
                end
            end

            %====================%
            % Calculation starts %
            %====================%

            tol_loc = consts_obj.ORTHTOL_DFT;
            if ~ismember('tol', ipObj.UsingDefaults)
                tol_loc = tol;
            end

            n = size(A, 2);

            % N.B. (20240304): In some cases, due to compiler bugs, we need to disable the test. We signify such
            % cases by setting ORTHTOL_DFT to REALMAX. For instance, NAG Fortran Compiler Release 7.1(Hanzomon)
            % Build 7143 is buggy concerning half-precision floating-point numbers. See the following:
            % https://fortran-lang.discourse.group/t/nagfor-7-1-supports-half-precision-floating-point-numbers-but-with-many-bugs
            is_orth = true;
            if n > size(A, 1)
                is_orth = false;
            elseif any(infnan_obj.is_nan(A), 'all')
                is_orth = false;
            elseif consts_obj.ORTHTOL_DFT < consts_obj.REALMAX
                is_orth = all(abs(obj.matprod22(A', A) - obj.eye1(n)) <= max(tol_loc, tol_loc * max(abs(A), [], 'all')), 'all');
            end

            %====================%
            %  Calculation ends  %
            %====================%
        end
        function y = project1(obj, x, v)
            %--------------------------------------------------------------------------------------------------%
            % This function returns the projection of X to SPAN(V).
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();

            % Inputs


            % Outputs
            y = NaN(numel(x), 1);
            % Local variables
            srname = "PROJECT1";
            u = NaN(numel(v), 1);
            tol = NaN;

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(x) == numel(v), "SIZE(X) == SIZE(V)", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            if all(abs(x) <= 0, 'all') || all(abs(v) <= 0, 'all')
                y = repmat(consts_obj.ZERO, size(y));
            elseif any(infnan_obj.is_nan(x), 'all') || any(infnan_obj.is_nan(v), 'all')
                y = repmat(sum(x, 'all') + sum(v, 'all'), size(y)); % Set Y to NaN

            elseif any(infnan_obj.is_inf(v), 'all')
                u = repmat(consts_obj.ZERO, size(u));
                u(obj.trueloc(infnan_obj.is_inf(v))) = abs(consts_obj.ONE) .* ((v(obj.trueloc(infnan_obj.is_inf(v))) >= 0) * 2 - 1);
                %%MATLAB: u = 0; u(isinf(v)) = sign(v(isinf(v)))
                u(:) = u ./ obj.p_norm(u);
                y(:) = obj.inprod(x, u) * u;
            else
                u(:) = v ./ obj.p_norm(v);
                y(:) = obj.inprod(x, u) * u;
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                if infnan_obj.is_finite(obj.p_norm(x)) && infnan_obj.is_finite(obj.p_norm(v))
                    tol = max(consts_obj.TEN ^ max(-10, -consts_obj.MAXPOW10), min(0.1, consts_obj.TEN ^ min(6, consts_obj.MAXPOW10) * consts_obj.EPS));
                    debug_obj.assert(obj.p_norm(y) <= (consts_obj.ONE + tol) * obj.p_norm(x), "NORM(Y) <= NORM(X)", srname);
                    debug_obj.assert(obj.p_norm(x - y) <= (consts_obj.ONE + tol) * obj.p_norm(x), "NORM(X - Y) <= NORM(X)", srname);
                    % The following test may not be passed.
                    debug_obj.assert(abs(obj.inprod(x - y, v)) <= max(tol, tol * max(obj.p_norm(x - y) * obj.p_norm(v), abs(obj.inprod(x, v)))), "X - Y is orthogonal to V", srname);
                end
            end
        end
        function y = project2(obj, x, V)
            %--------------------------------------------------------------------------------------------------%
            % This function returns the projection of X to RANGE(V).
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();

            % Inputs


            % Outputs
            y = NaN(numel(x), 1);
            % Local variables
            srname = "PROJECT2";
            U = NaN(size(V, 1), min(size(V, 1), size(V, 2)));
            V_loc = NaN(size(V, 1), size(V, 2));
            tol = NaN;

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(x) == size(V, 1), "SIZE(X) == SIZE(V, 1)", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            if size(V, 2) == 1
                y(:) = obj.project1(x, V(:, 1));
            elseif all(abs(x) <= 0, 'all') || all(abs(V) <= 0, 'all')
                y = repmat(consts_obj.ZERO, size(y));
            elseif any(infnan_obj.is_nan(x), 'all') || any(infnan_obj.is_nan(V), 'all')
                y = repmat(sum(x, 'all') + sum(V, 'all'), size(y)); % Set Y to NaN

            elseif any(infnan_obj.is_inf(V), 'all')
                mask00 = infnan_obj.is_inf(V); %Unsupported statement inside WHERE block: StmtLineBreak 1
                V_loc(mask00) = abs(consts_obj.ONE) .* ((V(mask00) >= 0) * 2 - 1); %Unsupported statement inside WHERE block: StmtLineBreak 1
                mask01 = ~mask00; %Unsupported statement inside WHERE block: StmtLineBreak 1
                V_loc(mask01) = consts_obj.ZERO; %Unsupported statement inside WHERE block: StmtLineBreak 1

                %%MATLAB: V_loc = 0; V_loc(isinf(V)) = sign(V);
                U = obj.qr(V_loc);
                y(:) = obj.matprod21(U, obj.matprod12(x, U));
            else
                U = obj.qr(V);
                y(:) = obj.matprod21(U, obj.matprod12(x, U));
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                if infnan_obj.is_finite(obj.p_norm(x)) && infnan_obj.is_finite(sum(V .^ 2, 'all'))
                    tol = max(consts_obj.TEN ^ max(-10, -consts_obj.MAXPOW10), min(0.1, consts_obj.TEN ^ min(6, consts_obj.MAXPOW10) * consts_obj.EPS));
                    debug_obj.assert(obj.p_norm(y) <= (consts_obj.ONE + tol) * obj.p_norm(x), "NORM(Y) <= NORM(X)", srname);
                    debug_obj.assert(obj.p_norm(x - y) <= (consts_obj.ONE + tol) * obj.p_norm(x), "NORM(X - Y) <= NORM(X)", srname);
                    % The following test may not be passed.
                    debug_obj.assert(obj.p_norm(obj.matprod12(x - y, V)) <= max(tol, tol * max(obj.p_norm(x - y) * obj.named_norm_mat(V, "fro"), obj.p_norm(obj.matprod12(x, V)))), "X - Y is orthogonal to V", srname);
                end
            end
        end
        function r = hypotenuse(~, x1, x2)
            %--------------------------------------------------------------------------------------------------%
            % HYPOTENUSE(X1, X2) returns SQRT(X1^2 + X2^2), handling over/underflow.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();

            % Inputs


            % Outputs
            r = NaN;
            % Local variables
            srname = "HYPOTENUSE";
            y = NaN(2, 1);

            %====================%
            % Calculation starts %
            %====================%

            if ~infnan_obj.is_finite(x1)
                r = abs(x1);
            elseif ~infnan_obj.is_finite(x2)
                r = abs(x2);
            else
                y(:) = abs([x1, x2]);
                y(:) = [min(y, [], 'all'), max(y, [], 'all')];
                if y(1) > sqrt(consts_obj.REALMIN) && y(2) < sqrt(consts_obj.REALMAX / 2.1)
                    r = sqrt(sum(y .^ 2, 'all'));
                elseif y(2) > 0
                    r = y(2) * sqrt((y(1) / y(2)) ^ 2 + consts_obj.ONE);
                else
                    r = consts_obj.ZERO;
                end
                % Without the following line, R > Y(1) + Y(2) or R < Y(2) may happen due to rounding errors.
                r = min(sum(y, 'all'), max(y(2), r));
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                if infnan_obj.is_nan_sp(x1) || infnan_obj.is_nan_sp(x2)
                    debug_obj.assert(infnan_obj.is_nan_sp(r), "R is NaN if X1 or X2 is NaN", srname);
                else
                    debug_obj.assert(r >= abs(x1) && r >= abs(x2) && r <= abs(x1) + abs(x2), "MAX{ABS(X1), ABS(X2)} <= R <= ABS(X1) + ABS(X2)", srname);
                end
            end
        end
        function G = planerot(obj, x)
            %--------------------------------------------------------------------------------------------------%
            % As in MATLAB, PLANEROT(X) returns a 2x2 Givens matrix G for X in R^2 so that Y = G*X has Y(2) = 0.
            % Roughly speaking, using a MATLAB-style formulation of matrices,
            % G = [X(1)/R, X(2)/R; -X(2)/R, X(1)/R] with R = SQRT(X(1)^2+X(2)^2), and G*X = [R; 0].
            % 0. We need to take care of the possibilities of R = 0, Inf, NaN, and over/underflow.
            % 1. The G defined above is continuous with respect to X except at 0. Following this definition,
            % G = [sign(X(1)), 0; 0, sign(X(1))] if X(2) = 0, G = [0, sign(X(2)); -sign(X(2)), 0] if X(2) = 0.
            % Yet some implementations ignore the signs, leading to discontinuity and numerical instability.
            % 2. Difference from MATLAB: if X contains NaN or consists of only Inf, MATLAB returns a NaN matrix,
            % but we return an identity matrix or a matrix of +/-SQRT(2). We intend to keep G always orthogonal.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();

            % Inputs

            % Outputs
            G = NaN(2);
            % Local variables
            srname = "PLANEROT";
            c = NaN;
            s = NaN;
            r = NaN;
            t = NaN;
            u = NaN;
            tol = NaN;

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(x) == 2, "SIZE(X) == 2", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            % Define C = X(1) / R and S = X(2) / R with R = HYPOT(X(1), X(2)). Handle Inf/NaN, over/underflow.
            if any(infnan_obj.is_nan(x), 'all')
                % In this case, MATLAB sets G to NaN(2, 2). We refrain from doing so to keep G orthogonal.
                c = consts_obj.ONE;
                s = consts_obj.ZERO;
            elseif all(infnan_obj.is_inf(x), 'all')
                % In this case, MATLAB sets G to NaN(2, 2). We refrain from doing so to keep G orthogonal.
                c = abs(1 / sqrt(2.0)) .* ((x(1) >= 0) * 2 - 1);
                s = abs(1 / sqrt(2.0)) .* ((x(2) >= 0) * 2 - 1);
            elseif abs(x(1)) <= 0 && abs(x(2)) <= 0
                % X(1) == 0 == X(2).
                c = consts_obj.ONE;
                s = consts_obj.ZERO;
            elseif abs(x(2)) <= consts_obj.EPS * abs(x(1))
                % N.B.:
                % 0. With <= instead of <, this case covers X(1) == 0 == X(2), which is treated above separately
                % to avoid the confusing SIGN(., 0) (see 1).
                % 1. SIGN(A, 0) = ABS(A) in Fortran but sign(0) = 0 in MATLAB, Python, Julia, and R!
                % 2. Taking SIGN(X(1)) into account ensures the continuity of G with respect to X except at 0.
                c = abs(consts_obj.ONE) .* ((x(1) >= 0) * 2 - 1); %%MATLAB: c = sign(x(1))
                s = consts_obj.ZERO;
            elseif abs(x(1)) <= consts_obj.EPS * abs(x(2))
                % N.B.: SIGN(A, X) = ABS(A) * sign of X /= A * sign of X ! Therefore, it is WRONG to define G
                % as SIGN(RESHAPE([ZERO, -ONE, ONE, ZERO], [2, 2]), X(2)). This mistake was committed on
                % 20211206 and took a whole day to debug! NEVER use SIGN on arrays unless you are really sure.
                c = consts_obj.ZERO;
                s = abs(consts_obj.ONE) .* ((x(2) >= 0) * 2 - 1); %%MATLAB: s = sign(x(2))

            else
                % Here is the normal case. It implements the Givens rotation in a stable & continuous way as in:
                % Bindel, D., Demmel, J., Kahan, W., and Marques, O. (2002). On computing Givens rotations
                % reliably and efficiently. ACM Transactions on Mathematical Software (TOMS), 28(2), 206-238.
                % N.B.: 1. Modern compilers compute SQRT(REALMIN) and SQRT(REALMAX/2.1) at compilation time.
                % 2. The direct calculation without involving T and U seems to work better; use it if possible.
                if all(abs(x) > sqrt(consts_obj.REALMIN) & abs(x) < sqrt(consts_obj.REALMAX / 2.1), 'all')
                    % Do NOT use HYPOTENUSE here; the best implementation for one may be suboptimal for the other
                    r = obj.p_norm(x);
                    c = x(1) / r;
                    s = x(2) / r;
                elseif abs(x(1)) > abs(x(2))
                    t = x(2) / x(1);
                    u = max([consts_obj.ONE, abs(t), sqrt(consts_obj.ONE + t ^ 2)], [], 'all'); % MAXVAL: precaution against rounding error.
                    u = abs(u) .* ((x(1) >= 0) * 2 - 1); %%MATLAB: u = sign(x(1))*sqrt(ONE + t**2)
                    c = consts_obj.ONE / u;
                    s = t / u;
                else
                    t = x(1) / x(2);
                    u = max([consts_obj.ONE, abs(t), sqrt(consts_obj.ONE + t ^ 2)], [], 'all'); % MAXVAL: precaution against rounding error.
                    u = abs(u) .* ((x(2) >= 0) * 2 - 1); %%MATLAB: u = sign(x(2))*sqrt(ONE + t**2)
                    c = t / u;
                    s = consts_obj.ONE / u;
                end
            end

            G(:, :) = reshape([c, -s, s, c], [2, 2]); %%MATLAB: G = [c, s; -s, c]

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(size(G, 1) == 2 && size(G, 2) == 2, "SIZE(G) == [2, 2]", srname);
                debug_obj.assert(all(infnan_obj.is_finite(G), 'all'), "G is finite", srname);
                debug_obj.assert(abs(G(1, 1) - G(2, 2)) + abs(G(1, 2) + G(2, 1)) <= 0, "G(1,1) == G(2,2), G(1,2) = -G(2,1)", srname);
                tol = max(consts_obj.TEN ^ max(-10, -consts_obj.MAXPOW10), min(0.1, 10.0 ^ min(6, consts_obj.MAXPOW10) * consts_obj.EPS));
                debug_obj.assert(obj.isorth(G, 'tol', tol), "G is orthonormal", srname);
                if all(infnan_obj.is_finite(x) & abs(x) < sqrt(consts_obj.REALMAX / 2.1), 'all')
                    r = obj.p_norm(x);
                    debug_obj.assert(max(abs(obj.matprod21(G, x) - [r, consts_obj.ZERO]), [], 'all') <= max(tol, tol * r), "G * X = [||X||, 0]", srname);
                end
            end
        end
        function A = symmetrize(obj, A)
            %--------------------------------------------------------------------------------------------------%
            % SYMMETRIZE(A) symmetrizes A.
            % N.B.: Here, we assume that A is a matrix that IS SUPPOSED TO BE symmetric in precise arithmetic,
            % and its asymmetry comes only from errors (e.g., rounding, noise).
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();

            % In-outputs

            % Local variables
            j = NaN;
            srname = "SYMMETRIZE";

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(size(A, 1) == size(A, 2), "A is square", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            % A is symmetrized by copying A(LOWER_TRI) to A(UPPER_TRI).
            % N.B.: The following assumes that A(LOWER_TRI) has been properly defined.
            for j = 1:size(A, 1)
                A(1:j - 1, j) = A(j, 1:j - 1);
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(obj.issymmetric(A), "A is symmetrized", srname);
            end
        end
        function is_minor = isminor0(~, x, ref)
            %--------------------------------------------------------------------------------------------------%
            % This function tests whether X is minor compared to REF. It is used by Powell, e.g., in COBYLA.
            % In precise arithmetic, ISMINOR(X, REF) is TRUE if and only if X == 0; in floating-point
            % arithmetic, ISMINOR(X, REF) is true if X is zero or its nonzero value can be attributed to
            % computer rounding errors according to REF.
            % Larger SENSITIVITY means the function is more strict/precise, the value TENTH being due to Powell.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();

            % Inputs


            % Outputs
            is_minor = false;
            % Local variables
            sensitivity = consts_obj.TENTH;
            refa = NaN;
            refb = NaN;

            %====================%
            % Calculation starts %
            %====================%

            refa = abs(ref) + sensitivity * abs(x);
            refb = abs(ref) + consts_obj.TWO * sensitivity * abs(x);
            is_minor = (abs(ref) >= refa || refa >= refb);

            %====================%
            %  Calculation ends  %
            %====================%

        end
        function is_minor = isminor1(obj, x, ref)
            %--------------------------------------------------------------------------------------------------%
            % This function tests whether X is minor compared to REF. It is used by Powell, e.g., in COBYLA.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();

            % Inputs


            % Outputs
            is_minor = false(numel(x), 1);
            % Local variables
            srname = "ISMINOR1";
            i = NaN;

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(x) == numel(ref), "SIZE(X) == SIZE(REF)", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            is_minor(:) = reshape(cell2mat(arrayfun(@(i) obj.isminor0(x(i), ref(i)), (1:fix(numel(x))), "UniformOutput", false)), [], 1);

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(is_minor) == numel(x), "SIZE(IS_MINOR) == SIZE(X)", srname);
            end
        end
        function is_symmetric = issymmetric(~, A, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This function tests whether A is symmetric up to TOL.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();

            % Inputs


            % Outputs
            is_symmetric = false;
            % Local variables
            srname = "ISSYMMETRIC";
            tol_loc = NaN;

            % Preconditions
            ipObj = inputParser();
            addParameter(ipObj, 'tol', NaN);
            parse(ipObj, varargin{:});
            tol = ipObj.Results.tol;
            if consts_obj.DEBUGGING
                if ~ismember('tol', ipObj.UsingDefaults)
                    debug_obj.assert(tol >= 0, "TOL >= 0", srname);
                end
            end

            %====================%
            % Calculation starts %
            %====================%

            tol_loc = consts_obj.SYMTOL_DFT;
            if ~ismember('tol', ipObj.UsingDefaults)
                tol_loc = tol;
            end

            % N.B.:
            % 0. It may be expensive to take TRANSPOSE(A), let alone doing it multiple times, but this is not an
            % issue in our project. We call ISSYMMETRIC only in the debugging mode, but never in production.
            % 1. In Fortran, the following instructions cannot be written as the following Boolean expression:
            % %IS_SYMMETRIC = (SIZE(A, 1)==SIZE(A, 2) .AND. &
            % % & ALL(IS_NAN(A) .EQV. IS_NAN(TRANSPOSE(A))) .AND. &
            % % & .NOT. ANY(ABS(A - TRANSPOSE(A)) > TOL_LOC * MAX(MAXVAL(ABS(A)), ONE)))
            % This is because Fortran may not perform short-circuit evaluation of this expression. If A is not
            % square, then IS_NAN(A) .EQV. IS_NAN(TRANSPOSE(A)) and A - TRANSPOSE(A) are invalid.
            % 2. In addition, since Inf - Inf is NaN, we cannot replace ANY(ABS(A - TRANSPOSE(A)) > TOL_LOC ...)
            % with .NOT. ALL(ABS(A - TRANSPOSE(A)) <= TOL_LOC ...).
            % 3. In some cases, due to compiler bugs / features, we need to disable the test. We signify such
            % cases by setting SYMTOL_DFT to REALMAX. For instance, when invoked with aggressive optimization
            % options (e.g., -fast-math), gfortran 11 is buggy with ALL and ANY: ALL returns .FALSE. on a vector
            % of .TRUE., while ANY returns .TRUE. on a vector of .FALSE.. In that case, we cannot test
            % ALL(IS_NAN(A) .EQV. IS_NAN(TRANSPOSE(A))).
            is_symmetric = true;
            if size(A, 1) ~= size(A, 2)
                is_symmetric = false;
            elseif consts_obj.SYMTOL_DFT < 0.9 * consts_obj.REALMAX
                is_symmetric = (~any(abs(A - A') > tol_loc * max(max(abs(A), [], 'all'), consts_obj.ONE), 'all')) && all(infnan_obj.is_nan(A) == infnan_obj.is_nan(A'), 'all');
            end

            %====================%
            %  Calculation ends  %
            %====================%
        end
        function y = p_norm(~, x, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This function calculates the P-norm of a vector X.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();

            % Inputs


            % Outputs
            y = NaN;
            % Local variables
            srname = "P_NORM";
            maxabs = NaN;
            p_loc = NaN;
            scaling = NaN;
            scalmax = NaN;
            scalmin = NaN;

            % Preconditions
            ipObj = inputParser();
            addParameter(ipObj, 'p', NaN);
            parse(ipObj, varargin{:});
            p = ipObj.Results.p;
            if ~ismember('p', ipObj.UsingDefaults)
                debug_obj.validate(p >= 0, "P >= 0", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            if ismember('p', ipObj.UsingDefaults)
                p_loc = consts_obj.TWO;
            else
                p_loc = p;
            end

            % If SIZE(X) = 0, then MAXVAL(ABS(X)) = -HUGE(X); since we handle such a case individually,
            % it is OK to write MAXVAL(ABS(X)) below, but we append 0 for robustness.
            maxabs = max([reshape(abs(x), [], 1); consts_obj.ZERO], [], 'all');

            if numel(x) == 0
                y = consts_obj.ZERO;
            elseif p_loc <= 0 && ~any(infnan_obj.is_nan(x), 'all')
                y = double(nnz(abs(x) > 0));
            elseif ~all(infnan_obj.is_finite(x), 'all')
                % If X contains NaN, then Y is NaN. Otherwise, Y is Inf when X contains +/-Inf unless P = 0.
                y = sum(abs(x), 'all');
            elseif maxabs <= 0
                % If MAXABS is zero, then Y is zero. Note that we do this only when X does not contain NaN.
                % Otherwise, MAXABS = 0 does not necessarily guarantee that X is all zero.
                y = consts_obj.ZERO;
            else                % Now P > 0 and X is a finite-valued nonzero vector, as we have handled the other cases above.
                if infnan_obj.is_posinf(p_loc)
                    y = maxabs;
                elseif abs(p_loc - consts_obj.ONE) <= 0
                    y = sum(abs(x), 'all');
                elseif abs(p_loc - consts_obj.TWO) <= 0
                    % N.B.: We may use the intrinsic NORM2. Here, we use the following naive implementation
                    % to get full control on the computation in a way similar to MATPROD and INPROD.

                    % To avoid over/underflow, we scale X by SCALING defined as follows if necessary. We make
                    % sure SCALMIN >= MAX(REALMIN, 1/REALMAX) and SCALMAX <= MIN(REALMAX, 1/REALMIN), or we may
                    % encounter over/underflow when dividing by SCALING, and even NaN if the compiler evaluates
                    % 1/SCALING first, which did happen with `flang -ffast-math` on REAL32 with LLVM flang 21.
                    %
                    % Given a numeric model for floating-point numbers,
                    % REALMIN = TINY(ZERO) = 2^{emin-1}, REALMAX = HUGE(ZERO) = (1-b^{-d})*b^{emax} >= b^{emax-1},
                    % where b = RADIX(ZERO) is the base, d = DIGITS(ZERO) > 0 is the number of base-b significant
                    % digits, emin = MINEXPONENT(ZERO) & emax = MAXEXPONENT(ZERO) are the min & max exponents.
                    % N.B.: IEEE 754 specifies emax and requires that emin = 1 - emax for "Binary interchange
                    % floating-point formats" binary32, binary64, and binary128 (Sec. 3.3 of IEEE Std 754-2019).
                    % However, mathematically, [emin, emax] defined in Fortran standards indeed corresponds to
                    % [emin+1, emax+1] in IEEE 754. In addition, Fortran compilers may not implement REAL32,
                    % REAL64, and REAL128 according to binary32, binary64, and binary128. For instance,
                    % nagfor 7 has d = 106, emin = -968 and emax = 1023 for REAL128, while
                    % IEEE 754 has d = 113, emin = -16382, and emax = 16383 for binary128. See
                    % http://fortran-lang.discourse.group/t/ieee-754-binary-interchange-floating-point-formats-versus-iso-fortran-env-real-kinds

                    y = sqrt(sum(x .^ 2, 'all'));
                    % The following code handles over/underflow naively.
                    if infnan_obj.is_posinf(y) || y <= 0
                        scalmin = double(radix(consts_obj.ZERO)) ^ max(minexponent(consts_obj.ZERO) - 1, 1 - maxexponent(consts_obj.ZERO));
                        scalmax = double(radix(consts_obj.ZERO)) ^ min(maxexponent(consts_obj.ZERO) - 1, 1 - minexponent(consts_obj.ZERO));
                        scaling = min(max(maxabs, scalmin), scalmax);
                        y = scaling * sqrt(sum((x ./ scaling) .^ 2, 'all'));
                    end
                else
                    y = sum(abs(x) .^ p_loc, 'all') ^ (consts_obj.ONE / p_loc);
                    % The following code handles over/underflow naively.
                    if infnan_obj.is_posinf(y) || y <= 0
                        scalmin = double(radix(consts_obj.ZERO)) ^ max(minexponent(consts_obj.ZERO) - 1, 1 - maxexponent(consts_obj.ZERO));
                        scalmax = double(radix(consts_obj.ZERO)) ^ min(maxexponent(consts_obj.ZERO) - 1, 1 - minexponent(consts_obj.ZERO));
                        scaling = min(max(maxabs, scalmin), scalmax);
                        y = scaling * sum(abs(x ./ scaling) .^ p_loc, 'all') ^ (consts_obj.ONE / p_loc);
                    end
                end
            end

            %====================%
            %  Calculation ends  %
            %====================%

            if consts_obj.DEBUGGING
                debug_obj.assert(y >= 0 || any(infnan_obj.is_nan(x), 'all'), "Y >= 0 unless X contains NaN", srname);
                debug_obj.assert(infnan_obj.is_nan_sp(y) == any(infnan_obj.is_nan(x), 'all'), "Y is NaN if and only if X contains NaN", srname);
                % Even with scaling, Y may still be 0 if all entries of X are zero or subnormal.
                debug_obj.assert(y > 0 || any(infnan_obj.is_nan(x), 'all') || all(abs(x) < consts_obj.REALMIN, 'all'), "Y > 0 unless X contains NaN or all its entries are below REALMIN", srname);
            end

        end
        function y = named_norm_vec(obj, x, nname)
            %--------------------------------------------------------------------------------------------------%
            % This function calculates named norms of a vector X.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();
            string_obj = string_mod();

            % Inputs


            % Outputs
            y = NaN;
            % Local variables
            srname = "NAMED_NORM_VEC";

            %====================%
            % Calculation starts %
            %====================%

            if numel(x) == 0
                y = consts_obj.ZERO;
            elseif ~all(infnan_obj.is_finite(x), 'all')
                % If X contains NaN, then Y is NaN. Otherwise, Y is Inf when X contains +/-Inf.
                y = sum(abs(x), 'all');
            elseif ~any(abs(x) > 0, 'all')
                % The following is incorrect without checking the last case, as X may be all NaN.
                y = consts_obj.ZERO;
            else
                switch string_obj.lower(string_obj.strip(nname))
                case "fro"
                    y = obj.p_norm(x); % 2-norm, which is the default case of P_NORM.
                case "inf"
                    % If SIZE(X) = 0, then MAXVAL(ABS(X)) = -HUGE(X); since we have handled such a case in the
                    % above, it is OK to write Y = MAXVAL(ABS(X)) below, but we append a 0 for robustness.
                    y = max([reshape(abs(x), [], 1); consts_obj.ZERO], [], 'all');
                otherwise
                    debug_obj.warning(srname, "Unknown name of norm: " + string_obj.strip(nname) + "; default to the L2-norm");
                    y = obj.p_norm(x); % 2-norm, which is the default case of P_NORM.
                end
            end

            %====================%
            %  Calculation ends  %
            %====================%
        end
        function y = named_norm_mat(~, x, nname)
            %--------------------------------------------------------------------------------------------------%
            % This function calculates named norms of a vector X.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();
            string_obj = string_mod();

            % Inputs


            % Outputs
            y = NaN;
            % Local variables
            srname = "NAMED_NORM_MAT";

            %====================%
            % Calculation starts %
            %====================%

            % N.B.: Ideally, we should also do a scaling similar to that in P_NORM to avoid over/underflow.

            if size(x, 1) * size(x, 2) == 0
                y = consts_obj.ZERO;
            elseif ~all(infnan_obj.is_finite(x), 'all')
                % If X contains NaN, then Y is NaN. Otherwise, Y is Inf when X contains +/-Inf.
                y = sum(abs(x), 'all');
            elseif ~any(abs(x) > 0, 'all')
                % The following is incorrect without checking the last case, as X may be all NaN.
                y = consts_obj.ZERO;
            else
                switch string_obj.lower(string_obj.strip(nname))
                case "fro"
                    y = sqrt(sum(x .^ 2, 'all'));
                case "inf"
                    % If SIZE(X) = 0, then MAXVAL(SUM(ABS(X), DIM=2)) = -HUGE(X); since we have handled such a
                    % case in the above, it is OK to write Y = MAXVAL(SUM(ABS(X), DIM=2)) below, but we append
                    % a 0 for robustness.
                    y = max([reshape(sum(abs(x), 2), 1, []), consts_obj.ZERO], [], 'all');
                otherwise
                    debug_obj.warning(srname, "Unknown name of norm: " + string_obj.strip(nname) + "; default to the Frobenius norm");
                    y = sqrt(sum(x .^ 2, 'all'));
                end
            end

            %====================%
            %  Calculation ends  %
            %====================%
        end
        function y = sort_i1(~, x, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This function sorts X according to DIRECTION, which should be 'ascend' (default) or 'descend'.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();

            % Inputs


            % Outputs
            y = NaN(numel(x), 1);
            % Local variables
            srname = "SORT_I1";

            i = NaN;
            n = NaN;
            newn = NaN;
            ascending = false;

            %====================%
            % Calculation starts %
            %====================%

            ascending = true;
            ipObj = inputParser();
            addParameter(ipObj, 'direction', "");
            parse(ipObj, varargin{:});
            direction = ipObj.Results.direction;
            if ~ismember('direction', ipObj.UsingDefaults)
                if direction == "descend" || direction == "DESCEND"
                    ascending = false;
                end
            end

            y(:) = x;
            n = fix(numel(y));
            while n > 1                % Bubble sort.
                newn = 0;
                for i = 2:n
                    if (y(i - 1) > y(i) && ascending) || (y(i - 1) < y(i) && ~ascending)
                        y([i - 1, i]) = y([i, i - 1]);
                        newn = i;
                    end
                end
                n = newn;
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                if ascending
                    debug_obj.assert(all(y(1:n - 1) <= y(2:n), 'all'), "Y is ascending", srname);
                else
                    debug_obj.assert(all(y(1:n - 1) >= y(2:n), 'all'), "Y is descending", srname);
                end
            end
        end
        function y = sort_i2(obj, x, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This function sorts a matrix X according to DIM (1 or 2) and DIRECTION ('ascend' or 'descend').
            %--------------------------------------------------------------------------------------------------%
            debug_obj = debug_mod();
            consts_obj = consts_mod();
            string_obj = string_mod();

            % Inputs



            % Outputs
            y = NaN(size(x, 1), size(x, 2));
            % Local variables
            srname = "SORT_I2";
            direction_loc = "";
            dim_loc = NaN;
            i = NaN;
            n = NaN;

            %====================%
            % Calculation starts %
            %====================%

            dim_loc = 1;
            ipObj = inputParser();
            addParameter(ipObj, 'dim', NaN);
            addParameter(ipObj, 'direction', "");
            parse(ipObj, varargin{:});
            dim = ipObj.Results.dim;
            direction = ipObj.Results.direction;
            if ~ismember('dim', ipObj.UsingDefaults)
                dim_loc = dim;
            end

            direction_loc = "ascend";
            if ~ismember('direction', ipObj.UsingDefaults)
                direction_loc = string_obj.strip(direction);
            end

            y(:, :) = x;
            if dim_loc == 1
                for i = 1:size(x, 2)
                    y(:, i) = obj.sort_i1(y(:, i), 'direction', direction_loc);
                end
            else
                for i = 1:size(x, 1)
                    y(i, :) = obj.sort_i1(y(i, :), 'direction', direction_loc);
                end
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                if dim_loc == 1
                    n = size(y, 1);
                    if direction_loc == "ascend" || direction_loc == "ASCEND"
                        debug_obj.assert(all(y(1:n - 1, :) <= y(2:n, :), 'all'), "Y is ascending along dimension 1", srname);
                    else
                        debug_obj.assert(all(y(1:n - 1, :) >= y(2:n, :), 'all'), "Y is descending along dimension 1", srname);
                    end
                else
                    n = size(y, 2);
                    if direction_loc == "ascend" || direction_loc == "ASCEND"
                        debug_obj.assert(all(y(:, 1:n - 1) <= y(:, 2:n), 'all'), "Y is ascending along dimension 2", srname);
                    else
                        debug_obj.assert(all(y(:, 1:n - 1) >= y(:, 2:n), 'all'), "Y is descending along dimension 2", srname);
                    end
                end
            end
        end
        function y = logical_to_int(~, x)
            %--------------------------------------------------------------------------------------------------%
            % LOGICAL_TO_INT(.TRUE.) = 1, LOGICAL_TO_INT(.FALSE.) = 0
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();

            % Inputs

            % Outputs


            y = fortran.merge('tsource', 1, 'fsource', 0, 'mask', x);
        end
        function loc = trueloc(obj, x)
            %--------------------------------------------------------------------------------------------------%
            % Similar to the `find` function in MATLAB, TRUELOC returns the indices where X is true in
            % the ASCENDING order.
            % The motivation for this function is the fact that Fortran does not support logical indexing. See,
            % for example, https:
            % 1. MATLAB, Python, Julia, and R support logical indexing, so that the Fortran code Y(TRUELOC(X))
            % can simply be translated to Y(X).
            % 2. If the return of TRUELOC is NOT used for indexing, its analogs in other languages are:
            % MATLAB -- find, Python -- numpy.argwhere, Julia -- findall, R -- which.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            memory_obj = memory_mod();

            % Inputs

            % Outputs
            loc = NaN(1); % INTEGER(IK) :: LOC(COUNT(X)) does not work with Absoft 22.0
            % Local variables
            srname = "TRUELOC";
            n = NaN;

            %====================%
            % Calculation starts %
            %====================%

            loc = memory_obj.alloc_ivector(loc, fix(nnz(x))); % Removable in F03.
            n = fix(numel(x));
            loc = feval(@(a, m) reshape(a(m & true(size(a))), [], 1), obj.linspace_i(1, n, n), x);

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(all(loc >= 1 & loc <= n, 'all'), "1 <= LOC <= N", srname);
                debug_obj.assert(numel(loc) == nnz(x), "SIZE(LOC) == COUNT(X)", srname);
                debug_obj.assert(all(x(loc), 'all'), "X(LOC) is all TRUE", srname);
                debug_obj.assert(all(loc(2:numel(loc)) > loc(1:numel(loc) - 1), 'all'), "LOC is strictly ascending", srname);
            end
        end
        function loc = falseloc(obj, x)
            %--------------------------------------------------------------------------------------------------%
            % FALSELOC = TRUELOC(.NOT. X)
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            memory_obj = memory_mod();

            % Inputs

            % Outputs
            loc = NaN(1); % INTEGER(IK) :: LOC(COUNT(.NOT.X)) does not work with Absoft 22.0
            % Local variables
            srname = "FALSELOC";

            %====================%
            % Calculation starts %
            %====================%

            loc = memory_obj.alloc_ivector(loc, fix(nnz(~x))); % Removable in F03.
            loc = obj.trueloc(~x);

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(all(loc >= 1 & loc <= numel(x), 'all'), "1 <= LOC <= N", srname);
                debug_obj.assert(numel(loc) == numel(x) - nnz(x), "SIZE(LOC) == SIZE(X) - COUNT(X)", srname);
                debug_obj.assert(all(~x(loc), 'all'), "X(LOC) is all FALSE", srname);
                debug_obj.assert(all(loc(2:numel(loc)) > loc(1:numel(loc) - 1), 'all'), "LOC is strictly ascending", srname);
            end
        end
        function y = minimum1(~, x)
            %--------------------------------------------------------------------------------------------------%
            % This function returns NaN if X contains NaN; otherwise, it returns MINVAL(X). Vector version.
            % F2018 does not specify MINVAL(X) when X contains NaN, which motivates this function. The behavior
            % of this function is the same as the following functions in various languages:
            % MATLAB: min(x, [], 'includenan')
            % Python: numpy.min(x)
            % Julia: minimum(x)
            % R: min(x)
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();

            % Inputs

            % Outputs
            y = NaN;
            % Local variables
            srname = "MINIMUM1";
            nan_test = NaN;

            %====================%
            % Calculation starts %
            %====================%

            %y = merge(tsource=sum(x), fsource=minval(x), mask=any(is_nan(x)))
            nan_test = sum(abs(x), 'all'); % 1. Assume: X has NaN iff NAN_TEST = NaN. 2. Avoid enormous calls to IS_NAN
            y = fortran.merge('tsource', nan_test, 'fsource', min(x, [], 'all'), 'mask', infnan_obj.is_nan_sp(nan_test));

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(~any(x < y, 'all'), "No entry of X is smaller than Y", srname);
                debug_obj.assert((~infnan_obj.is_nan_sp(y)) || any(infnan_obj.is_nan(x), 'all'), "Y is not NaN unless X contains NaN", srname);
                debug_obj.assert(infnan_obj.is_nan_sp(y) || ~any(infnan_obj.is_nan(x), 'all'), "Y is NaN if X contains NaN", srname);
            end
        end
        function y = minimum2(~, x)
            %--------------------------------------------------------------------------------------------------%
            % This function returns NaN if X contains NaN; otherwise, it returns MINVAL(X). Matrix version.
            % F2018 does not specify MINVAL(X) when X contains NaN, which motivates this function. The behavior
            % of this function is the same as the following functions in various languages:
            % MATLAB: min(x, [], 'all', 'includenan')
            % Python: numpy.min(x)
            % Julia: minimum(x)
            % R: min(x)
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();

            % Inputs

            % Outputs
            y = NaN;
            % Local variables
            srname = "MINIMUM2";
            nan_test = NaN;

            %====================%
            % Calculation starts %
            %====================%

            %y = merge(tsource=sum(x), fsource=minval(x), mask=any(is_nan(x)))
            nan_test = sum(abs(x), 'all'); % 1. Assume: X has NaN iff NAN_TEST = NaN. 2. Avoid enormous calls to IS_NAN
            y = fortran.merge('tsource', nan_test, 'fsource', min(x, [], 'all'), 'mask', infnan_obj.is_nan_sp(nan_test));

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(~any(x < y, 'all'), "No entry of X is smaller than Y", srname);
                debug_obj.assert((~infnan_obj.is_nan_sp(y)) || any(infnan_obj.is_nan(x), 'all'), "Y is not NaN unless X contains NaN", srname);
                debug_obj.assert(infnan_obj.is_nan_sp(y) || ~any(infnan_obj.is_nan(x), 'all'), "Y is NaN if X contains NaN", srname);
            end
        end
        function y = maximum1(~, x)
            %--------------------------------------------------------------------------------------------------%
            % This function returns NaN if X contains NaN; otherwise, it returns MAXVAL(X). Vector version.
            % F2018 does not specify MAXVAL(X) when X contains NaN, which motivates this function. The behavior
            % of this function is the same as the following functions in various languages:
            % MATLAB: max(x, [], 'includenan')
            % Python: numpy.max(x)
            % Julia: maximum(x)
            % R: max(x)
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();

            % Inputs

            % Outputs
            y = NaN;
            % Local variables
            srname = "MAXIMUM1";
            nan_test = NaN;

            %====================%
            % Calculation starts %
            %====================%

            %y = merge(tsource=sum(x), fsource=maxval(x), mask=any(is_nan(x)))
            nan_test = sum(abs(x), 'all'); % 1. Assume: X has NaN iff NAN_TEST = NaN. 2. Avoid enormous calls to IS_NAN
            y = fortran.merge('tsource', nan_test, 'fsource', max(x, [], 'all'), 'mask', infnan_obj.is_nan_sp(nan_test));

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(~any(x > y, 'all'), "No entry of X is larger than Y", srname);
                debug_obj.assert((~infnan_obj.is_nan_sp(y)) || any(infnan_obj.is_nan(x), 'all'), "Y is not NaN unless X contains NaN", srname);
                debug_obj.assert(infnan_obj.is_nan_sp(y) || ~any(infnan_obj.is_nan(x), 'all'), "Y is NaN if X contains NaN", srname);
            end
        end
        function y = maximum2(~, x)
            %--------------------------------------------------------------------------------------------------%
            % This function returns NaN if X contains NaN; otherwise, it returns MAXVAL(X). Matrix version.
            % F2018 does not specify MAXVAL(X) when X contains NaN, which motivates this function. The behavior
            % of this function is the same as the following functions in various languages:
            % MATLAB: max(x, [], 'all', 'includenan')
            % Python: numpy.max(x)
            % Julia: maximum(x)
            % R: max(x)
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();

            % Inputs

            % Outputs
            y = NaN;
            % Local variables
            srname = "MAXIMUM2";
            nan_test = NaN;

            %====================%
            % Calculation starts %
            %====================%

            %y = merge(tsource=sum(x), fsource=maxval(x), mask=any(is_nan(x)))
            nan_test = sum(abs(x), 'all'); % 1. Assume: X has NaN iff NAN_TEST = NaN. 2. Avoid enormous calls to IS_NAN
            y = fortran.merge('tsource', nan_test, 'fsource', max(x, [], 'all'), 'mask', infnan_obj.is_nan_sp(nan_test));

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(~any(x > y, 'all'), "No entry of X is larger than Y", srname);
                debug_obj.assert((~infnan_obj.is_nan_sp(y)) || any(infnan_obj.is_nan(x), 'all'), "Y is not NaN unless X contains NaN", srname);
                debug_obj.assert(infnan_obj.is_nan_sp(y) || ~any(infnan_obj.is_nan(x), 'all'), "Y is NaN if X contains NaN", srname);
            end
        end
        function x = linspace_r(~, xstart, xstop, n)
            %--------------------------------------------------------------------------------------------------%
            % Similar to the function `linspace` in MATLAB and Python, this function generates N evenly spaced
            % numbers, the space between the consecutive points being (XSTOP-XSTART)/(N-1).
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();

            % Inputs



            % Outputs
            x = NaN(max(n, 0), 1);
            % Local variables
            srname = "LINSPACE_R";
            i = NaN;
            nm = NaN;
            xunit = NaN;

            %====================%
            % Calculation starts %
            %====================%

            if n <= 0
                % Quick return when N <= 0.
                return
            end

            nm = n - 1;

            if n == 1 || (xstart <= xstop && xstop <= xstart)
                x = repmat(xstop, size(x));
            elseif abs(xstart) <= abs(xstop) && abs(xstop) <= abs(xstart)
                xunit = xstop / double(nm);
                x(:) = xunit * double(reshape((-nm:2:nm), [], 1));
                if mod(nm, 2) == 0
                    x(1 + nm / 2) = consts_obj.ZERO;
                end
            else
                xunit = (xstop - xstart) / double(nm);
                x(:) = xstart + xunit * double(reshape((0:nm), [], 1));
            end

            if n >= 1
                % Indeed, N < 1 cannot happen due to the quick return when N <= 0.
                x(1) = xstart;
                x(n) = xstop;
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(x) == max(n, 0), "SIZE(X) == MAX(N, 0)", srname);
            end
        end
        function x = linspace_i(obj, xstart, xstop, n)
            %--------------------------------------------------------------------------------------------------%
            % This function returns INT(LINSPACE_R(REAL(XSTART, RP), REAL(XSTOP, RP), N), IK).
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();

            % Inputs



            % Outputs
            x = NaN(max(n, 0), 1);
            % Local variables
            srname = "LINSPACE_I";

            %====================%
            % Calculation starts %
            %====================%

            x(:) = round(obj.linspace_r(double(xstart), double(xstop), n)); % Rounded to the closest integer.

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(x) == max(n, 0), "SIZE(X) == MAX(N, 0)", srname);
            end
        end
        function [A, tdiag, tsubdiag] = hessenberg_hhd_trid(obj, A, tdiag, tsubdiag)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine applies Householder transformations to obtain a tridiagonal matrix that is similar
            % to a SYMMETRIC matrix A. The tridiagonal matrix is the Hessenberg form of A; its diagonal will be
            % stored in TDIAD, and the subdiagonal in TSUBDIAG. At the return, the matrix A will be DESTROYED
            % and its lower triangular part will store the Householder vectors. The code is retrieved from
            % Powell's trust region subproblem solver in UOBYQA.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();

            % In-outputs

            % Outputs


            % Local variables
            srname = "HESSENBERG_HHD_TRID";
            i = NaN;
            j = NaN;
            k = NaN;
            n = NaN;
            Asubd = NaN;
            colsq = NaN;
            scaling = NaN;
            w = NaN(size(A, 1), 1);
            wz = NaN;
            z = NaN(size(A, 1), 1);
            scaled = false;

            % Sizes
            n = size(A, 1);

            % Preconditions
            if consts_obj.DEBUGGING
                % Even though we only need the lower triangular part of A, we assume that, in our project,
                % something is wrong if this subroutine is invoked with a non-symmetric matrix A.
                debug_obj.assert(obj.issymmetric(A), "A is symmetric", srname);
                debug_obj.assert(numel(tdiag) == n, "SIZE(TDIAG) == N", srname);
                debug_obj.assert(numel(tsubdiag) == max(0, n - 1), "SIZE(TDIAG) == MAX(0, N-1)", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            if n <= 0
                % Quick return when N <= 0. Of course, N < 0 is impossible.
                return
            end

            % According to a test on 20220508, scaling enhances the stability and slightly improves the
            % performance of UOBYQA. Indeed, when A contains huge values, NaN can occur if no scaling is applied.
            scaling = max(abs(A), [], 'all');
            scaled = false;
            if scaling <= 0
                tdiag = repmat(consts_obj.ZERO, size(tdiag));
                tsubdiag = repmat(consts_obj.ZERO, size(tsubdiag));
                return
            elseif scaling > 1.0e8 || scaling < 1.0e-4
                % The thresholds are empirical.
                A(:, :) = A ./ scaling;
                scaled = true;
            end

            tdiag(:) = obj.diag(A);

            for k = 1:n - 1
                colsq = sum(A(k + 2:n, k) .^ 2, 'all');
                if colsq <= 0
                    tsubdiag(k) = A(k + 1, k); % A(K+1, K) may have been updated in previous loops.
                    A(k + 1, k) = consts_obj.ZERO;
                    continue
                end

                Asubd = A(k + 1, k);
                tsubdiag(k) = abs(sqrt(colsq + Asubd ^ 2)) .* ((Asubd >= 0) * 2 - 1);

                A(k + 1, k) = -colsq / (Asubd + tsubdiag(k));
                w(k + 1:n) = sqrt(consts_obj.TWO / (colsq + A(k + 1, k) ^ 2)) * A(k + 1:n, k);
                %----------------------------------------------------------------------------------------------%
                % The two lines above are from Powell. They are equivalent to the following two lines.
                % %A(K + 1, K) = A(K + 1, K) - ASUBD
                % %W(K + 1:N) = sqrt(TWO) * A(K + 1:N, K) / NORM(A(K + 1:N, K))
                %----------------------------------------------------------------------------------------------%
                A(k + 1:n, k) = w(k + 1:n);

                z(k + 1:n) = tdiag(k + 1:n) .* w(k + 1:n);
                for j = k + 1:n - 1
                    z(j + 1:n) = z(j + 1:n) + A(j + 1:n, j) * w(j);
                    for i = j + 1:n
                        z(j) = z(j) + A(i, j) * w(i);
                    end
                end
                wz = obj.inprod(w(k + 1:n), z(k + 1:n));

                tdiag(k + 1:n) = tdiag(k + 1:n) + w(k + 1:n) .* (wz * w(k + 1:n) - consts_obj.TWO * z(k + 1:n));
                for j = k + 1:n
                    A(j + 1:n, j) = A(j + 1:n, j) - w(j + 1:n) * z(j) - w(j) * (z(j + 1:n) - wz * w(j + 1:n));
                end
            end

            if scaled
                tdiag(:) = tdiag * scaling;
                tsubdiag(:) = tsubdiag * scaling;
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(tdiag) == n, "SIZE(TDIAG) == N", srname);
                debug_obj.assert(numel(tsubdiag) == max(0, n - 1), "SIZE(TDIAG) == MAX(0, N-1)", srname);
            end
        end
        function [H, Q] = hessenberg_full(obj, A, H, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine finds a Hessenberg matrix H (all entries below the subdiagonal are 0) such that
            % H = Q^T*A*Q, where Q is a orthogonal matrix that may also be returned. A will stay unchanged.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();

            % Inputs

            % Outputs


            % Local variables
            srname = "HESSENBERG_FULL";
            i = NaN;
            j = NaN;
            n = NaN;
            colsq = NaN;
            subd = NaN;
            v = NaN(size(A, 1), 1);
            w = NaN(size(A, 1), 1);
            scaling = NaN;
            scaled = false;

            % Debugging variables
            tol = NaN;

            % Sizes
            n = size(A, 1);

            %====================%
            % Calculation starts %
            %====================%

            % Preconditions
            ipObj = inputParser();
            addParameter(ipObj, 'Q', NaN);
            parse(ipObj, varargin{:});
            Q = ipObj.Results.Q;
            if consts_obj.DEBUGGING
                debug_obj.assert(size(A, 1) == size(A, 2), "A is square", srname);
                debug_obj.assert(size(H, 1) == n && size(H, 2) == n, "SIZE(H) == [N, N]", srname);
                if nargout >= 2
                    debug_obj.assert(size(Q, 1) == n && size(Q, 2) == n, "SIZE(Q) == [N, N]", srname);
                end
            end

            if n <= 0
                % Quick return when N <= 0. Of course, N < 0 is impossible.
                return
            end

            H(:, :) = A;
            if nargout >= 2
                Q(:, :) = obj.eye1(n);
            end

            % According to a test on 20220508, scaling enhances the stability and slightly improves the
            % performance of UOBYQA. Indeed, when A contains huge values, NaN can occur if no scaling is applied.
            scaling = max(abs(H), [], 'all');
            scaled = false;
            if scaling <= 0
                return
            elseif scaling > 1000000.0 || scaling < 1.0e-6
                % 1.0E6 and 1.0E-6 are heuristic.
                H(:, :) = H ./ scaling;
                scaled = true;
            end

            for j = 1:n - 1
                colsq = sum(H(j + 2:n, j) .^ 2, 'all');
                if colsq <= 0
                    continue
                end

                v(j + 1:n) = H(j + 1:n, j);
                subd = abs(sqrt(v(j + 1) ^ 2 + colsq)) .* ((v(j + 1) >= 0) * 2 - 1);

                %----------------------------------------------------------------------------------------------%
                v(j + 1) = -colsq / (v(j + 1) + subd);
                v(j + 1:n) = sqrt(consts_obj.TWO / (colsq + v(j + 1) ^ 2)) * v(j + 1:n);
                % The two lines above are from Powell. They are equivalent to the following two lines.
                % %V(J + 1) = V(J + 1) - SUBD
                % %V(J + 1:N) = sqrt(TWO) * V(J + 1:N) / NORM(V(J + 1:N))
                %----------------------------------------------------------------------------------------------%

                for i = j + 1:n
                    H(j + 1:n, i) = H(j + 1:n, i) - obj.inprod(H(j + 1:n, i), v(j + 1:n)) * v(j + 1:n);
                end
                H(j + 1, j) = subd;
                H(j + 2:n, j) = consts_obj.ZERO;

                w(:) = obj.matprod21(H(:, j + 1:n), v(j + 1:n));
                for i = j + 1:n
                    H(:, i) = H(:, i) - w * v(i);
                end

                if nargout >= 2
                    w(:) = obj.matprod21(Q(:, j + 1:n), v(j + 1:n));
                    for i = j + 1:n
                        Q(:, i) = Q(:, i) - w * v(i);
                    end
                end
            end

            if scaled
                H(:, :) = H * scaling;
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(size(H, 1) == n && size(H, 2) == n, "SIZE(H) == [N, N]", srname);
                debug_obj.assert(obj.isbanded(H, 1, n - 1), "H is a Hessenberg matrix", srname);
                tol = max(consts_obj.TEN ^ max(-8, -consts_obj.MAXPOW10), min(0.1, consts_obj.TEN ^ min(10, consts_obj.MAXPOW10) * consts_obj.EPS * double(n)));
                debug_obj.assert(obj.issymmetric(H, 'tol', tol) || ~obj.issymmetric(A), "H is symmetric if so is A", srname);
                if nargout >= 2
                    debug_obj.assert(size(Q, 1) == n && size(Q, 2) == n, "SIZE(Q) == [N, N]", srname);
                    debug_obj.assert(obj.isorth(Q, 'tol', tol), "Q is orthogonal", srname);
                    debug_obj.assert(all(abs(obj.matprod22(Q, H) - obj.matprod22(A, Q)) <= tol * max(abs(A), [], 'all'), 'all'), "Q*H = A*Q", srname);
                end
            end
        end
        function eig_min = eigmin_sym_trid(obj, td, tn, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This function approximates the smallest eigenvalue EIG_MIN of a symmetric tridiagonal matrix by a
            % bisection method, in which process EMINLB is a lower bound on EIG_MIN and EMINUB an upper bound.
            % TD is the diagonal of the tridiagonal matrix, and TN is the subdiagonal and superdiagonal. EMINUB
            % is occasionally adjusted by the rule of false position (https:
            % which attempts to accelerate the bisection by linear interpolation. The code is retrieved from
            % Powell's trust region subproblem solver in UOBYQA.
            %
            % The bisection algorithm for eigenvalues (not only the smallest) of symmetric tridiagonal matrices
            % can be found in
            % Barth, Martin, and Wilkinson, Calculation of the eigenvalues of a symmetric tridiagonal matrix by
            % the method of bisection, Numerische Mathematik 9, 386--393 (1967).
            % The algorithm is based on the sign changes of the Sturm sequence {P_i(LAMBDA)} defined in (1)--(2)
            % of the above mentioned paper (P_i(LAMBDA) is the determinant of the i-th principle submatrix of
            % the matrix minus LAMBDA*I), or the number of negative values of the Sturm-ratio sequence
            % {Q_i(LAMBDA) = P_i(LAMBDA)/P_{i-1}(LAMBDA)} in (3)--(5) of the paper. The theoretical basis is the
            % following fact: for any symmetric tridiagonal n-by-n matrix A, the number of negative eigenvalues
            % is equal to the number of sign changes in the Sturm sequence 1, det(A^(1)), det(A^(2)), ...,
            % det(A^(n)), where A^(i) is the i-th principle submatrix of A, where a "sign change" is a transition
            % from nonpositive to positive or from nonnegative to negative (see, e.g., pages 228--229 of
            % Trefethen-Bau 1997, Numerical Linear Algebra, or pages 300--301 of Wilkinson 1965, The Algebraic
            % Eigenvalue Problem).
            %
            % In MATLAB/Python/Julia/R, to get the smallest eigenvalue, we should use the eigenvalue computation
            % function built in the languages or standard libraries. For example, in MATLAB, we can do
            % %tridh = spdiags([[tn; 0], td, [0; tn]], -1:1, n, n);
            % %crvmin = eigs(tridh, 1, 'smallestreal');
            % %% It is critical for the efficiency to use `spdiags` to construct `tridh` in the sparse form.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();

            % Inputs



            % Outputs
            eig_min = NaN;
            % Local variables
            srname = "EIGMIN";
            iter = NaN;
            k = NaN;
            ksav = NaN;
            maxiter = NaN;
            n = NaN;
            eminlb = NaN;
            eminub = NaN;
            piv = NaN(numel(td), 1);
            pivksv = NaN;
            pivnew = NaN(numel(td), 1);
            tol_loc = NaN;

            % Sizes
            n = fix(numel(td));

            % Preconditions
            ipObj = inputParser();
            addParameter(ipObj, 'tol', NaN);
            parse(ipObj, varargin{:});
            tol = ipObj.Results.tol;
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(tn) == n - 1, "SIZE(TN) == N - 1", srname);
                if ~ismember('tol', ipObj.UsingDefaults)
                    debug_obj.assert(tol >= 0, "TOL >= 0", srname);
                end
            end

            %====================%
            % Calculation starts %
            %====================%

            maxiter = 100;
            tol_loc = consts_obj.TEN ^ max(-6, -consts_obj.MAXPOW10);
            if ~ismember('tol', ipObj.UsingDefaults)
                tol_loc = tol;
            end

            % The following loop calculates the Sturm ratios [Q_1(0), ..., Q_n(0)]. These ratios are all positive
            % iff all the eigenvalues of the matrix are positive definite. Note that these ratios are also the
            % pivots of the Cholesky factorization of the matrix (i.e., the square of the diagonal of L in LL^T,
            % or the diagonal of D in LDL^T). All the pivots are positive iff there exists a Cholesky
            % factorization with a positive diagonal, i.e., the matrix is positive definite.
            piv = repmat(-consts_obj.ONE, size(piv));
            piv(1) = td(1);
            for k = 1:n - 1
                if piv(k) > 0
                    piv(k + 1) = td(k + 1) - tn(k) ^ 2 / piv(k);
                else
                    break
                end
            end

            if all(piv >= 0, 'all')
                % The matrix is positive semidefinite.
                eminub = min(piv, [], 'all');
                eminlb = consts_obj.ZERO;
            else
                eminub = min(td, [], 'all');
                eminlb = -max(abs([consts_obj.ZERO; reshape(tn, [], 1)]) + abs(td) + abs([reshape(tn, [], 1); consts_obj.ZERO]), [], 'all');
            end

            ksav = 0;
            pivksv = consts_obj.ZERO; % This initial value will not be used, but Fortran compilers may complain without it.
            for iter = 1:maxiter                % Powell's code is essentially a DO WHILE loop. We impose an explicit MAXITER.
                if eminub - eminlb <= tol_loc * max(abs(eminlb), abs(eminub))
                    break
                end
                eig_min = consts_obj.HALF * (eminlb + eminub);

                % The following loop calculates the Sturm ratios [Q_1(EIG_MIN), ..., Q_n(EIG_MIN)]. These ratios
                % are all positive iff all the eigenvalues of the matrix are larger than EIG_MIN, i.e., EIG_MIN
                % underestimates the smallest eigenvalue. Note that these ratios are also the pivots of the
                % Cholesky factorization of the matrix minus EIG_MIN*I (i.e., the square of the diagonal of L in
                % LL^T, or the diagonal of D in LDL^T). All the pivots are positive iff there exists a Cholesky
                % factorization with a positive diagonal, i.e., the matrix minus LAMBDA*I is positive definite.
                pivnew = repmat(-consts_obj.ONE, size(pivnew));
                pivnew(1) = td(1) - eig_min;
                for k = 1:n - 1
                    if pivnew(k) > 0
                        pivnew(k + 1) = td(k + 1) - eig_min - tn(k) ^ 2 / pivnew(k);
                    else
                        break
                    end
                end

                if all(pivnew > 0, 'all')
                    piv(:) = pivnew;
                    eminlb = eig_min;
                    continue
                end

                % We arrive here iff PIVNEW contains nonpositive entries and EIG_MIN is no less than the smallest
                % eigenvalue. We set EMINUB to EIG_MIN except a possible adjustment by the rule of false position.
                k = min(obj.trueloc(~(pivnew > 0)), [], 'all');
                piv(1:k - 1) = pivnew(1:k - 1);

                % KSAV was initialized to 0, triggering the ELSE when ALL(PIVNEW > 0) fails for the first time.
                if k == ksav && pivksv < 0 && piv(k) - pivnew(k) >= pivnew(k) - pivksv
                    pivksv = consts_obj.ZERO;
                    eminub = (eig_min * piv(k) - eminlb * pivnew(k)) / (piv(k) - pivnew(k));
                else
                    ksav = k;
                    pivksv = pivnew(k); % PIVKSAV <= 0.
                    eminub = eig_min;
                end

                %----------------------------------------------------------------------------------------------%
                % Powell's original code contains the following, why? It seems to cause wrong outputs sometimes.
                % Zaikun 20220511: Does this affect the adjustment by the rule of false position?
                % %IF (K < KSAV .OR. (K == KSAV .AND. PIVKSV == 0)) EXIT
                %----------------------------------------------------------------------------------------------%
            end

            eig_min = eminlb;

            %====================%
            %  Calculation ends  %
            %====================%
        end
        function smat = vec2smat(obj, vec)
            %--------------------------------------------------------------------------------------------------%
            % This function transforms a vector VEC to a symmetric matrix SMAT with the vector storing the upper
            % triangular part of the matrix column by column.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            % Inputs

            % Outputs
            smat = NaN((round(sqrt(double(8 * numel(vec) + 1))) - 1) / 2);
            % Local variables
            srname = "SMAT2VEC";
            ih = NaN;
            j = NaN;
            n = NaN;

            % Sizes
            n = size(smat, 1);

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(vec) == n * (n + 1) / 2, "SIZE(VEC) = N*(N+1)/2", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            for j = 1:n
                ih = (j - 1) * j / 2;
                smat(1:j, j) = vec(ih + 1:ih + j);
                smat(j, 1:j - 1) = smat(1:j - 1, j);
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(obj.issymmetric(smat), "SMAT is symmetric", srname);
            end
        end
        function vec = smat2vec(obj, smat)
            %--------------------------------------------------------------------------------------------------%
            % This function transforms a symmetric matrix SMAT to a vector VEC that stores the upper triangular
            % part of the matrix column by column.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            % Inputs

            % Outputs
            vec = NaN((size(smat, 1) * (size(smat, 1) + 1)) / 2, 1);
            % Local variables
            srname = "SMAT2VEC";
            ih = NaN;
            n = NaN;
            j = NaN;

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(obj.issymmetric(smat), "SMAT is symmetric", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            n = size(smat, 1);
            for j = 1:n
                ih = (j - 1) * j / 2;
                vec(ih + 1:ih + j) = smat(1:j, j);
            end

            %====================%
            % Calculation ends   %
            %====================%

        end
        function y = smat_mul_vec(obj, smatv, x)
            %--------------------------------------------------------------------------------------------------%
            % This function calculates the product of a symmetric matrix and a vector X, with the upper
            % triangular part of the matrix stored in the vector SMATV column by column.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            % Inputs


            % Outputs
            y = NaN(numel(x), 1);
            % Local variables
            srname = "SMAT_MUL_VEC";
            ih = NaN;
            n = NaN;
            j = NaN;

            % Sizes
            n = fix(numel(x));

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(smatv) == n * (n + 1) / 2, "SIZE(SMATV) = N*(N+1)/2", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            for j = 1:n
                ih = (j - 1) * j / 2;
                y(j) = obj.inprod(smatv(ih + 1:ih + j), x(1:j));
                y(1:j - 1) = y(1:j - 1) + x(j) * smatv(ih + 1:ih + j - 1);
            end

            %====================%
            % Calculation ends   %
            %====================%

        end

    end
end
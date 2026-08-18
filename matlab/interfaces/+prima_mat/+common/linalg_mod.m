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
        function A = r1_sym(obj, A, alpha, x)
            %--------------------------------------------------------------------------------------------------%
            % R1_SYM sets
            % A = A + ALPHA*( X*X^T ),
            % where A is an NxN matrix, ALPHA is a scalar, and X is an N-dimensional vector.
            %--------------------------------------------------------------------------------------------------%


            % A(SIZE(X), SIZE(X))


            n = numel(x);

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


        end
        function A = r1(obj, A, alpha, x, y)
            %--------------------------------------------------------------------------------------------------%
            % R1 sets
            % A = A + ALPHA*( X*Y^T ),
            % where A is an MxN matrix, ALPHA is a real scalar, X is an M-dimensional vector, and Y is an
            % N-dimensional vector.
            %--------------------------------------------------------------------------------------------------%


            % A(SIZE(X), SIZE(Y))


            %====================%
            % Calculation starts %
            %====================%

            % N.B.: The use of OUTPROD is expensive memory-wise, but it is not our concern in this implementation.
            A(:, :) = A + alpha * x * y.';
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


            % A(SIZE(X), SIZE(X))


            n = numel(x);

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


        end
        function A = r2(obj, A, alpha, x, y, beta, u, v)
            %--------------------------------------------------------------------------------------------------%
            % R2 sets
            % A = A + ( ALPHA*( X*Y^T ) + BETA*( U*V^T ) ),
            % where A is an MxN matrix, ALPHA and BETA are real scalars, X and U are M-dimensional vectors,
            % Y and V are N-dimensional vectors.
            %--------------------------------------------------------------------------------------------------%


            % U(SIZE(X))
            % V(SIZE(Y))

            % A(SIZE(X), SIZE(Y))


            %====================%
            % Calculation starts %
            %====================%

            % N.B.: The use of OUTPROD is expensive memory-wise, but it is not our concern in this implementation.
            A(:, :) = A + alpha * x * y.' + beta * u * v.';
            %A = A + (alpha * outprod(x, y) + beta * outprod(u, v))

            %====================%
            %  Calculation ends  %
            %====================%
        end
        function is_inv = isinv(obj, A, B, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This procedure tests whether A = B^{-1} up to the tolerance TOL.
            %--------------------------------------------------------------------------------------------------%


            is_inv = false;

            n = size(A, 1);

            ipObj = inputParser();
            addParameter(ipObj, 'tol', NaN);
            parse(ipObj, varargin{:});
            tol = ipObj.Results.tol;

            %====================%
            % Calculation starts %
            %====================%

            if ismember('tol', ipObj.UsingDefaults)
                tol_loc = min(1.0e-3, 100.0 * eps(1.0) * double(max(size(A, 1), size(A, 2))));
            else
                tol_loc = tol;
            end
            tol_loc = max([tol_loc, tol_loc * max(abs(A), [], 'all'), tol_loc * max(abs(B), [], 'all')], [], 'all');
            is_inv = all(abs(A * B - eye(n)) <= tol_loc, 'all') || all(abs(B * A - eye(n)) <= tol_loc, 'all');

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


            Q_loc = NaN(size(A, 1));
            T = NaN(size(A, 2), size(A, 1));

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

            m = size(A, 1);
            n = size(A, 2);

            %====================%
            % Calculation starts %
            %====================%

            pivot = (nargout >= 3);
            Q_loc(:, :) = eye(m);
            T(:, :) = A.'; % T is the transpose of R. We consider T in order to work on columns.
            if pivot
                P = (1:n).';
            end

            for j = 1:n
                if pivot
                    k = fortran.maxloc(sum(T(j:n, j:m) .^ 2, 2), 'dim', 1);
                    if k > 1 && k <= n - j + 1
                        k = k + j - 1;
                        P([j, k]) = P([k, j]);
                        T([j, k], :) = T([k, j], :);
                    end
                end
                for i = m:-1:j + 1
                    G = obj.planerot(T(j, [j, i]).').';
                    T(j, [j, i]) = [obj.hypotenuse(T(j, j), T(j, i)), 0.0]; %T(j, [j, i]) = [sqrt(T(j, j)**2 + T(j, i)**2), ZERO]
                    T(j + 1:n, [j, i]) = T(j + 1:n, [j, i]) * G;
                    Q_loc(:, [j, i]) = Q_loc(:, [j, i]) * G;
                end
            end

            if nargout >= 1
                Q = Q_loc(:, 1:size(Q, 2));
            end
            if nargout >= 2
                R = T(:, 1:size(R, 1)).';
            end

            %====================%
            %  Calculation ends  %
            %====================%


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


            % A(M, N)
            % B(M)
            % Q(M, :), SIZE(Q, 2) = M or MIN(M, N)
            % Rdiag(MIN(M, N))

            x = NaN(size(A, 2), 1);

            P = NaN(size(A, 2), 1);

            Q_loc = NaN(size(A, 1), min(size(A, 1), size(A, 2)));
            Rdiag_loc = NaN(min(size(A, 1), size(A, 2)), 1);

            y = NaN(numel(b), 1);

            m = size(A, 1);
            n = size(A, 2);

            ipObj = inputParser();
            addParameter(ipObj, 'Q', NaN);
            addParameter(ipObj, 'Rdiag', NaN);
            parse(ipObj, varargin{:});
            Q = ipObj.Results.Q;
            Rdiag = ipObj.Results.Rdiag;

            %====================%
            % Calculation starts %
            %====================%

            if n <= 0
                % Of course, N < 0 should never happen.
                return
            end

            if ismember('Q', ipObj.UsingDefaults)
                [Q_loc, ~, P] = obj.qr(A);
                Rdiag_loc(:) = arrayfun(@(i) sum(Q_loc(:, i) .* A(:, P(i)), 'all'), (1:min(m, n))');
                %%MATLAB: Rdiag_loc = sum(Q_loc(:, 1:min(m,n)) .* A(:, P(1:min(m,n))), 1); % Row vector
                rank = max([0; find(abs(Rdiag_loc) > 0)], [], 'all');
                pivot = true;
            else
                Q_loc(:, :) = Q(:, 1:size(Q_loc, 2));
                if ismember('Rdiag', ipObj.UsingDefaults)
                    Rdiag_loc(:) = arrayfun(@(i) sum(Q_loc(:, i) .* A(:, i), 'all'), (1:min(m, n))');
                    %%MATLAB: Rdiag_loc = sum(Q_loc(:, 1:min(m,n)) .* A(:, 1:min(m,n)), 1); % Row vector
                else
                    Rdiag_loc(:) = Rdiag;
                end
                rank = min(m, n);
                pivot = false;
            end

            x(:) = 0.0;
            y(:) = b; % Local copy of B; B is INTENT(IN) and should not be modified.

            for i = rank:-1:1
                if pivot
                    j = P(i);
                else
                    j = i;
                end
                % The following IF comes from Powell. It forces X(J) = 0 if deviations from this value can be
                % attributed to computer rounding errors. This is a favorable choice in the context of COBYLA.
                yq = sum(y .* Q_loc(:, i), 'all');
                yqa = sum(abs(y) .* abs(Q_loc(:, i)), 'all');
                if obj.isminor0(yq, yqa)
                    x(j) = 0.0;
                else
                    x(j) = yq / Rdiag_loc(i);
                    y = y - x(j) * A(:, j);
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


            % B(M)
            % Q(M, N)
            % R(N, N)

            x = NaN(size(R, 2), 1);

            n = size(R, 2);

            %====================%
            % Calculation starts %
            %====================%

            if n <= 0
                % Of course, N < 0 should never happen.
                return
            end

            x(:) = Q.' * b;
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
        function is_banded = isbanded(~, A, lwidth, uwidth, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This function tests whether the matrix A banded within the bandwidth specified by LWIDTH and
            % UWIDTH up to the tolerance TOL.
            %--------------------------------------------------------------------------------------------------%


            is_banded = false;

            ipObj = inputParser();
            addParameter(ipObj, 'tol', NaN);
            parse(ipObj, varargin{:});
            tol = ipObj.Results.tol;

            %====================%
            % Calculation starts %
            %====================%

            tol_loc = 0.0;
            if ~ismember('tol', ipObj.UsingDefaults)
                tol_loc = max(tol, tol * max(abs(A), [], 'all'));
            end
            if isnan(tol_loc)
                tol_loc = 0.0;
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


            is_tril = false;

            ipObj = inputParser();
            addParameter(ipObj, 'tol', NaN);
            parse(ipObj, varargin{:});
            tol = ipObj.Results.tol;

            %====================%
            % Calculation starts %
            %====================%

            if ismember('tol', ipObj.UsingDefaults)
                tol_loc = 0.0;
            else
                tol_loc = tol;
            end
            width = max(0, size(A, 1) - 1);
            is_tril = obj.isbanded(A, width, 0, 'tol', tol_loc);

            %====================%
            %  Calculation ends  %
            %====================%
        end
        function is_triu = istriu(obj, A, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This function tests whether the matrix A is upper triangular up to the tolerance TOL.
            %--------------------------------------------------------------------------------------------------%


            is_triu = false;

            ipObj = inputParser();
            addParameter(ipObj, 'tol', NaN);
            parse(ipObj, varargin{:});
            tol = ipObj.Results.tol;

            %====================%
            % Calculation starts %
            %====================%

            if ismember('tol', ipObj.UsingDefaults)
                tol_loc = 0.0;
            else
                tol_loc = tol;
            end
            width = max(0, size(A, 2) - 1);
            is_triu = obj.isbanded(A, 0, width, 'tol', tol_loc);

            %====================%
            %  Calculation ends  %
            %====================%
        end
        function is_orth = isorth(obj, A, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This function tests whether the matrix A has orthonormal columns up to the tolerance TOL.
            %--------------------------------------------------------------------------------------------------%


            is_orth = false;

            ipObj = inputParser();
            addParameter(ipObj, 'tol', NaN);
            parse(ipObj, varargin{:});
            tol = ipObj.Results.tol;

            %====================%
            % Calculation starts %
            %====================%

            tol_loc = realmax;
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
            elseif any(isnan(A), 'all')
                is_orth = false;
            elseif realmax < realmax
                is_orth = all(abs(A.' * A - eye(n)) <= max(tol_loc, tol_loc * max(abs(A), [], 'all')), 'all');
            end

            %====================%
            %  Calculation ends  %
            %====================%
        end
        function y = project1(obj, x, v)
            %--------------------------------------------------------------------------------------------------%
            % This function returns the projection of X to SPAN(V).
            %--------------------------------------------------------------------------------------------------%


            y = NaN(numel(x), 1);

            u = NaN(numel(v), 1);

            %====================%
            % Calculation starts %
            %====================%

            if all(abs(x) <= 0, 'all') || all(abs(v) <= 0, 'all')
                y(:) = 0.0;
            elseif any(isnan(x), 'all') || any(isnan(v), 'all')
                y(:) = sum(x, 'all') + sum(v, 'all'); % Set Y to NaN

            elseif any(isinf(v), 'all')
                u(:) = 0.0;
                u(isinf(v)) = 1.0 .* ((v(isinf(v)) > 0) .* 2 - 1);
                %%MATLAB: u = 0; u(isinf(v)) = sign(v(isinf(v)))
                u = u ./ norm(u);
                y(:) = sum(x .* u, 'all') * u;
            else
                u(:) = v ./ norm(v);
                y(:) = sum(x .* u, 'all') * u;
            end

            %====================%
            %  Calculation ends  %
            %====================%


        end
        function y = project2(obj, x, V)
            %--------------------------------------------------------------------------------------------------%
            % This function returns the projection of X to RANGE(V).
            %--------------------------------------------------------------------------------------------------%


            y = NaN(numel(x), 1);

            V_loc = NaN(size(V, 1), size(V, 2));

            %====================%
            % Calculation starts %
            %====================%

            if size(V, 2) == 1
                y = obj.project1(x, V(:, 1));
            elseif all(abs(x) <= 0, 'all') || all(abs(V) <= 0, 'all')
                y(:) = 0.0;
            elseif any(isnan(x), 'all') || any(isnan(V), 'all')
                y(:) = sum(x, 'all') + sum(V, 'all'); % Set Y to NaN

            elseif any(isinf(V), 'all')
                mask00 = isinf(V);
                V_loc(mask00) = 1.0 .* ((V(mask00) > 0) .* 2 - 1);
                mask01 = ~mask00;
                V_loc(mask01) = 0.0;

                %%MATLAB: V_loc = 0; V_loc(isinf(V)) = sign(V);
                U = obj.qr(V_loc);
                y(:) = U * (U.' * x);
            else
                U = obj.qr(V);
                y(:) = U * (U.' * x);
            end

            %====================%
            %  Calculation ends  %
            %====================%


        end
        function r = hypotenuse(~, x1, x2)
            %--------------------------------------------------------------------------------------------------%
            % HYPOTENUSE(X1, X2) returns SQRT(X1^2 + X2^2), handling over/underflow.
            %--------------------------------------------------------------------------------------------------%


            r = NaN;

            y = NaN(2, 1);

            %====================%
            % Calculation starts %
            %====================%

            if ~isfinite(x1)
                r = abs(x1);
            elseif ~isfinite(x2)
                r = abs(x2);
            else
                y(:) = abs([x1, x2]);
                y(:) = [min(y, [], 'all'), max(y, [], 'all')];
                if y(1) > sqrt(realmin) && y(2) < sqrt(realmax / 2.1)
                    r = sqrt(sum(y .^ 2, 'all'));
                elseif y(2) > 0
                    r = y(2) * sqrt((y(1) / y(2)) ^ 2 + 1.0);
                else
                    r = 0.0;
                end
                % Without the following line, R > Y(1) + Y(2) or R < Y(2) may happen due to rounding errors.
                r = min(sum(y, 'all'), max(y(2), r));
            end

            %====================%
            %  Calculation ends  %
            %====================%


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


            G = NaN(2);

            %====================%
            % Calculation starts %
            %====================%

            % Define C = X(1) / R and S = X(2) / R with R = HYPOT(X(1), X(2)). Handle Inf/NaN, over/underflow.
            if any(isnan(x), 'all')
                % In this case, MATLAB sets G to NaN(2, 2). We refrain from doing so to keep G orthogonal.
                c = 1.0;
                s = 0.0;
            elseif all(isinf(x), 'all')
                % In this case, MATLAB sets G to NaN(2, 2). We refrain from doing so to keep G orthogonal.
                c = 1 / sqrt(2.0) .* ((x(1) > 0) .* 2 - 1);
                s = 1 / sqrt(2.0) .* ((x(2) > 0) .* 2 - 1);
            elseif abs(x(1)) <= 0 && abs(x(2)) <= 0
                % X(1) == 0 == X(2).
                c = 1.0;
                s = 0.0;
            elseif abs(x(2)) <= eps(1.0) * abs(x(1))
                % N.B.:
                % 0. With <= instead of <, this case covers X(1) == 0 == X(2), which is treated above separately
                % to avoid the confusing SIGN(., 0) (see 1).
                % 1. SIGN(A, 0) = ABS(A) in Fortran but sign(0) = 0 in MATLAB, Python, Julia, and R!
                % 2. Taking SIGN(X(1)) into account ensures the continuity of G with respect to X except at 0.
                c = 1.0 .* ((x(1) > 0) .* 2 - 1); %%MATLAB: c = sign(x(1))
                s = 0.0;
            elseif abs(x(1)) <= eps(1.0) * abs(x(2))
                % N.B.: SIGN(A, X) = ABS(A) * sign of X /= A * sign of X ! Therefore, it is WRONG to define G
                % as SIGN(RESHAPE([ZERO, -ONE, ONE, ZERO], [2, 2]), X(2)). This mistake was committed on
                % 20211206 and took a whole day to debug! NEVER use SIGN on arrays unless you are really sure.
                c = 0.0;
                s = 1.0 .* ((x(2) > 0) .* 2 - 1); %%MATLAB: s = sign(x(2))

            else
                % Here is the normal case. It implements the Givens rotation in a stable & continuous way as in:
                % Bindel, D., Demmel, J., Kahan, W., and Marques, O. (2002). On computing Givens rotations
                % reliably and efficiently. ACM Transactions on Mathematical Software (TOMS), 28(2), 206-238.
                % N.B.: 1. Modern compilers compute SQRT(REALMIN) and SQRT(REALMAX/2.1) at compilation time.
                % 2. The direct calculation without involving T and U seems to work better; use it if possible.
                if all(abs(x) > sqrt(realmin) & abs(x) < sqrt(realmax / 2.1), 'all')
                    % Do NOT use HYPOTENUSE here; the best implementation for one may be suboptimal for the other
                    r = norm(x);
                    c = x(1) / r;
                    s = x(2) / r;
                elseif abs(x(1)) > abs(x(2))
                    t = x(2) / x(1);
                    u = max([1.0, abs(t), sqrt(1.0 + t ^ 2)], [], 'all'); % MAXVAL: precaution against rounding error.
                    u = u .* ((x(1) > 0) .* 2 - 1); %%MATLAB: u = sign(x(1))*sqrt(ONE + t**2)
                    c = 1.0 / u;
                    s = t / u;
                else
                    t = x(1) / x(2);
                    u = max([1.0, abs(t), sqrt(1.0 + t ^ 2)], [], 'all'); % MAXVAL: precaution against rounding error.
                    u = u .* ((x(2) > 0) .* 2 - 1); %%MATLAB: u = sign(x(2))*sqrt(ONE + t**2)
                    c = t / u;
                    s = 1.0 / u;
                end
            end

            G = reshape([c, -s, s, c], [2, 2]); %%MATLAB: G = [c, s; -s, c]

            %====================%
            %  Calculation ends  %
            %====================%


        end
        function A = symmetrize(~, A)
            %--------------------------------------------------------------------------------------------------%
            % SYMMETRIZE(A) symmetrizes A.
            % N.B.: Here, we assume that A is a matrix that IS SUPPOSED TO BE symmetric in precise arithmetic,
            % and its asymmetry comes only from errors (e.g., rounding, noise).
            %--------------------------------------------------------------------------------------------------%


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


        end
        function is_minor = isminor0(~, x, ref)
            %--------------------------------------------------------------------------------------------------%
            % This function tests whether X is minor compared to REF. It is used by Powell, e.g., in COBYLA.
            % In precise arithmetic, ISMINOR(X, REF) is TRUE if and only if X == 0; in floating-point
            % arithmetic, ISMINOR(X, REF) is true if X is zero or its nonzero value can be attributed to
            % computer rounding errors according to REF.
            % Larger SENSITIVITY means the function is more strict/precise, the value TENTH being due to Powell.
            %--------------------------------------------------------------------------------------------------%


            is_minor = false;

            sensitivity = 0.1;

            %====================%
            % Calculation starts %
            %====================%

            refa = abs(ref) + sensitivity * abs(x);
            refb = abs(ref) + 2.0 * sensitivity * abs(x);
            is_minor = (abs(ref) >= refa || refa >= refb);

            %====================%
            %  Calculation ends  %
            %====================%

        end
        function is_minor = isminor1(obj, x, ref)
            %--------------------------------------------------------------------------------------------------%
            % This function tests whether X is minor compared to REF. It is used by Powell, e.g., in COBYLA.
            %--------------------------------------------------------------------------------------------------%


            is_minor = false(numel(x), 1);

            %====================%
            % Calculation starts %
            %====================%

            is_minor(:) = arrayfun(@(i) obj.isminor0(x(i), ref(i)), (1:numel(x))');

            %====================%
            %  Calculation ends  %
            %====================%


        end
        function is_symmetric = issymmetric(~, A, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This function tests whether A is symmetric up to TOL.
            %--------------------------------------------------------------------------------------------------%


            is_symmetric = false;

            ipObj = inputParser();
            addParameter(ipObj, 'tol', NaN);
            parse(ipObj, varargin{:});
            tol = ipObj.Results.tol;

            %====================%
            % Calculation starts %
            %====================%

            tol_loc = 1.0e-10;
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
            elseif 1.0e-10 < 0.9 * realmax
                is_symmetric = (~any(abs(A - A.') > tol_loc * max(max(abs(A), [], 'all'), 1.0), 'all')) && all(isnan(A) == isnan(A.'), 'all');
            end

            %====================%
            %  Calculation ends  %
            %====================%
        end
        function [A, tdiag, tsubdiag] = hessenberg_hhd_trid(obj, A, tdiag, tsubdiag)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine applies Householder transformations to obtain a tridiagonal matrix that is similar
            % to a SYMMETRIC matrix A. The tridiagonal matrix is the Hessenberg form of A; its diagonal will be
            % stored in TDIAD, and the subdiagonal in TSUBDIAG. At the return, the matrix A will be DESTROYED
            % and its lower triangular part will store the Householder vectors. The code is retrieved from
            % Powell's trust region subproblem solver in UOBYQA.
            %--------------------------------------------------------------------------------------------------%


            i = NaN;
            j = NaN;

            Asubd = NaN;

            w = NaN(size(A, 1), 1);
            wz = NaN;
            z = NaN(size(A, 1), 1);

            n = size(A, 1);

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
                tdiag(:) = 0.0;
                tsubdiag(:) = 0.0;
                return
            elseif scaling > 1.0e8 || scaling < 1.0e-4
                % The thresholds are empirical.
                A(:, :) = A ./ scaling;
                scaled = true;
            end

            tdiag(:) = diag(A);

            for k = 1:n - 1
                colsq = sum(A(k + 2:n, k) .^ 2, 'all');
                if colsq <= 0
                    tsubdiag(k) = A(k + 1, k); % A(K+1, K) may have been updated in previous loops.
                    A(k + 1, k) = 0.0;
                    continue
                end

                Asubd = A(k + 1, k);
                tsubdiag(k) = sqrt(colsq + Asubd ^ 2) .* ((Asubd > 0) .* 2 - 1);

                A(k + 1, k) = -colsq / (Asubd + tsubdiag(k));
                w(k + 1:n) = sqrt(2.0 / (colsq + A(k + 1, k) ^ 2)) * A(k + 1:n, k);
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
                wz = sum(w(k + 1:n) .* z(k + 1:n), 'all');

                tdiag(k + 1:n) = tdiag(k + 1:n) + w(k + 1:n) .* (wz * w(k + 1:n) - 2.0 * z(k + 1:n));
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


        end
        function [H, Q] = hessenberg_full(obj, A, H, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine finds a Hessenberg matrix H (all entries below the subdiagonal are 0) such that
            % H = Q^T*A*Q, where Q is a orthogonal matrix that may also be returned. A will stay unchanged.
            %--------------------------------------------------------------------------------------------------%


            i = NaN;

            subd = NaN;
            v = NaN(size(A, 1), 1);
            w = NaN(size(A, 1), 1);

            % Debugging variables


            n = size(A, 1);

            %====================%
            % Calculation starts %
            %====================%


            ipObj = inputParser();
            addParameter(ipObj, 'Q', NaN);
            parse(ipObj, varargin{:});
            Q = ipObj.Results.Q;

            if n <= 0
                % Quick return when N <= 0. Of course, N < 0 is impossible.
                return
            end

            H(:, :) = A;
            if nargout >= 2
                Q = eye(n);
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
                subd = sqrt(v(j + 1) ^ 2 + colsq) .* ((v(j + 1) > 0) .* 2 - 1);

                %----------------------------------------------------------------------------------------------%
                v(j + 1) = -colsq / (v(j + 1) + subd);
                v(j + 1:n) = sqrt(2.0 / (colsq + v(j + 1) ^ 2)) * v(j + 1:n);
                % The two lines above are from Powell. They are equivalent to the following two lines.
                % %V(J + 1) = V(J + 1) - SUBD
                % %V(J + 1:N) = sqrt(TWO) * V(J + 1:N) / NORM(V(J + 1:N))
                %----------------------------------------------------------------------------------------------%

                for i = j + 1:n
                    H(j + 1:n, i) = H(j + 1:n, i) - sum(H(j + 1:n, i) .* v(j + 1:n), 'all') * v(j + 1:n);
                end
                H(j + 1, j) = subd;
                H(j + 2:n, j) = 0.0;

                w(:) = H(:, j + 1:n) * v(j + 1:n);
                for i = j + 1:n
                    H(:, i) = H(:, i) - w * v(i);
                end

                if nargout >= 2
                    w(:) = Q(:, j + 1:n) * v(j + 1:n);
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


            eig_min = NaN;

            k = NaN;

            piv = NaN(numel(td), 1);

            pivnew = NaN(numel(td), 1);

            n = numel(td);

            ipObj = inputParser();
            addParameter(ipObj, 'tol', NaN);
            parse(ipObj, varargin{:});
            tol = ipObj.Results.tol;

            %====================%
            % Calculation starts %
            %====================%

            maxiter = 100;
            tol_loc = 10.0 ^ max(-6, -308);
            if ~ismember('tol', ipObj.UsingDefaults)
                tol_loc = tol;
            end

            % The following loop calculates the Sturm ratios [Q_1(0), ..., Q_n(0)]. These ratios are all positive
            % iff all the eigenvalues of the matrix are positive definite. Note that these ratios are also the
            % pivots of the Cholesky factorization of the matrix (i.e., the square of the diagonal of L in LL^T,
            % or the diagonal of D in LDL^T). All the pivots are positive iff there exists a Cholesky
            % factorization with a positive diagonal, i.e., the matrix is positive definite.
            piv(:) = -1.0;
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
                eminlb = 0.0;
            else
                eminub = min(td, [], 'all');
                eminlb = -max(abs([0.0; tn]) + abs(td) + abs([tn; 0.0]), [], 'all');
            end

            ksav = 0;
            pivksv = 0.0; % This initial value will not be used, but Fortran compilers may complain without it.
            for iter = 1:maxiter                % Powell's code is essentially a DO WHILE loop. We impose an explicit MAXITER.
                if eminub - eminlb <= tol_loc * max(abs(eminlb), abs(eminub))
                    break
                end
                eig_min = 0.5 * (eminlb + eminub);

                % The following loop calculates the Sturm ratios [Q_1(EIG_MIN), ..., Q_n(EIG_MIN)]. These ratios
                % are all positive iff all the eigenvalues of the matrix are larger than EIG_MIN, i.e., EIG_MIN
                % underestimates the smallest eigenvalue. Note that these ratios are also the pivots of the
                % Cholesky factorization of the matrix minus EIG_MIN*I (i.e., the square of the diagonal of L in
                % LL^T, or the diagonal of D in LDL^T). All the pivots are positive iff there exists a Cholesky
                % factorization with a positive diagonal, i.e., the matrix minus LAMBDA*I is positive definite.
                pivnew(:) = -1.0;
                pivnew(1) = td(1) - eig_min;
                for k = 1:n - 1
                    if pivnew(k) > 0
                        pivnew(k + 1) = td(k + 1) - eig_min - tn(k) ^ 2 / pivnew(k);
                    else
                        break
                    end
                end

                if all(pivnew > 0, 'all')
                    piv = pivnew;
                    eminlb = eig_min;
                    continue
                end

                % We arrive here iff PIVNEW contains nonpositive entries and EIG_MIN is no less than the smallest
                % eigenvalue. We set EMINUB to EIG_MIN except a possible adjustment by the rule of false position.
                k = min(find(~(pivnew > 0)), [], 'all');
                piv(1:k - 1) = pivnew(1:k - 1);

                % KSAV was initialized to 0, triggering the ELSE when ALL(PIVNEW > 0) fails for the first time.
                if k == ksav && pivksv < 0 && piv(k) - pivnew(k) >= pivnew(k) - pivksv
                    pivksv = 0.0;
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
        function smat = vec2smat(~, vec)
            %--------------------------------------------------------------------------------------------------%
            % This function transforms a vector VEC to a symmetric matrix SMAT with the vector storing the upper
            % triangular part of the matrix column by column.
            %--------------------------------------------------------------------------------------------------%


            smat = NaN((round(sqrt(double(8 * numel(vec) + 1))) - 1) / 2);

            n = size(smat, 1);

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


        end
        function vec = smat2vec(~, smat)
            %--------------------------------------------------------------------------------------------------%
            % This function transforms a symmetric matrix SMAT to a vector VEC that stores the upper triangular
            % part of the matrix column by column.
            %--------------------------------------------------------------------------------------------------%


            vec = NaN((size(smat, 1) * (size(smat, 1) + 1)) / 2, 1);

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


            y = NaN(numel(x), 1);

            n = numel(x);

            %====================%
            % Calculation starts %
            %====================%

            for j = 1:n
                ih = (j - 1) * j / 2;
                y(j) = sum(smatv(ih + 1:ih + j) .* x(1:j), 'all');
                y(1:j - 1) = y(1:j - 1) + x(j) * smatv(ih + 1:ih + j - 1);
            end

            %====================%
            % Calculation ends   %
            %====================%

        end

    end
end
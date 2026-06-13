classdef powalg_mod
    %--------------------------------------------------------------------------------------------------%
    % This module provides some Powell-style linear algebra procedures.
    %
    % TODO:
    % Divide the module into three submodules:
    % - QR: procedures concerning QR factorization
    % - QUADRATIC: procedures concerning quadratic polynomials represented by [GQ, PQ, HQ] so that
    %   Q(Y) = <Y, GQ> + 0.5*<Y, HESSIAN*Y>,
    %   HESSIAN consists of an explicit part HQ and an implicit part PQ in Powell's way:
    %   HESSIAN = HQ + sum_K=1^NPT PQ(K)*(XPT(:, K)*XPT(:, K)^T) .
    % - LAGINT: procedures concerning quadratic LAGrange INTerpolation.
    %
    % Zaikun (20230321): In a test on 20230321 on problems of at most 200 variables, it affects (not
    % necessarily worsens) the performance of NEWUOA/LINCOA quite marginally if the update of IDZ is
    % completely disabled (and hence IDZ remains one forever). Therefore, in the first implementation of
    % an algorithm based on the derivative-free PSB, it seems same to ignore IDZ. It is similar for the
    % RESCUE technique of BOBYQA.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net).
    %
    % Started: July 2020
    %
    % Last Modified: Thu 14 Aug 2025 07:36:04 AM CST
    %--------------------------------------------------------------------------------------------------%

    methods
        function obj = powalg_mod()
            %--------------------------------------------------------------------------------------------------%
            % QR:

            %--------------------------------------------------------------------------------------------------%
            % QUADRATIC:


            %--------------------------------------------------------------------------------------------------%
            % LAGINT (quadratic LAGrange INTerpolation):

        end
        function varargout = qradd(obj, varargin)
            if numel(varargin) == 4 && isvector(varargin{3})
                [varargout{1:nargout}] = obj.qradd_Rdiag(varargin{:});
            else
                [varargout{1:nargout}] = obj.qradd_Rfull(varargin{:});
            end
        end
        function varargout = qrexc(obj, varargin)
            if numel(varargin) == 3 && isinteger(varargin{3}) && isscalar(varargin{3})
                [varargout{1:nargout}] = obj.qrexc_Rfull(varargin{:});
            else
                [varargout{1:nargout}] = obj.qrexc_Rdiag(varargin{:});
            end
        end
        function varargout = quadinc(obj, varargin)
            if numel(varargin) == 3 && isvector(varargin{2})
                [varargout{1:nargout}] = obj.quadinc_ghv(varargin{:});
            else
                [varargout{1:nargout}] = obj.quadinc_d0(varargin{:});
            end
        end
        function varargout = calvlag(obj, varargin)
            if numel(varargin) >= 5 && numel(varargin) <= 6 && isinteger(varargin{1}) && isscalar(varargin{1}) && (~isvector(varargin{2}) && ~isscalar(varargin{2})) && isfloat(varargin{4}) && (~isvector(varargin{4}) && ~isscalar(varargin{4}))
                [varargout{1:nargout}] = obj.calvlag_lfqint(varargin{:});
            else
                [varargout{1:nargout}] = obj.calvlag_qint(varargin{:});
            end
        end
        function [Q, Rdiag, n] = qradd_Rdiag(~, c, Q, Rdiag, n)            % Used in COBYLA
            %--------------------------------------------------------------------------------------------------%
            % This subroutine updates the QR factorization of an MxN matrix A of full column rank, attempting to
            % add a new column C is to this matrix as the LAST column while maintaining the full-rankness.
            % Case 1. If C is not in range(A) (theoretically, it implies N < M), then the new matrix is [A, C];
            % Case 2. If C is in range(A), then the new matrix is [A(:, 1:N-1), C].
            % Zaikun 2023903: It may happen in the second case that C is in the range of A(:, 1:N-1), and the
            % new matrix does not have full column rank any more. Indeed, Powell wrote in comments that "set
            % IOUT to the index of the constraint (here, column of A --- Zaikun) to be deleted, but branch if no
            % suitable index can be found". The idea is to replace a column of A by C so that the new matrix
            % still has full rank (such a column must exist unless C = 0). But his code sets IOUT = N always.
            % Maybe he found this worked well enough in practice. Meanwhile, Powell's code includes a snippet
            % that can never be reached, which was probably intended to deal with the case with IOUT =/= N.
            % N.B.:
            % 0. Instead of R, this subroutine updates RDIAG, which is diag(R), with a size at most M and at
            % least MIN(M, N+1). The number is MIN(M, N+1) rather than MIN(M, N) as N may be augmented by 1 in
            % the subroutine.
            % 1. With the two cases specified as above, this function does not need A as an input.
            % 2. The subroutine changes only Q(:, NSAVE+1:M) (NSAVE is the original value of N)
            % and R(:, N) (N takes the updated value).
            %--------------------------------------------------------------------------------------------------%

            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();
            linalg_obj = linalg_mod();

            % Inputs
            % C(M)

            % In-outputs

            % Q(M, M)
            % MIN(M, N+1) <= SIZE(Rdiag) <= M

            % Local variables
            srname = "QRADD_RDIAG";
            k = NaN;
            m = NaN;
            nsave = NaN;
            cq = NaN(size(Q, 2), 1);
            cqa = NaN(size(Q, 2), 1);
            G = NaN(2);
            %------------------------------------------------------------%
            Qsave = NaN(size(Q, 1), n); % Debugging only
            Rdsave = NaN(n, 1); % Debugging only
            tol = NaN; % Debugging only
            %------------------------------------------------------------%

            % Sizes
            m = size(Q, 2);

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(n >= 0 && n <= m, "0 <= N <= M", srname); % N = 0 is possible.
                debug_obj.assert(numel(c) == m, "SIZE(C) == M", srname);
                debug_obj.assert(numel(Rdiag) >= min(m, n + 1) && numel(Rdiag) <= m, "MIN(M, N+1) <= SIZE(Rdiag) <= M", srname);
                debug_obj.assert(size(Q, 1) == m && size(Q, 2) == m, "SIZE(Q) == [M, M]", srname);
                tol = max(consts_obj.TEN ^ max(-8, -consts_obj.MAXPOW10), min(0.1, consts_obj.TEN ^ min(12, consts_obj.MAXPOW10) * consts_obj.EPS * double(m + 1)));
                debug_obj.assert(linalg_obj.isorth(Q, 'tol', tol), "The columns of Q are orthonormal", srname); % Costly!
                Qsave(:, :) = Q(:, 1:n); % For debugging only
                Rdsave(:) = Rdiag(1:n); % For debugging only

            end

            %====================%
            % Calculation starts %
            %====================%

            nsave = n; % Needed for debugging (only).

            % As in Powell's COBYLA, CQ is set to 0 at the positions with CQ being negligible as per ISMINOR.
            % This may not be the best choice if the subroutine is used in other contexts, e.g., LINCOA.
            cq(:) = linalg_obj.matprod12(c, Q);
            cqa(:) = linalg_obj.matprod12(abs(c), abs(Q));
            cq(linalg_obj.trueloc(linalg_obj.isminor1(cq, cqa))) = consts_obj.ZERO; %%MATLAB: cq(isminor(cq, cqa)) = zero

            % Update Q so that the columns of Q(:, N+2:M) are orthogonal to C. This is done by applying a 2D
            % Givens rotation to Q(:, [K, K+1]) from the right to zero C'*Q(:, K+1) out for K = N+1, ..., M-1
            % in the reverse order. Nothing will be done if N >= M-1.
            for k = m - 1:-1:n + 1
                if abs(cq(k + 1)) > 0
                    % Powell wrote CQ(K+1) /= 0 instead of ABS(CQ(K+1)) > 0. The two differ if CQ(K+1) is NaN.
                    % If we apply the rotation below when CQ(K+1) = 0, then CQ(K) will get updated to |CQ(K)|.
                    G(:, :) = linalg_obj.planerot(cq([k, k + 1]));
                    Q(:, [k, k + 1]) = linalg_obj.matprod22(Q(:, [k, k + 1]), G');
                    cq(k) = linalg_obj.hypotenuse(cq(k), cq(k + 1)); %cq(k) = sqrt(cq(k)**2 + cq(k + 1)**2)

                end
            end

            % Augment N by 1 if C is not in range(A).
            % The two IFs cannot be merged as Fortran may evaluate CQ(N+1) even if N>=M, leading to a SEGFAULT.
            if n < m
                % Powell's condition for the following IF: CQ(N+1) /= 0.
                if abs(cq(n + 1)) > consts_obj.EPS ^ 2 && ~linalg_obj.isminor0(cq(n + 1), cqa(n + 1))
                    n = n + 1;
                end
            end

            % Update RDIAG so that RDIAG(N) = CQ(N) = INPROD(C, Q(:, N)). Note that N may have been augmented.
            % Zaikun 20230903: Different from QRADD_RFULL, Powell did not maintain the positiveness of RDIAG.
            if n >= 1 && n <= m
                % Indeed, N > M should not happen unless the input is wrong.
                Rdiag(n) = cq(n); % Indeed, RDIAG(N) = INPROD(C, Q(:, N))

            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(n >= nsave && n <= min(nsave + 1, m), "NSAV <= N <= MIN(NSAV + 1, M)", srname);
                debug_obj.assert(numel(Rdiag) >= n && numel(Rdiag) <= m, "N <= SIZE(Rdiag) <= M", srname);
                debug_obj.assert(size(Q, 1) == m && size(Q, 2) == m, "SIZE(Q) == [M, M]", srname);
                debug_obj.assert(linalg_obj.isorth(Q, 'tol', tol), "The columns of Q are orthonormal", srname); % Costly!

                debug_obj.assert(all(abs(Q(:, 1:nsave) - Qsave(:, 1:nsave)) <= 0, 'all'), "Q(:, 1:NSAVE) is unchanged", srname);
                debug_obj.assert(all(abs(Rdiag(1:n - 1) - Rdsave(1:n - 1)) <= 0, 'all'), "Rdiag(1:N-1) is unchanged", srname);

                if n < m && infnan_obj.is_finite(linalg_obj.p_norm(c))
                    debug_obj.assert(linalg_obj.p_norm(linalg_obj.matprod12(c, Q(:, n + 1:m))) <= max(tol, tol * linalg_obj.p_norm(c)), "C^T*Q(:, N+1:M) == 0", srname);
                end
                if n >= 1
                    % N = 0 is possible.
                    debug_obj.assert(abs(linalg_obj.inprod(c, Q(:, n)) - Rdiag(n)) <= max(tol, tol * linalg_obj.inprod(abs(c), abs(Q(:, n)))) || ~infnan_obj.is_finite(Rdiag(n)), "C^T*Q(:, N) == Rdiag(N)", srname);
                end
            end
        end
        function [Q, R, n] = qradd_Rfull(~, c, Q, R, n)            % Used in LINCOA
            %--------------------------------------------------------------------------------------------------%
            % This subroutine updates the QR factorization of an MxN matrix A = Q*R(:, 1:N) when a new column C
            % is appended to this matrix A as the LAST column.
            % N.B.:
            % 0. Different from QRADD_RDIAG, QRADD_RFULL always append C to A, and always increase N by 1. This
            % is because it is for sure that C is not in the column space of A in LINCOA.
            % 1. At entry, Q is a MxM orthonormal matrix, and R is a MxL upper triangular matrix with N < L <= M.
            % 2. The subroutine changes only Q(:, N+1:M) and R(:, N+1) with N taking the original value.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            linalg_obj = linalg_mod();

            % Inputs
            % C(M)

            % In-outputs

            % Q(M, M)
            % R(M, :), N+1 <= SIZE(R, 2) <= M

            % Local variables
            srname = "QRADD_RFULL";
            k = NaN;
            m = NaN;
            cq = NaN(size(Q, 2), 1);
            G = NaN(2);
            %------------------------------------------------------------%
            Anew = NaN(size(Q, 1), n + 1); % Debugging only
            Qsave = NaN(size(Q, 1), n); % Debugging only
            Rsave = NaN(size(R, 1), n); % Debugging only
            tol = NaN; % Debugging only
            %------------------------------------------------------------%

            % Sizes
            m = size(Q, 1);

            if consts_obj.DEBUGGING
                debug_obj.assert(n >= 0 && n <= m - 1, "0 <= N <= M - 1", srname);
                debug_obj.assert(numel(c) == m, "SIZE(C) == M", srname);
                debug_obj.assert(size(Q, 1) == m && size(Q, 2) == m, "SIZE(Q) = [M, M]", srname);
                debug_obj.assert(size(Q, 2) == size(R, 1), "SIZE(Q, 2) == SIZE(R, 1)", srname);
                debug_obj.assert(size(R, 2) >= n + 1 && size(R, 2) <= m, "N+1 <= SIZE(R, 2) <= M", srname);
                tol = max(consts_obj.TEN ^ max(-8, -consts_obj.MAXPOW10), min(0.1, consts_obj.TEN ^ min(8, consts_obj.MAXPOW10) * consts_obj.EPS * double(m + 1)));
                debug_obj.assert(linalg_obj.isorth(Q, 'tol', tol), "The columns of Q are orthogonal", srname);
                debug_obj.assert(linalg_obj.istriu(R), "R is upper triangular", srname);
                debug_obj.assert(all(linalg_obj.diag(R(:, 1:n)) > 0, 'all'), "DIAG(R(:, 1:N)) > 0", srname);
                Anew(:, :) = reshape([reshape(linalg_obj.matprod22(Q, R(:, 1:n)), 1, []), reshape(c, 1, [])], size(Anew));
                Qsave(:, :) = Q(:, 1:n); % For debugging only.
                Rsave(:, :) = R(:, 1:n); % For debugging only.

            end

            cq(:) = linalg_obj.matprod12(c, Q);

            % Update Q so that the columns of Q(:, N+2:M) are orthogonal to C. This is done by applying a 2D
            % Givens rotation to Q(:, [K, K+1]) from the right to zero C'*Q(:, K+1) out for K = N+1, ..., M-1.
            % Nothing will be done if N >= M-1.
            for k = m - 1:-1:n + 1
                if abs(cq(k + 1)) > 0
                    % Powell: IF (ABS(CQ(K + 1)) > 1.0D-20 * ABS(CQ(K))) THEN
                    G(:, :) = linalg_obj.planerot(cq([k, k + 1]));
                    Q(:, [k, k + 1]) = linalg_obj.matprod22(Q(:, [k, k + 1]), G');
                    cq(k) = fortran.sqrt(cq(k) ^ 2 + cq(k + 1) ^ 2);
                end
            end

            R(1:n, n + 1) = linalg_obj.matprod12(c, Q(:, 1:n));

            % Maintain the positiveness of the diagonal entries of R.
            if cq(n + 1) < 0
                Q(:, n + 1) = -Q(:, n + 1);
            end
            R(n + 1, n + 1) = abs(cq(n + 1));

            n = n + 1;

            if consts_obj.DEBUGGING
                debug_obj.assert(n >= 1 && n <= m, "1 <= N <= M", srname);
                debug_obj.assert(size(Q, 1) == m && size(Q, 2) == m, "SIZE(Q) = [M, M]", srname);
                debug_obj.assert(size(Q, 2) == size(R, 1), "SIZE(Q, 2) == SIZE(R, 1)", srname);
                debug_obj.assert(size(R, 2) >= n && size(R, 2) <= m, "N <= SIZE(R, 2) <= M", srname);
                debug_obj.assert(linalg_obj.isorth(Q, 'tol', tol), "The columns of Q are orthogonal", srname);
                debug_obj.assert(linalg_obj.istriu(R), "R is upper triangular", srname);
                debug_obj.assert(all(linalg_obj.diag(R(:, 1:n)) > 0, 'all'), "DIAG(R(:, 1:N)) > 0", srname);

                % %call assert(.not. any(abs(Q(:, 1:n - 1) - Qsave(:, 1:n - 1)) > 0), 'Q(:, 1:N-1) is unchanged', srname)
                % %call assert(.not. any(abs(R(:, 1:n - 1) - Rsave(:, 1:n - 1)) > 0), 'R(:, 1:N-1) is unchanged', srname)
                % If we can ensure that Q and R do not contain NaN or Inf, use the following lines instead of the last two.
                debug_obj.assert(all(abs(Q(:, 1:n - 1) - Qsave(:, 1:n - 1)) <= 0, 'all'), "Q(:, 1:N-1) is unchanged", srname);
                debug_obj.assert(all(abs(R(:, 1:n - 1) - Rsave(:, 1:n - 1)) <= 0, 'all'), "R(:, 1:N-1) is unchanged", srname);

                % The following test may fail.
                debug_obj.assert(all(abs(Anew - linalg_obj.matprod22(Q, R(:, 1:n))) <= max(tol, tol * max(abs(Anew), [], 'all')), 'all'), "Anew = Q*R", srname);
            end
        end
        function [Q, Rdiag] = qrexc_Rdiag(~, A, Q, Rdiag, i)            % Used in COBYLA
            %--------------------------------------------------------------------------------------------------%
            % This subroutine updates the QR factorization for an MxN matrix A = Q*R so that the updated Q and
            % R form a QR factorization of [A_1, ..., A_{I-1}, A_{I+1}, ..., A_N, A_I], which is the matrix
            % obtained by rearranging columns [I, I+1, ..., N] of A to [I+1, ..., N, I]. Here, Q is a matrix
            % whose columns are orthogonal, and R, which is not present, is an upper triangular matrix whose
            % diagonal entries are nonzero. Q and R need not to be square.
            % N.B.:
            % 0. Instead of R, this subroutine updates RDIAG, which is diag(R), the size being N.
            % 1. With L = SIZE(Q, 2) = SIZE(R, 1), we have M >= L >= N. Most often, L = M or N.
            % 2. The subroutine changes only Q(:, I:N) and RDIAG(I:N).
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            linalg_obj = linalg_mod();

            % Inputs
            % A(M, N)

            % In-outputs
            % Q(M, :), N <= SIZE(Q, 2) <= M
            % Rdiag(N)


            % Local variables
            srname = "QREXC_RDIAG";
            k = NaN;
            m = NaN;
            n = NaN;
            G = NaN(2);
            %------------------------------------------------------------%
            Anew = NaN(size(A, 1), size(A, 2)); % Debugging only
            Qsave = NaN(size(Q, 1), size(Q, 2)); % Debugging only
            QtAnew = NaN(size(Q, 2), size(A, 2)); % Debugging only
            Rdsave = NaN(i, 1); % Debugging only
            tol = NaN; % Debugging only
            %------------------------------------------------------------%

            % Sizes
            m = size(A, 1);
            n = size(A, 2);

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(n >= 1 && n <= m, "1 <= N <= M", srname);
                debug_obj.assert(i >= 1 && i <= n, "1 <= i <= N", srname);
                debug_obj.assert(numel(Rdiag) == n, "SIZE(Rdiag) == N", srname);
                debug_obj.assert(size(Q, 1) == m && size(Q, 2) >= n && size(Q, 2) <= m, "SIZE(Q, 1) == M, N <= SIZE(Q, 2) <= M", srname);
                tol = max(consts_obj.TEN ^ max(-8, -consts_obj.MAXPOW10), min(0.1, consts_obj.TEN ^ min(8, consts_obj.MAXPOW10) * consts_obj.EPS * double(m + 1)));
                debug_obj.assert(linalg_obj.isorth(Q, 'tol', tol), "The columns of Q are orthonormal", srname); % Costly!
                Qsave(:, :) = Q; % For debugging only.
                Rdsave(:) = Rdiag(1:i); % For debugging only.

            end

            %====================%
            % Calculation starts %
            %====================%

            if i <= 0 || i >= n
                % Only I == N is really needed, as 1 <= I <= N unless the input is wrong.
                return
            end

            % Let R be the upper triangular matrix in the QR factorization, namely R = Q^T*A.
            % For each K, find the Givens rotation G with G*R([K, K+1], :) = [HYPT, 0], and update Q(:, [K,K+1])
            % to Q(:, [K, K+1])*G^T. Then R = Q^T*A is an upper triangular matrix as long as A(:, [K, K+1]) is
            % updated to A(:, [K+1, K]). Indeed, this new upper triangular matrix can be obtained by first
            % updating R([K, K+1], :) to G*R([K, K+1], :) and then exchanging its columns K and K+1; at the same
            % time, entries K and K+1 of R's diagonal RDIAG become [HYPT, -(RDIAG(K+1) / HYPT) * RDIAG(K)].
            % After this is done for each K = 1, ..., N-1, we obtain the QR factorization of the matrix that
            % rearranges columns [I, I+1, ..., N] of A as [I+1, ..., N, I].
            % Powell's code, however, is slightly different: before everything, he first exchanged columns K and
            % K+1 of Q (as well as rows K and K+1 of R). This makes sure that the entires of the update RDIAG
            % are all positive if it is the case for the original RDIAG.
            % Zaikun 20230903: It turns out that Powell's code does not ensure that the original RDIAG is
            % positive (see QRADD_RDIAG), and hence the updated RDIAG may contain negative values.
            for k = i:n - 1
                G(:, :) = linalg_obj.planerot([Rdiag(k + 1), linalg_obj.inprod(Q(:, k), A(:, k + 1))]);
                Q(:, [k, k + 1]) = linalg_obj.matprod22(Q(:, [k + 1, k]), G');
                % Powell's code updates RDIAG in the following way:
                % %HYPT = SQRT(RDIAG(K + 1)**2 + INPROD(Q(:, K), A(:, K + 1))**2)
                % %RDIAG([K, K + 1_IK]) = [HYPT, (RDIAG(K + 1) / HYPT) * RDIAG(K)]
                % Note that RDIAG(N) inherits all rounding in RDIAG(I:N-1) and Q(:, I:N-1) and hence contain
                % significant errors. Thus we may modify Powell's code to set only RDIAG(K) = HYPT here and then
                % calculate RDIAG(N) by an inner product after the loop. Nevertheless, we simply calculate RDIAG
                % from scratch we do below.
            end

            % Calculate RDIAG(I:N) from scratch.
            Rdiag(i:n - 1) = reshape(cell2mat(arrayfun(@(k) linalg_obj.inprod(Q(:, k), A(:, k + 1)), (i:n - 1), "UniformOutput", false)), [], 1);
            %%MATLAB: Rdiag(i:n-1) = sum(Q(:, i:n-1) .* A(:, i+1:n), 1);  % Row vector
            Rdiag(n) = linalg_obj.inprod(Q(:, n), A(:, i)); % Calculate RDIAG(N) from scratch. See the comments above.

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(Rdiag) == n, "SIZE(Rdiag) == N", srname);
                debug_obj.assert(size(Q, 1) == m && size(Q, 2) >= n && size(Q, 2) <= m, "SIZE(Q, 1) == M, N <= SIZE(Q, 2) <= M", srname);
                debug_obj.assert(linalg_obj.isorth(Q, 'tol', tol), "The columns of Q are orthonormal", srname); % Costly!

                Qsave(:, i:n) = Q(:, i:n);
                debug_obj.assert(all(abs(Q - Qsave) <= 0, 'all'), "Q is unchanged except Q(:, I:N)", srname);
                debug_obj.assert(all(abs(Rdiag(1:i - 1) - Rdsave(1:i - 1)) <= 0, 'all'), "Rdiag(1:I-1) is unchanged", srname);

                Anew(:, :) = reshape([reshape(A(:, 1:i - 1), 1, []), reshape(A(:, i + 1:n), 1, []), reshape(A(:, i), 1, [])], size(Anew));
                QtAnew(:, :) = linalg_obj.matprod22(Q', Anew);
                debug_obj.assert(linalg_obj.istriu(QtAnew, 'tol', tol), "Q^T*Anew is upper triangular", srname);
                % The following test may fail if RDIAG is not calculated from scratch.
                debug_obj.assert(linalg_obj.p_norm(linalg_obj.diag(QtAnew) - Rdiag) <= max(tol, tol * linalg_obj.p_norm(reshape(cell2mat(arrayfun(@(k) linalg_obj.inprod(abs(Q(:, k)), abs(Anew(:, k))), (1:n), "UniformOutput", false)), [], 1))), "Rdiag == diag(Q^T*Anew)", srname);
                %%MATLAB: norm(diag(QtAnew) - Rdiag) <= max(tol, tol * norm(sum(abs(Q(:, 1:n)) .* abs(Anew), 1)))

            end
        end
        function [Q, R] = qrexc_Rfull(~, Q, R, i)            % Used in LINCOA
            %--------------------------------------------------------------------------------------------------%
            % This subroutine updates the QR factorization for an MxN matrix A = Q*R so that the updated Q and
            % R form a QR factorization of [A_1, ..., A_{I-1}, A_{I+1}, ..., A_N, A_I], which is the matrix
            % obtained by rearranging columns [I, I+1, ..., N] of A to [I+1, ..., N, I]. At entry, A = Q*R,
            % Q is a matrix whose columns are orthogonal, and R is an upper triangular matrix whose diagonal
            % entries are all nonzero. Q and R need not to be square.
            % N.B.:
            % 1. With L = SIZE(Q, 2) = SIZE(R, 1), we have M >= L >= N. Most often, L = M or N.
            % 2. The subroutine changes only Q(:, I:N) and R(:, I:N).
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            linalg_obj = linalg_mod();

            % Inputs


            % In-outputs
            % Q(M, :), SIZE(Q, 2) <= M
            % R(:, N), SIZE(R, 1) >= N

            % Local variables
            srname = "QREXC_RFULL";
            k = NaN;
            m = NaN;
            n = NaN;
            G = NaN(2);
            hypt = NaN;
            %------------------------------------------------------------%
            Anew = NaN(size(Q, 1), size(R, 2)); % Debugging only
            Qsave = NaN(size(Q, 1), size(Q, 2)); % Debugging only
            Rsave = NaN(size(R, 1), i); % Debugging only
            tol = NaN; % Debugging only
            %------------------------------------------------------------%

            % Sizes
            m = size(Q, 1);
            n = size(R, 2);

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(n >= 1 && n <= m, "1 <= N <= M", srname);
                debug_obj.assert(i >= 1 && i <= n, "1 <= I <= N", srname);
                debug_obj.assert(size(Q, 2) == size(R, 1), "SIZE(Q, 2) == SIZE(R, 1)", srname);
                debug_obj.assert(size(Q, 2) >= n && size(Q, 2) <= m, "N <= SIZE(Q, 2) <= M", srname);
                debug_obj.assert(size(R, 1) >= n && size(R, 1) <= m, "N <= SIZE(R, 1) <= M", srname);
                tol = max(consts_obj.TEN ^ max(-8, -consts_obj.MAXPOW10), min(0.1, consts_obj.TEN ^ min(8, consts_obj.MAXPOW10) * consts_obj.EPS * double(m + 1)));
                debug_obj.assert(linalg_obj.isorth(Q, 'tol', tol), "The columns of Q are orthogonal", srname);
                debug_obj.assert(linalg_obj.istriu(R), "R is upper triangular", srname);
                debug_obj.assert(all(linalg_obj.diag(R(:, 1:n)) > 0, 'all'), "DIAG(R(:, 1:N)) > 0", srname);
                Anew(:, :) = linalg_obj.matprod22(Q, R);
                Anew(:, :) = reshape([reshape(Anew(:, 1:i - 1), 1, []), reshape(Anew(:, i + 1:n), 1, []), reshape(Anew(:, i), 1, [])], size(Anew));
                Qsave(:, :) = Q; % For debugging only.
                Rsave(:, :) = R(:, 1:i); % For debugging only.

            end

            %====================%
            % Calculation starts %
            %====================%

            if i <= 0 || i >= n
                % Only I == N is really needed, as 1 <= I <= N unless the input is wrong.
                return
            end

            % For each K, find the Givens rotation G with G*R([K, K+1], K+1) = [HYPT, 0]. Then make two updates.
            % First, update Q(:, [K, K+1]) to Q(:, [K, K+1])*G^T, and R([K, K+1], :) to G*R[K+1, K], :), which
            % keeps Q*R unchanged and maintains the orthogonality of Q's columns. Second, exchange columns K and
            % K+1 of R. Then R becomes upper triangular, and the new product Q*R exchanges columns K and K+1 of
            % the original one. After this is done for each K = 1, ..., N-1, we obtain the QR factorization of
            % the matrix that rearranges columns [I, I+1, ..., N] of A as [I+1, ..., N, I].
            % Powell's code, however, is slightly different: before everything, he first exchanged columns K and
            % K+1 of Q as well as rows K and K+1 of R. This makes sure that the diagonal entries of the updated
            % R are all positive if it is the case for the original R.
            for k = i:n - 1
                G(:, :) = linalg_obj.planerot(R([k + 1, k], k + 1));
                % HYPT must be calculated before R is updated.
                hypt = linalg_obj.hypotenuse(R(k + 1, k + 1), R(k, k + 1)); %hypt = sqrt(R(k, k + 1)**2 + R(k + 1, k + 1)**2)

                % Update Q(:, [K, K+1]).
                Q(:, [k, k + 1]) = linalg_obj.matprod22(Q(:, [k + 1, k]), G');

                % Update R([K, K+1], :).
                R([k, k + 1], k:n) = linalg_obj.matprod22(G, R([k + 1, k], k:n));
                R(1:k + 1, [k, k + 1]) = R(1:k + 1, [k + 1, k]);
                % N.B.: The above two lines implement the following while noting that R is upper triangular.
                % %R([K, K + 1_IK], :) = MATPROD(G, R([K + 1_IK, K], :))  ! No need for R([K, K+1], 1:K-1) = 0
                % %R(:, [K, K + 1_IK]) = R(:, [K + 1_IK, K])  ! No need for R(K+2:, [K, K+1]) = 0

                % Revise R([K, K+1], K). Changes nothing in theory but seems good for the practical performance.
                R([k, k + 1], k) = [hypt, consts_obj.ZERO];

                %----------------------------------------------------------------------------------------------%
                % The following code performs the update without exchanging columns K and K+1 of Q or rows K and
                % K+1 of R beforehand. If the diagonal entries of the original R are positive, then all the
                % updated ones become negative.
                %
                % %G = planerot(R([k, k + 1_IK], k + 1))
                % %hypt = hypotenuse(R(k + 1, k + 1), R(k, k + 1)) !hypt = sqrt(R(k, k + 1)**2 + R(k + 1, k + 1)**2)
                % %
                % %Q(:, [k, k + 1_IK]) = matprod(Q(:, [k, k + 1_IK]), transpose(G))
                % %
                % %R([k, k + 1_IK], k:n) = matprod(G, R([k, k + 1_IK], k:n))
                % %R(1:k + 1, [k, k + 1_IK]) = R(1:k + 1, [k + 1_IK, k])
                % %R([k, k + 1_IK], k) = [hypt, ZERO]
                %----------------------------------------------------------------------------------------------%
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(size(Q, 2) == size(R, 1), "SIZE(Q, 2) == SIZE(R, 1)", srname);
                debug_obj.assert(size(Q, 2) >= n && size(Q, 2) <= m, "N <= SIZE(Q, 2) <= M", srname);
                debug_obj.assert(size(R, 1) >= n && size(R, 1) <= m, "N <= SIZE(R, 1) <= M", srname);
                debug_obj.assert(linalg_obj.isorth(Q, 'tol', tol), "The columns of Q are orthogonal", srname);
                debug_obj.assert(linalg_obj.istriu(R), "R is upper triangular", srname);
                debug_obj.assert(all(linalg_obj.diag(R(:, 1:n)) > 0, 'all'), "DIAG(R(:, 1:N)) > 0", srname);

                Qsave(:, i:n) = Q(:, i:n);
                % %call assert(.not. any(abs(Q - Qsave) > 0), 'Q is unchanged except Q(:, I:N)', srname)
                % %call assert(.not. any(abs(R(:, 1:i - 1) - Rsave(:, 1:i - 1)) > 0), 'R(:, 1:I-1) is unchanged', srname)
                % If we can ensure that Q and R do not contain NaN or Inf, use the following lines instead of the last two.
                debug_obj.assert(all(abs(Q - Qsave) <= 0, 'all'), "Q is unchanged except Q(:, I:N)", srname);
                debug_obj.assert(all(abs(R(:, 1:i - 1) - Rsave(:, 1:i - 1)) <= 0, 'all'), "R(:, 1:I-1) is unchanged", srname);

                % The following test may fail.
                debug_obj.assert(all(abs(Anew - linalg_obj.matprod22(Q, R)) <= max(tol, tol * max(abs(Anew), [], 'all')), 'all'), "Anew = Q*R", srname);
            end

        end
        function qinc = quadinc_d0(~, d, xpt, gq, pq, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This function evaluates QINC = Q(D) - Q(0) with Q being the quadratic function defined
            % via [GQ, HQ, PQ] by
            % Q(Y) = <Y, GQ> + 0.5*<Y, HESSIAN*Y>,
            % where HESSIAN consists of an explicit part HQ and an implicit part PQ in Powell's way:
            % HESSIAN = HQ + sum_K=1^NPT PQ(K)*(XPT(:, K)*XPT(:, K)^T) .
            % N.B.: QUADINC_D0(D, XPT, GQ, PQ, HQ) = QUADINC_DX(D, ZEROS(SIZE(D)), XPT, GQ, PQ, HQ)
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();
            linalg_obj = linalg_mod();

            % Inputs
            % D(N)
            % XPT(N, NPT)
            % GQ(N)
            % PQ(NPT)
            % HQ(N, N)

            % Output
            qinc = NaN;

            % Local variable
            srname = "QUADINC_D0";
            n = NaN;
            npt = NaN;
            dxpt = NaN(numel(pq), 1);

            % Sizes
            n = size(xpt, 1);
            npt = size(xpt, 2);

            % Preconditions
            ipObj = inputParser();
            addParameter(ipObj, 'hq', NaN);
            parse(ipObj, varargin{:});
            hq = ipObj.Results.hq;
            if consts_obj.DEBUGGING
                debug_obj.assert(n >= 1, "N >= 1", srname);
                debug_obj.assert(numel(d) == n && all(infnan_obj.is_finite(d), 'all'), "SIZE(D) == N, D is finite", srname);
                debug_obj.assert(all(infnan_obj.is_finite(xpt), 'all'), "XPT is finite", srname);
                debug_obj.assert(numel(gq) == n, "SIZE(GQ) = N", srname);
                debug_obj.assert(numel(pq) == npt, "SIZE(PQ) = NPT", srname);
                if ~ismember('hq', ipObj.UsingDefaults)
                    debug_obj.assert(size(hq, 1) == n && linalg_obj.issymmetric(hq), "HQ is an NxN symmetric matrix", srname);
                end
            end

            %====================%
            % Calculation starts %
            %====================%

            %--------------------------------------------------------------------------------------------------%
            % The following is Powell's scheme in LINCOA.
            % %% First-order term and explicit second-order term
            % %qinc = ZERO
            % %do j = 1, n
            % %    qinc = qinc + d(j) * gq(j)
            % %    do i = 1, j
            % %        t = d(i) * d(j)
            % %        if (i == j) then
            % %            t = HALF * t
            % %        end if
            % %        if (present(hq)) then
            % %            qinc = qinc + t * hq(i, j)
            % %        end if
            % %    end do
            % %end do
            % %
            % %% Implicit second-order term
            % %dxpt = matprod(d, xpt)
            % %do i = 1, npt
            % %    qinc = qinc + HALF * pq(i) * dxpt(i) * dxpt(i)  ! In BOBYQA, it is QINC - HALF * PQ(I) * DXPT(I)**2.
            % %end do
            %--------------------------------------------------------------------------------------------------%

            %--------------------------------------------------------------------------------------------------%
            % The following is a loop-free implementation, which should be applied in MATLAB/Python/R/Julia.
            % N.B.: INPROD(DXPT, PQ * DXPT) = INPROD(D, HESS_MUL(D, XPT, PQ))
            %--------------------------------------------------------------------------------------------------%
            dxpt(:) = linalg_obj.matprod12(d, xpt);
            if ismember('hq', ipObj.UsingDefaults)
                qinc = linalg_obj.inprod(d, gq) + consts_obj.HALF * linalg_obj.inprod(dxpt, pq .* dxpt);
            else
                qinc = linalg_obj.inprod(d, gq + consts_obj.HALF * linalg_obj.matprod21(hq, d)) + consts_obj.HALF * linalg_obj.inprod(dxpt, pq .* dxpt);
            end
            %%MATLAB:
            %%if nargin >= 5
            %%    qinc = d'*(gq + 0.5*hq*d) + 0.5*dxpt'*(pq*dxpt);
            %%else
            %%    qinc = d'*gq + 0.5*dxpt'*(pq*dxpt);
            %%end
            %--------------------------------------------------------------------------------------------------%

            %====================%
            %  Calculation ends  %
            %====================%

        end
        function qinc = quadinc_ghv(~, ghv, d, x)
            %--------------------------------------------------------------------------------------------------%
            % This function evaluates QINC = Q(X+D) - Q(X) with Q being the quadratic function defined via GHV by
            % Q(Y) = <Y, GQ> + 0.5*<Y, HESSIAN*Y>,
            % where GQ is GHV(1:N), and HESSIAN is the symmetric matrix whose upper triangular part is stored in
            % GHV(N+1:N*(N+3)/2) column by column.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            linalg_obj = linalg_mod();
            % Inputs



            % Outputs
            qinc = NaN;
            % Local variables
            srname = "QUADINC_GHV";
            ih = NaN;
            n = NaN;
            j = NaN;
            s = NaN(numel(x), 1);
            w = NaN(numel(ghv), 1);

            % Sizes
            n = fix(numel(x));

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(d) == n, "SIZE(D) = N", srname);
                debug_obj.assert(numel(ghv) == n * (n + 3) / 2, "SIZE(GHV) = N*(N+3)/2", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            s(:) = x + d;

            w(1:n) = d;
            for j = 1:n
                ih = n + (j - 1) * j / 2;
                w(ih + 1:ih + j) = d(1:j) * s(j) + d(j) * x(1:j);
                w(ih + j) = consts_obj.HALF * w(ih + j);
            end

            qinc = linalg_obj.inprod(ghv, w);

            %====================%
            % Calculation ends   %
            %====================%

        end
        function err = errquad(obj, fval, xpt, gq, pq, hq, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This function calculates the maximal relative error of Q in interpolating FVAL on XPT.
            % Here, Q is the quadratic function defined via [GQ, HQ, PQ] by
            % Q(Y) = <Y, GQ> + 0.5*<Y, HESSIAN*Y> if KREF is absent,
            % Q(Y) = <Y-XREF, GQ> + 0.5*<Y-XREF, HESSIAN*(Y-XREF)> with XREF = XPT(:, KREF), if KREF is present.
            % Here, HESSIAN consists of an explicit part HQ and an implicit part PQ in Powell's way:
            % HESSIAN = HQ + sum_K=1^NPT PQ(K)*(XPT(:, K)*XPT(:, K)^T).
            % N.B.: If KREF is absent, then GQ = nabla Q(0); otherwise, GQ = nabla Q(XREF).
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();
            linalg_obj = linalg_mod();

            % Inputs
            % FVAL(NPT)
            % XPT(N, NPT)
            % GQ(N)
            % PQ(NPT)
            % HQ(N, N)


            % Outputs
            err = NaN;

            % Local variables
            srname = "ERRQUAD";
            k = NaN;
            n = NaN;
            npt = NaN;
            fmq = NaN(size(xpt, 2), 1);
            qval = NaN(size(xpt, 2), 1);

            % Sizes
            n = size(xpt, 1);
            npt = size(xpt, 2);

            % Preconditions
            ipObj = inputParser();
            addParameter(ipObj, 'kref', NaN);
            parse(ipObj, varargin{:});
            kref = ipObj.Results.kref;
            if consts_obj.DEBUGGING
                debug_obj.assert(n >= 1, "N >= 1", srname);
                debug_obj.assert(numel(fval) == npt, "SIZE(FVAL) == NPT", srname);
                debug_obj.assert(~any(infnan_obj.is_nan(fval) | infnan_obj.is_posinf(fval), 'all'), "FVAL is not NaN/+Inf", srname);
                debug_obj.assert(all(infnan_obj.is_finite(xpt), 'all'), "XPT is finite", srname);
                debug_obj.assert(numel(gq) == n, "SIZE(GQ) == N", srname);
                debug_obj.assert(numel(pq) == npt, "SIZE(PQ) == NPT", srname);
                debug_obj.assert(size(hq, 1) == n && linalg_obj.issymmetric(hq), "HQ is an NxN symmetric matrix", srname);
                if ~ismember('kref', ipObj.UsingDefaults)
                    debug_obj.assert(kref >= 1 && kref <= npt, "1 <= KREF <= NPT", srname);
                end
            end

            %====================%
            % Calculation starts %
            %====================%

            if ismember('kref', ipObj.UsingDefaults)
                qval(:) = reshape(cell2mat(arrayfun(@(k) obj.quadinc_d0(xpt(:, k), xpt, gq, pq, 'hq', hq), (1:npt), "UniformOutput", false)), [], 1);
            else
                qval(:) = reshape(cell2mat(arrayfun(@(k) obj.quadinc_d0(xpt(:, k) - xpt(:, kref), xpt, gq, pq, 'hq', hq), (1:npt), "UniformOutput", false)), [], 1);
            end
            %%MATLAB:
            %%if nargin >= 5
            %%    qval = cellfun(@(x) quadinc(x, xpt, gq, pq, hq), num2cell(xpt - xpt(:, kref), 1));  % Row vector
            %%    % xpt - xpt(:, kref): Implicit expansion
            %%else
            %%    qval = cellfun(@(x) quadinc(x, xpt, gq, pq, hq), num2cell(xpt, 1));  % Row vector
            %%end
            if all(infnan_obj.is_finite(qval), 'all')
                fmq(:) = fval - qval;
                err = (max(fmq, [], 'all') - min(fmq, [], 'all')) / max([consts_obj.ONE; reshape(abs(fval), [], 1)], [], 'all');
            else
                err = consts_obj.REALMAX;
            end

            %====================%
            %  Calculation ends  %
            %====================%
            %
        end
        function y = hess_mul(~, x, xpt, pq, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This function calculates HESSIAN*X, with HESSIAN consisting of an explicit part HQ (0 if absent)
            % and an implicit part PQ in Powell's way: HESSIAN = HQ + sum_K=1^NPT PQ(K)*(XPT(:, K)*XPT(:, K)^T).
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();
            linalg_obj = linalg_mod();

            % Inputs
            % X(N)
            % XPT(N, NPT)
            % PQ(NPT)
            % HQ(N, N)

            % Outputs
            y = NaN(numel(x), 1);

            % Local variables
            srname = "HESS_MUL";
            j = NaN;
            n = NaN;
            npt = NaN;

            % Sizes
            n = size(xpt, 1);
            npt = size(xpt, 2);

            % Preconditions
            ipObj = inputParser();
            addParameter(ipObj, 'hq', NaN);
            parse(ipObj, varargin{:});
            hq = ipObj.Results.hq;
            if consts_obj.DEBUGGING
                debug_obj.assert(n >= 1, "N >= 1", srname);
                debug_obj.assert(numel(x) == n, "SIZE(Y) == N", srname);
                debug_obj.assert(all(infnan_obj.is_finite(xpt), 'all'), "XPT is finite", srname);
                debug_obj.assert(numel(pq) == npt, "SIZE(PQ) == NPT", srname);
                if ~ismember('hq', ipObj.UsingDefaults)
                    debug_obj.assert(size(hq, 1) == n && linalg_obj.issymmetric(hq), "HQ is an NxN symmetric matrix", srname);
                end
            end

            %====================%
            % Calculation starts %
            %====================%

            %--------------------------------------------------------------------------------%
            %----------! y = matprod(hq, x) + matprod(xpt, pq * matprod(x, xpt)) !-----------%
            %--------------------------------------------------------------------------------%
            y(:) = linalg_obj.matprod21(xpt, pq .* linalg_obj.matprod12(x, xpt));
            if ~ismember('hq', ipObj.UsingDefaults)
                for j = 1:n
                    y(:) = y + hq(:, j) * x(j);
                end
            end

            %====================%
            %  Calculation ends  %
            %====================%

        end
        function y = omega_col(~, idz, zmat, k)
            %--------------------------------------------------------------------------------------------------%
            % This function calculates Y = column K of OMEGA. As Powell did in NEWUOA, BOBYQA, and LINCOA,
            % OMEGA = sum_{i=1}^{K} S_i*ZMAT(:, i)*ZMAT(:, i)^T if S_i = -1 when i < IDZ and S_i = 1 if i >= IDZ
            % OMEGA is the leading NPT-by-NPT block of the matrix H in (3.12) of the NEWUOA paper.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            linalg_obj = linalg_mod();

            % Inputs



            % Outputs
            y = NaN(size(zmat, 1), 1);

            % Local variables
            srname = "OMEGA_COL";
            zk = NaN(size(zmat, 2), 1);

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(idz >= 1 && idz <= size(zmat, 2) + 1, "1 <= IDZ <= SIZE(ZMAT, 2) + 1", srname);
                debug_obj.assert(k >= 1 && idz <= size(zmat, 1), "1 <= K <= SIZE(ZMAT, 1)", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            zk(:) = zmat(k, :);
            zk(1:idz - 1) = -zk(1:idz - 1);
            y(:) = linalg_obj.matprod21(zmat, zk);

            %====================%
            %  Calculation ends  %
            %====================%

        end
        function y = omega_mul(~, idz, zmat, x)
            %--------------------------------------------------------------------------------------------------%
            % This function calculates Y = OMEGA*X. As Powell did in NEWUOA, BOBYQA, and LINCOA,
            % OMEGA = sum_{i=1}^{K} S_i*ZMAT(:, i)*ZMAT(:, i)^T if S_i = -1 when i < IDZ and S_i = 1 if i >= IDZ
            % OMEGA is the leading NPT-by-NPT block of the matrix H in (3.12) of the NEWUOA paper.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            linalg_obj = linalg_mod();

            % Inputs



            % Outputs
            y = NaN(size(zmat, 1), 1);

            % Local variables
            srname = "OMEGA_MUL";
            xz = NaN(size(zmat, 2), 1);

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(idz >= 1 && idz <= size(zmat, 2) + 1, "1 <= IDZ <= SIZE(ZMAT, 2) + 1", srname);
                debug_obj.assert(numel(x) == size(zmat, 1), "SIZE(X) == SIZE(ZMAT, 1)", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            xz(:) = linalg_obj.matprod12(x, zmat);
            xz(1:idz - 1) = -xz(1:idz - 1);
            y(:) = linalg_obj.matprod21(zmat, xz);

            %====================%
            %  Calculation ends  %
            %====================%

        end
        function p = omega_inprod(~, idz, zmat, x, y)
            %--------------------------------------------------------------------------------------------------%
            % This function calculates P = X^T*OMEGA*Y. As Powell did in NEWUOA, BOBYQA, and LINCOA,
            % OMEGA = sum_{i=1}^{K} S_i*ZMAT(:, i)*ZMAT(:, i)^T if S_i = -1 when i < IDZ and S_i = 1 if i >= IDZ
            % OMEGA is the leading NPT-by-NPT block of the matrix H in (3.12) of the NEWUOA paper.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            linalg_obj = linalg_mod();

            % Inputs



            % Outputs
            p = NaN;

            % Local variables
            srname = "OMEGA_INPROD";
            xz = NaN(size(zmat, 2), 1);
            yz = NaN(size(zmat, 2), 1);

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(idz >= 1 && idz <= size(zmat, 2) + 1, "1 <= IDZ <= SIZE(ZMAT, 2) + 1", srname);
                debug_obj.assert(numel(x) == size(zmat, 1), "SIZE(X) == SIZE(ZMAT, 1)", srname);
                debug_obj.assert(numel(y) == size(zmat, 1), "SIZE(Y) == SIZE(ZMAT, 1)", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            xz(:) = linalg_obj.matprod12(x, zmat);
            xz(1:idz - 1) = -xz(1:idz - 1);
            yz(:) = linalg_obj.matprod12(y, zmat);
            p = linalg_obj.inprod(xz, yz);

            %====================%
            %  Calculation ends  %
            %====================%

        end
        function err = errh(~, idz, bmat, zmat, xpt)
            %--------------------------------------------------------------------------------------------------%
            % This function calculates the error in H as the inverse of W. See (3.12) of the NEWUOA paper.
            % N.B.: The (NPT+1)th column (row) of H is not contained in [BMAT, ZMAT]. It is [r; t(1); s] below.
            % In the complete form, using MATLAB-style notation, the W and H in the NEWUOA paper are as follows.
            % W = [A, ONES(NPT, 1), XPT^T; ONES(1, NPT), ZERO, ZEROS(1, N); XPT, ZEROS(N, 1), ZEROS(N, N)]
            % H = [Omega, r, BMAT(:, 1:NPT)^T; r^T, t(1), s^T, BMAT(:, 1:NPT), s, BMAT(:, NPT+1:NPT+N)]
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();
            linalg_obj = linalg_mod();

            % Inputs



            % Outputs
            err = NaN;

            % Local variables
            srname = "ERRH";
            n = NaN;
            npt = NaN;
            A = NaN(size(xpt, 2));
            e = NaN(3);
            maxabs = NaN;
            Omega = NaN(size(xpt, 2));
            U = NaN(size(xpt, 2));
            V = NaN(size(xpt, 1), size(xpt, 2));
            r = NaN(size(xpt, 2), 1);
            s = NaN(size(xpt, 1), 1);
            t = NaN(size(xpt, 2), 1);

            % Sizes
            n = size(xpt, 1);
            npt = size(xpt, 2);

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(n >= 1, "N >= 1", srname);
                debug_obj.assert(npt >= n + 2, "NPT >= N + 2", srname);
                debug_obj.assert(idz >= 1 && idz <= npt - n, "1 <= IDZ <= NPT-N", srname);
                debug_obj.assert(size(bmat, 1) == n && size(bmat, 2) == npt + n, "SIZE(BMAT)==[N, NPT+N]", srname);
                debug_obj.assert(linalg_obj.issymmetric(bmat(:, npt + 1:npt + n)), "BMAT(:, NPT+1:NPT+N) is symmetric", srname);
                debug_obj.assert(size(zmat, 1) == npt && size(zmat, 2) == npt - n - 1, "SIZE(ZMAT) == [NPT, NPT-N-1]", srname);
                debug_obj.assert(all(infnan_obj.is_finite(xpt), 'all'), "XPT is finite", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            A(:, :) = consts_obj.HALF * fortran.power(linalg_obj.matprod22(xpt', xpt), 2);
            Omega(:, :) = -linalg_obj.matprod22(zmat(:, 1:idz - 1), zmat(:, 1:idz - 1)') + linalg_obj.matprod22(zmat(:, idz:npt - n - 1), zmat(:, idz:npt - n - 1)');
            maxabs = max([consts_obj.ONE, max(abs(A), [], 'all'), max(abs(Omega), [], 'all'), max(abs(bmat), [], 'all')], [], 'all');
            U(:, :) = linalg_obj.eye1(npt) - linalg_obj.matprod22(A, Omega) - linalg_obj.matprod22(xpt', bmat(:, 1:npt));
            V(:, :) = -linalg_obj.matprod22(bmat(:, 1:npt), A) - linalg_obj.matprod22(bmat(:, npt + 1:npt + n), xpt);
            r(:) = sum(U, 1) ./ double(npt);
            s(:) = sum(V, 2) ./ double(npt);
            t(:) = -linalg_obj.matprod21(A, r) - linalg_obj.matprod12(s, xpt);
            e(1, 1) = max(max(U, [], 1) - min(U, [], 1), [], 'all');
            e(1, 2) = max(t, [], 'all') - min(t, [], 'all');
            e(1, 3) = max(max(V, [], 2) - min(V, [], 2), [], 'all');
            e(2, 1) = max(abs(sum(Omega, 1)), [], 'all');
            e(2, 2) = abs(sum(r, 'all') - consts_obj.ONE);
            e(2, 3) = max(abs(sum(bmat(:, 1:npt), 2)), [], 'all');
            e(3, 1) = max(abs(linalg_obj.matprod22(xpt, Omega)), [], 'all');
            e(3, 2) = max(abs(linalg_obj.matprod21(xpt, r)), [], 'all');
            e(3, 3) = max(abs(linalg_obj.matprod22(xpt, bmat(:, 1:npt)') - linalg_obj.eye1(n)), [], 'all');
            err = max(e, [], 'all') / (maxabs * double(n + npt));

            %====================%
            %  Calculation ends  %
            %====================%

        end
        function [idz, bmat, zmat, info] = updateh(obj, knew, kref, d, xpt, idz, bmat, zmat, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine updates arrays [BMAT, ZMAT, IDZ], in order to replace the interpolation point
            % XPT(:, KNEW) by XNEW = XPT(:, KREF) + D, where KREF usually equals KOPT in practice. See Section 4
            % of the NEWUOA paper. [BMAT, ZMAT, IDZ] describes the matrix H in the NEWUOA paper (eq. 3.12),
            % which is the inverse of the coefficient matrix of the KKT system for the least-Frobenius norm
            % interpolation problem: ZMAT holds a factorization of the leading NPT*NPT submatrix OMEGA of H, the
            % factorization being OMEGA = ZMAT*Diag(S)*ZMAT^T with S(1:IDZ-1)= -1 and S(IDZ : NPT-N-1) = +1;
            % BMAT holds the last N ROWs of H except for the (NPT+1)th column. Note that the (NPT + 1)th row and
            % (NPT + 1)th column of H are not stored as they are unnecessary for the calculation. The matrix
            % H is also formulated in (2.7) of the BOBYQA paper. Thanks to the RESCUE method (see Section 5 of
            % the BOBYQA paper), BOBYQA does not have IDZ (equivalent to IDZ = 1).
            %
            % N.B.:
            % 1. What is H? As mentioned above, it is the inverse of the coefficient matrix of the KKT system
            % for the lest-Frobenius norm interpolation problem. Moreover, we should note that the K-th column
            % of H contain the coefficients of the K-th Lagrange function LFUNC_K for this interpolation problem,
            % where 1 <= K <= NPT. More specifically, the first NPT entries of H(:, K) provide the parameters
            % for the Hessian of LFUNC_K so that nabla^2 LFUNC_K = sum_{I=1}^NPT H(I, K) XPT(:, I)*XPT(:, I)^T;
            % the last N entries of H(:, K) constitute precisely the gradient of LFUNC_K at the base point XBASE.
            % Recalling that H is represented by OMEGA and BMAT in the block form elaborated above, we can see
            % OMEGA(:, K) contains the leading NPT entries of H(:, K), while BMAT(:, K) contains the last N.
            % 2. Powell's code normally invokes this subroutine with KREF set to KOPT, which is the index of the
            % current best interpolation point (also the current center of the trust region). In theory, the
            % update should however be independent of KREF. The most natural version (not necessarily the best
            % one in practice) of UPDATEH should work based on [KNEW, XNEW - XPT(:,KNEW)] rather than
            % [KNEW, KREF, D]. UPDATEH needs KREF only for calculating VLAG and BETA, where XPT(:, KREF) is used
            % as a reference point that can be any column of XPT in precise arithmetic. Using XPT(:, KNEW) as
            % the reference point, VLAG and BETA can be calculated by
            % %VLAG = CALVLAG(KNEW, BMAT, XNEW - XPT(:, KNEW), XPT, ZMAT, IDZ)
            % %BETA = CALBETA(KNEW, BMAT, XNEW - XPT(:, KNEW), XPT, ZMAT, IDZ)
            % Theoretically (but not numerically), they should return the same VLAG and BETA as the calls below.
            % However, as observed on 20220412, such an implementation can lead to significant errors in H!
            %--------------------------------------------------------------------------------------------------%

            % Common modules
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();
            infos_obj = infos_mod();
            linalg_obj = linalg_mod(); %, r2update
            string_obj = string_mod();

            % Inputs


            % D(N)
            % XPT(N, NPT)

            % In-outputs

            % BMAT(N, NPT + N)
            % ZMAT(NPT, NPT - N - 1)

            % Outputs


            % Local variables
            srname = "UPDATEH";
            j = NaN;
            ja = NaN;
            jb = NaN;
            jl = NaN;
            n = NaN;
            npt = NaN;
            alpha = NaN;
            beta = NaN;
            denom = NaN;
            grot = NaN(2);
            hcol = NaN(size(bmat, 2), 1);
            scala = NaN;
            scalb = NaN;
            sqrtdn = NaN;
            tau = NaN;
            temp = NaN;
            tempa = NaN;
            tempb = NaN;
            v1 = NaN(size(bmat, 1), 1);
            v2 = NaN(size(bmat, 1), 1);
            vlag = NaN(size(bmat, 2), 1);

            % Debugging variables
            %real(RP) :: beta_test
            %real(RP) :: tol
            %real(RP), allocatable :: vlag_test(:)
            %real(RP), allocatable :: xpt_test(:, :)

            % Sizes
            n = size(xpt, 1);
            npt = size(xpt, 2);

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(n >= 1 && npt >= n + 2, "N >= 1, NPT >= N + 2", srname);
                debug_obj.assert(knew >= 0 && knew <= npt, "0 <= KNEW <= NPT", srname);
                debug_obj.assert(kref >= 1 && kref <= npt, "1 <= KREF <= NPT", srname);
                debug_obj.assert(idz >= 1 && idz <= size(zmat, 2) + 1, "1 <= IDZ <= SIZE(ZMAT, 2) + 1", srname);
                debug_obj.assert(numel(d) == n && all(infnan_obj.is_finite(d), 'all'), "SIZE(D) == N, D is finite", srname);
                debug_obj.assert(size(bmat, 1) == n && size(bmat, 2) == npt + n, "SIZE(BMAT)==[N, NPT+N]", srname);
                debug_obj.assert(linalg_obj.issymmetric(bmat(:, npt + 1:npt + n)), "BMAT(:, NPT+1:NPT+N) is symmetric", srname);
                debug_obj.assert(size(zmat, 1) == npt && size(zmat, 2) == npt - n - 1, "SIZE(ZMAT) == [NPT, NPT - N - 1]", srname);

                for j = 1:npt
                    hcol(1:npt) = obj.omega_col(idz, zmat, j);
                    hcol(npt + 1:npt + n) = bmat(:, j);
                    debug_obj.assert(floor(-log10(eps(class(0.0)))) < floor(-log10(eps(class(0.0)))) || sum(abs(hcol), 'all') > 0, "Column " + string_obj.int2str(j) + " of H is nonzero", srname);
                end

                debug_obj.assert(all(infnan_obj.is_finite(xpt), 'all'), "XPT is finite", srname);

                % Theoretically, CALVLAG and CALBETA should be independent of the reference point XPT(:, KREF).
                % So we test the following. By the implementation of CALVLAG and CALBETA, we are indeed testing
                % H*[w(X_KREF) - w(X_KNEW)] = e_KREF - e_KNEW. Thus H = W^{-1} is also tested to some extend.
                % However, this is expensive to check.
                %if (knew >= 1) then
                %    tol = 1.0E-2_RP  ! W and H are quite ill-conditioned, so we do not test a high precision.
                %    call safealloc(vlag_test, npt + n)
                %    vlag_test = calvlag(knew, bmat, d + (xpt(:, kref) - xpt(:, knew)), xpt, zmat, idz)
                %    call wassert(all(abs(vlag_test - calvlag(kref, bmat, d, xpt, zmat, idz)) <= &
                %        & tol * maxval([ONE, abs(vlag_test)])) .or. precision(0.0_RP) < precision(0.0D0), 'VLAG_TEST == VLAG', srname)
                %    deallocate (vlag_test)
                %    beta_test = calbeta(knew, bmat, d + (xpt(:, kref) - xpt(:, knew)), xpt, zmat, idz)
                %    call wassert(abs(beta_test - calbeta(kref, bmat, d, xpt, zmat, idz)) <= &
                %        & tol * max(ONE, abs(beta_test)) .or. precision(0.0_RP) < precision(0.0D0), 'BETA_TEST == BETA', srname)
                %end if

                % The following is too expensive to check.
                %call wassert(errh(idz, bmat, zmat, xpt) <= tol .or. precision(0.0_RP) < precision(0.0D0), &
                %    & 'H = W^{-1} in (3.12) of the NEWUOA paper', srname)

            end

            %====================%
            % Calculation starts %
            %====================%

            ipObj = inputParser();
            addParameter(ipObj, 'info', NaN);
            parse(ipObj, varargin{:});
            info = ipObj.Results.info;
            if nargout >= 4
                info = infos_obj.INFO_DFT;
            end

            % We must not do anything if KNEW is 0. This can only happen sometimes after a trust-region step.
            if knew <= 0
                % KNEW < 0 is impossible if the input is correct.
                return
            end

            % Set the first NPT components of HCOL to the leading elements of the KNEW-th column of H. Powell's
            % code does this after ZMAT is rotated blow, which saves flops but also introduces rounding errors.
            hcol(1:npt) = obj.omega_col(idz, zmat, knew);
            hcol(npt + 1:npt + n) = bmat(:, knew);

            % Calculate VLAG and BETA according to D.
            % VLAG contains the components of the vector H*w of the updating formula (4.11) in the NEWUOA paper,
            % and BETA holds the value of the parameter that has this name.
            % N.B.: Powell's original comments mention that VLAG is "the vector THETA*WCHECK + e_b of the
            % updating formula (6.11)", which does not match the published version of the NEWUOA paper.
            vlag(:) = obj.calvlag_lfqint(kref, bmat, d, xpt, zmat, 'idz', idz);
            beta = obj.calbeta(kref, bmat, d, xpt, zmat, 'idz', idz); % Nonnegative in precise arithmetic.

            % Calculate the parameters of the updating formula (4.18)--(4.20) in the NEWUOA paper.
            alpha = hcol(knew); % Nonnegative in precise arithmetic.
            tau = vlag(knew); % Nonzero due to the definition of KNEW.
            denom = alpha * beta + tau ^ 2; % Positive in precise arithmetic.

            % After the following line, VLAG = H*w - e_KNEW in the NEWUOA paper (where t = KNEW).
            vlag(knew) = vlag(knew) - consts_obj.ONE;

            % Quite rarely, due to rounding errors, VLAG or BETA may not be finite, and ABS(DENOM) may not be
            % positive. In such cases, [BMAT, ZMAT] would be destroyed by the update, and hence we would rather
            % not update them at all. Or should we simply terminate the algorithm?
            if ~(infnan_obj.is_finite(sum(abs(hcol), 'all') + sum(abs(vlag), 'all') + abs(beta)) && abs(denom) > 0)
                if nargout >= 4
                    info = infos_obj.DAMAGING_ROUNDING;
                end
                return
            end

            % Update the matrix BMAT. It implements the last N rows of (4.11) in the NEWUOA paper.
            v1(:) = (alpha * vlag(npt + 1:npt + n) - tau * hcol(npt + 1:npt + n)) ./ denom;
            v2(:) = (-beta * hcol(npt + 1:npt + n) - tau * vlag(npt + 1:npt + n)) ./ denom;
            bmat(:, :) = bmat + linalg_obj.outprod(v1, vlag) + linalg_obj.outprod(v2, hcol); %call r2update(bmat, ONE, v1, vlag, ONE, v2, hcol)
            % N.B.: The use of OUTPROD is expensive memory-wise, but it is not our concern in this implementation.
            % Numerically, the update above does not guarantee BMAT(:, NPT+1 : NPT+N) to be symmetric.
            A_slice = linalg_obj.symmetrize(bmat(:, npt + 1:npt + n)); bmat(:, npt + 1:npt + n) = A_slice;

            % Apply Givens rotations to put zeros in the KNEW-th row of ZMAT and set JL. After this,
            % ZMAT(KNEW, :) contains at most two nonzero entries ZMAT(KNEW, 1) and ZMAT(KNEW, JL), one
            % corresponding to all the columns of ZMAT that has a coefficient -1 in the factorization of
            % OMEGA (if any), and the other corresponding to all the columns with +1. In specific,
            % 1. If IDZ = 1 (all coefficients are +1 for the columns of ZMAT in the factorization of OMEGA ) or
            % NPT - N (all the coefficients are -1), then JL = 1, and ZMAT(KNEW, 1) is L2-norm of ZMAT(KNEW, :);
            % 2. If 2 <= IDZ <= NPT - N -1, then JL = IDZ, and ZMAT(KNEW, 1) is L2-norm of ZMAT(KNEW, 1 : IDZ-1),
            % while ZMAT(KNEW, JL) is L2 norm of ZMAT(KNEW, IDZ : NPT-N-1).
            % See (4.15)--(4.17) of the NEWUOA paper and the elaboration around them.
            jl = 1; % In the loop below, if 2 <= J < IDZ, then JL = 1; if IDZ < J <= NPT-N-1, then JL = IDZ.
            for j = 2:npt - n - 1
                if j == idz
                    jl = idz; % Do nothing but changing JL from 1 to IDZ. It occurs at most once along the loop.
                    continue
                end

                % Powell's condition in NEWUOA/LINCOA for the IF ... THEN below: IF (ZMAT(KNEW, J) /= 0) THEN
                % A possible alternative: IF (ABS(ZMAT(KNEW, J)) > 1.0E-20 * ABS(ZMAT(KNEW, JL))) THEN
                if abs(zmat(knew, j)) > 1.0e-20 * max(abs(zmat), [], 'all')
                    % Threshold comes from Powell's BOBYQA
                    % Multiply a Givens rotation to ZMAT from the right so that ZMAT(KNEW, [JL,J]) becomes [*,0].
                    grot(:, :) = linalg_obj.planerot(zmat(knew, [jl, j])); %%MATLAB: grot = planerot(zmat(knew, [jl, j])')
                    zmat(:, [jl, j]) = linalg_obj.matprod22(zmat(:, [jl, j]), grot');
                end
                zmat(knew, j) = consts_obj.ZERO;
            end

            sqrtdn = fortran.sqrt(abs(denom));

            if jl == 1
                % Complete the updating of ZMAT when there is only 1 nonzero in ZMAT(KNEW, :) after the rotation.
                % This is the normal case, as IDZ = 1 in precise arithmetic; it also covers the rare case that
                % IDZ = NPT-N, meaning that OMEGA = -ZMAT*ZMAT^T. See (4.18) of the NEWUOA paper for details.
                % Note that (4.18) updates Z_{NPT-N-1}, but the code here updates ZMAT(:, 1). Correspondingly,
                % we implicitly update S_1 to SIGN(DENOM)*S_1 according to (4.18). If IDZ = NPT-N before the
                % update, then IDZ is reduced by 1, and we need to switch ZMAT(:, 1) and ZMAT(:, IDZ) to maintain
                % that S_J = -1 iff 1 <= J < IDZ, which is done after the END IF together with another case.

                %----------------------------------------------------------------------------------------------%
                % Up to now, TEMPA = ZMAT(KNEW, 1) if IDZ = 1 and TEMPA = -ZMAT(KNEW, 1) if IDZ >= 2. However,
                % according to (4.18) of the NEWUOA paper, TEMPB should always be ZMAT(KNEW, 1)/SQRTDN
                % regardless of IDZ. Therefore, the following definition of TEMPB is inconsistent with (4.18).
                % This is probably a BUG. See also Lemma 4 and (5.13) of Powell's paper "On updating the inverse
                % of a KKT matrix". However, the inconsistency is hardly observable in practice, because JL = 1
                % implies IDZ = 1 in precise arithmetic.
                %--------------------------------------------%
                % %tempb = tempa/sqrtdn
                % %tempa = tau/sqrtdn
                %--------------------------------------------%
                % Here is the corrected version (only TEMPB is changed).
                tempa = tau / sqrtdn;
                tempb = zmat(knew, 1) / sqrtdn;
                %----------------------------------------------------------------------------------------------%

                % The following line updates ZMAT(:, 1) according to (4.18) of the NEWUOA paper.
                zmat(:, 1) = tempa * zmat(:, 1) - tempb * vlag(1:npt);

                %----------------------------------------------------------------------------------------------%
                % Zaikun 20220411: The update of IDZ is decoupled from the update of ZMAT, located after END IF.
                %----------------------------------------------------------------------------------------------%
                % The following six lines from Powell's NEWUOA code are obviously problematic --- SQRTDN is
                % always nonnegative. According to (4.18) of the NEWUOA paper, "SQRTDN < 0" and "SQRTDN >= 0"
                % below should be both revised to "DENOM < 0". See also the corresponding part of the LINCOA
                % code. Note that the NEWUOA paper uses SIGMA to denote DENOM. Check also Lemma 4 and (5.13) of
                % Powell's paper "On updating the inverse of a KKT matrix". Note that the BOBYQA code does not
                % have this part, as it does not have IDZ at all.
                % %if (idz == 1 .and. sqrtdn < 0) then
                % %    idz = 2
                % %end if
                % %if (idz >= 2 .and. sqrtdn >= 0) then
                % %    reduce_idz = .true.
                % %end if
                % This is the corrected version, copied from LINCOA.
                % %if (denom < 0) then
                % %    if (idz == 1) then
                % %        idz = 2
                % %    else
                % %        reduce_idz = .true.
                % %    end if
                % %end if
                %----------------------------------------------------------------------------------------------%

            else
                % Complete the updating of ZMAT in the alternative case: ZMAT(KNEW, :) has 2 nonzeros. See (4.19)
                % and (4.20) of the NEWUOA paper.
                % First, set JA and JB so that ZMAT(: [JA, JB]) corresponds to [Z_1, Z_2] in (4.19) when BETA>=0,
                % and corresponds to [Z2, Z1] in (4.20) when BETA<0. In this way, the update of ZMAT(:, [JA, JB])
                % follows the same scheme regardless of BETA. Indeed, since S_1 = 1 and S_2 = -1 in (4.19)-(4.20)
                % as elaborated above the equations, ZMAT(:, [1, JL]) always correspond to [Z_2, Z_1].
                if beta >= 0
                    % ZMAT(:, [JA, JB]) corresponds to [Z_1, Z_2] in (4.19)
                    ja = jl;
                    jb = 1;
                else                    % ZMAT(:, [JA, JB]) corresponds to [Z_2, Z_1] in (4.20)
                    ja = 1;
                    jb = jl;
                end
                % Now update ZMAT(:, [ja, jb]) according to (4.19)--(4.20) of the NEWUOA paper.
                temp = zmat(knew, jb) / denom;
                %tempa = temp * beta
                %tempb = temp * tau
                tempa = (beta / denom) * zmat(knew, jb);
                tempb = (tau / denom) * zmat(knew, jb);
                temp = zmat(knew, ja);
                scala = consts_obj.ONE / fortran.sqrt(abs(beta) * temp ^ 2 + tau ^ 2); % 1/SQRT(ZETA) in (4.19)-(4.20) of NEWUOA paper
                scalb = scala * sqrtdn;
                zmat(:, ja) = scala * (tau * zmat(:, ja) - temp * vlag(1:npt));
                zmat(:, jb) = scalb * (zmat(:, jb) - tempa * hcol(1:npt) - tempb * vlag(1:npt));

                %----------------------------------------------------------------------------------------------%
                % Zaikun 20220411: The update of IDZ is decoupled from the update of ZMAT, located after END IF.
                %----------------------------------------------------------------------------------------------%
                % If and only if DENOM < 0, IDZ will be revised according to the sign of BETA.
                % See (4.19)--(4.20) of the NEWUOA paper.
                % %if (denom < 0) then
                % %    if (beta < 0) then
                % %        idz = idz + 1_IK
                % %    else
                % %        reduce_idz = .true.
                % %    end if
                % %end if
                %----------------------------------------------------------------------------------------------%
            end

            %--------------------------------------------------------------------------------------------------%
            % Zaikun 20220411: The update of IDZ is decoupled from the update of ZMAT, located right below.
            %--------------------------------------------------------------------------------------------------%
            % IDZ is reduced in the following case. Then exchange ZMAT(:, 1) and ZMAT(:, IDZ).
            % %if (reduce_idz) then
            % %    idz = idz - 1_IK
            % %    if (idz > 1) then
            % %        zmat(:, [1_IK, idz]) = zmat(:, [idz, 1_IK])
            % %    end if
            % %end if
            %--------------------------------------------------------------------------------------------------%

            % According to (4.18) and (4.19)--(4.20) of the NEWUOA paper, the coefficients {S_J} need update iff
            % DENOM < 0, in which case one of the S_J will flip the sign when multiplied by SIGN(DENOM), leading
            % to an increase of IDZ (if S_J flipped from 1 to -1) or a decrease (if S_J flipped from -1 to 1).
            if denom < 0
                if idz == 1 || (idz < npt - n && beta < 0)
                    % (4.18), (4.20) of the NEWUOA paper
                    idz = idz + 1;
                elseif idz == npt - n || (idz > 1 && beta >= 0)
                    % (4.18), (4.19) of the NEWUOA paper
                    idz = idz - 1;
                    % Exchange ZMAT(:, 1) and ZMAT(:. IDZ) if IDZ > 1. Why? No matter whether the update is
                    % given by (4.18) (IDZ = NPT-N) or (4.19) (1 < IDZ < NPT-N and BETA >= 0), we have S_1 = +1
                    % and S_{IDZ} = -1 at this moment (unless IDZ = 1). Thus we need to exchange ZMAT(:, 1) with
                    % ZMAT(:, IDZ) and implicitly S_1 with S_IDZ to maintain that S_J = -1 iff 1 <= J < IDZ.
                    % Note that, in the case of (4.18), ZMAT(:, 1) (and implicitly S_1) rather than
                    % ZMAT(:, NPT-N-1) was updated by the code above.
                    if idz > 1
                        zmat(:, [1, idz]) = zmat(:, [idz, 1]);
                    end
                end
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(idz >= 1 && idz <= size(zmat, 2) + 1, "1 <= IDZ <= SIZE(ZMAT, 2) + 1", srname);
                debug_obj.assert(size(bmat, 1) == n && size(bmat, 2) == npt + n, "SIZE(BMAT)==[N, NPT+N]", srname);
                debug_obj.assert(linalg_obj.issymmetric(bmat(:, npt + 1:npt + n)), "BMAT(:, NPT+1:NPT+N) is symmetric", srname);
                debug_obj.assert(size(zmat, 1) == npt && size(zmat, 2) == npt - n - 1, "SIZE(ZMAT) == [NPT, NPT-N-1]", srname);

                for j = 1:npt
                    hcol(1:npt) = obj.omega_col(idz, zmat, j);
                    hcol(npt + 1:npt + n) = bmat(:, j);
                    debug_obj.assert(floor(-log10(eps(class(0.0)))) < floor(-log10(eps(class(0.0)))) || sum(abs(hcol), 'all') > 0, "Column " + string_obj.int2str(j) + " of H is nonzero", srname);
                end

                % The following is too expensive to check.
                %call safealloc(xpt_test, n, npt)
                %xpt_test = xpt
                %xpt_test(:, knew) = xpt(:, kref) + d
                %call wassert(errh(idz, bmat, zmat, xpt_test) <= tol .or. precision(0.0_RP) < precision(0.0D0), &
                %    & 'H = W^{-1} in (3.12) of the NEWUOA paper', srname)
                %deallocate (xpt_test)

            end
        end
        %--------------------------------------------------------------------------------------------------%
        % CALVLAG, CALBETA, and CALDEN are subroutine that calculate VLAG, BETA, and DEN for a given step D.
        % VLAG(K), BETA, and DEN(K) are critical for the updating procedure of H when the interpolation set
        % replaces XPT(:, K) with D. The updating formula of H is detailed in (4.11) of the NEWUOA paper,
        % where the point being replaced is XPT(:, t). See (4.12) for the definition of BETA; VLAG is indeed
        % H*w without the (NPT+1)the entry; DEN(t) is SIGMA in (4.12). (4.25)--(4.26) formulate the actual
        % calculating scheme of VLAG and BETA.
        %
        % Zaikun 20250806: In theory, ALPHA and BETA are both nonnegative (see Lemma 2 of Powell 2002, Least
        % Frobenius norm updating of quadratic models that satisfy interpolation conditions), and DEN is
        % positive. Numerically, however, they can be negative due to rounding errors. Powell's code handles
        % such cases carefully. See the UPDATEH subroutine and Section 4 of the NEWUOA paper, especially the
        % discussions around (4.21)--(4.22).
        %
        % In languages like MATLAB/Python/Julia/R, CALVLAG and CALBETA should be implemented into one single
        % function, as they share most of the calculation. We separate them in Fortran (at the expense of
        % repeating some calculation) because Fortran functions can only have one output.
        %
        % Explanation on the matrix H in (3.12) and w(X) in (6.3) of the NEWUOA paper and WCHECK in the code:
        % 0. As defined in (6.3) of the paper w(X)(K) = 0.5*[(X-X_0)^T*(X_K-X_0)]^2 for K = 1, ..., NPT,
        % w(X)(NPT+1) = 1, and w(X)(NPT+2:NPT+N+1) = X - X_0. As in (3.12) of the paper, w(X_K) is the K-th
        % column in the coefficient matrix W of the KKT system for the interpolation problem
        % Minimize ||nabla^2 Q||_F s.t. Q(X_K) = Y_K, K = 1, ..., NPT.
        % This is why w(X) is ubiquitous in the code.
        % 1. By (3.9) of the paper, the solution to the above interpolation problem is Q(X) = Y^T*H*w(X),
        % with H = W^{-1}, and Y being the vector [Y_1; ...; Y_NPT; 0; ...; 0] with N trailing zeros. In
        % particular, the K-th Lagrange function of this interpolation problem is e_K^T*H*w(X), namely the
        % K-th entry of the vector H*w(X). This is why H*w(X) appears as the vector VLAG in the code.
        % As a consequence, SUM(VLAG(1:NPT)) = 1 in theory.
        % 2. As above, H can provide us interpolants and the Lagrange functions. Thus the code maintains H.
        % Indeed, the K-th column of H contain the coefficients of the K-th Lagrange function LFUNC_K,
        % where 1 <= K <= NPT. More specifically, the first NPT entries of H(:, K) provide the parameters
        % for the Hessian of LFUNC_K so that nabla^2 LFUNC_K = sum_{I=1}^NPT H(I, K) XPT(:, I)*XPT(:, I)^T;
        % the last N entries of H(:, K) constitute precisely the gradient of LFUNC_K at the base point X_0.
        % Recalling that H (except for the (NPT+1)the row and column) is represented by OMEGA and BMAT in
        % the block form [OMEGA, BMAT(:, 1:NPT)'; BMAT(:, 1:NPT), BMAT(:, NPT+1:NPT+N)], we can see
        % OMEGA(:, K) contains the leading NPT entries of H(:, K), while BMAT(:, K) contains the last N.
        % Hence, if X corresponds to XOPT + D, then for K /= KOPT, the K-th entry of VLAG = H*w(X) equals
        % LFUNC_K(X_0 + XOPT + D) - LFUNC_K(X_0 + XOPT) = QUADINC_DX(D, X, XPT, BMAT(:, K), OMEGA(:, K)),
        % because LFUNC_K(X_0 + XOPT) = 0; for K = KOPT, it equals
        % LFUNC_K(X_0 + XOPT + D)-LFUNC_K(X_0 + XOPT)+1 = QUADINC_DX(D, X, XPT, BMAT(:, K), OMEGA(:, K)) + 1,
        % as LFUNC_K(X_0 + XOPT) = 1 in this case.
        % 3. Since the matrix H is W^{-1} as defined in (3.12) of the paper, we have H*w(X_K) = e_K for
        % any K in {1, ..., NPT}.
        % 4. When the interpolation set is updated by replacing X_K with X, W is correspondingly updated by
        % changing the K-th column from w(X_K) to w(X). This is why the update of H = W^{-1} must involve
        % H*[w(X) - w(X_K)] = H*w(X) - e_K.
        % 5. As explained above, the vector H*w(X) is essential to the algorithm. The quantity w(X)^T*H*w(X)
        % is also needed in the update of H (particularly by BETA). However, they can be tricky to calculate,
        % because much cancellation can happen when X_0 is far away from the interpolation set, as explained
        % in (7.8)--(7.10) of the paper and the discussions around. To overcome the difficulty, we take
        % an integer KREF in {1, ..., NPT}, use XPT(:, KREF) as a reference point, and note that
        % H*w(X) = H*[w(X) - w(X_KREF)] + H*w(X_KREF) = H*[w(X) - w(X_KREF)] + e_KREF, and
        % w(x)^T*H*w(X) = [w(X) - w(X_KREF)]^T*H*[w(X) - w(X_KREF)] + 2*w(X)(KREF) - w(X_KREF)(KREF),
        % The dependence of w(X)-w(X_KREF) on X_0 is weaker, which reduces (but does not resolve) the
        % difficulty. In theory, these formulas are invariant with respect to KREF. In the code, this means
        % CALVLAG(KREF, BMAT, X - XPT(:, KREF), XPT, ZMAT, IDZ) and
        % CALBETA(KREF, BMAT, X - XPT(:, KREF), XPT, ZMAT, IDZ)
        % are invariant with respect to KREF. Powell's code normally uses KREF = KOPT.
        % 6. Since the (NPT+1)-th entry of w(X) - w(X_KREF) is 0, the above formulas do not require the
        % (NPT+1)-th column of H, which is not stored in the code.
        % 7. In the code, WCHECK contains the first NPT entries of w-v for the vectors w and v in (4.10) and
        % (4.24) of the NEWUOA paper, with w = w(X) and v = w(X_KREF) (KREF = KOPT in the paper); it is
        % also hat{w} in (6.5) of M. J. D. Powell, Least Frobenius norm updating of quadratic models that
        % satisfy interpolation conditions. Math. Program., 100:183--215, 2004 (KREF = b in the paper).
        % 8. Assume that the ||D|| ~ DELTA, ||XPT|| ~ ||XREF||, and DELTA < ||XREF||. Then WCHECK is of the
        % order DELTA*||XREF||^3, which can be huge at the beginning of the algorithm and quickly become tiny.
        %--------------------------------------------------------------------------------------------------%
        function vlag = calvlag_lfqint(obj, kref, bmat, d, xpt, zmat, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This function calculates VLAG = H*w for a given step D with respect to XREF = XPT(:, KREF). This
            % subroutine is usually invoked with KREF = KOPT, which correspond to the current best interpolation
            % point as well as the center of the trust region. See (4.25) of the NEWUOA paper.
            %--------------------------------------------------------------------------------------------------%

            % Common modules
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();
            linalg_obj = linalg_mod();


            % Inputs

            % BMAT(N, NPT + N)
            % D(N)
            % XPT(N, NPT)
            % ZMAT(NPT, NPT - N - 1)
            % Absent in BOBYQA, being equivalent to IDZ = 1

            % Outputs
            vlag = NaN(size(xpt, 1) + size(xpt, 2), 1); % VLAG(NPT + N)

            % Local variables
            srname = "CALVLAG";
            idz_loc = NaN;
            n = NaN;
            npt = NaN;
            tol = NaN; % For debugging only
            wcheck = NaN(size(zmat, 1), 1);
            xref = NaN(size(xpt, 1), 1);

            % Sizes
            n = size(xpt, 1);
            npt = size(xpt, 2);

            % Read IDZ, which is not present in BOBYQA, being equivalent to IDZ = 1.
            idz_loc = 1;
            ipObj = inputParser();
            addParameter(ipObj, 'idz', NaN);
            parse(ipObj, varargin{:});
            idz = ipObj.Results.idz;
            if ~ismember('idz', ipObj.UsingDefaults)
                idz_loc = idz;
            end

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(n >= 1 && npt >= n + 2, "N >= 1, NPT >= N + 2", srname);
                debug_obj.assert(idz_loc >= 1 && idz_loc <= size(zmat, 2) + 1, "1 <= ID <= SIZE(ZMAT, 2) + 1", srname);
                debug_obj.assert(kref >= 1 && kref <= npt, "1 <= KREF <= NPT", srname);
                debug_obj.assert(size(bmat, 1) == n && size(bmat, 2) == npt + n, "SIZE(BMAT)==[N, NPT+N]", srname);
                debug_obj.assert(linalg_obj.issymmetric(bmat(:, npt + 1:npt + n)), "BMAT(:, NPT+1:NPT+N) is symmetric", srname);
                debug_obj.assert(size(zmat, 1) == npt && size(zmat, 2) == npt - n - 1, "SIZE(ZMAT) == [NPT, NPT - N - 1]", srname);
                debug_obj.assert(numel(d) == n && all(infnan_obj.is_finite(d), 'all'), "SIZE(D) == N, D is finite", srname);
                debug_obj.assert(all(infnan_obj.is_finite(xpt), 'all'), "XPT is finite", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            xref(:) = xpt(:, kref); % Read XREF.

            % Set WCHECK to the first NPT entries of (w-v) for w and v in (4.10) and (4.24) of the NEWUOA paper.
            wcheck(:) = linalg_obj.matprod12(d, xpt);
            wcheck(:) = wcheck .* (consts_obj.HALF * wcheck + linalg_obj.matprod12(xref, xpt));

            % The following two lines set VLAG to H*(w-v).
            vlag(1:npt) = obj.omega_mul(idz_loc, zmat, wcheck) + linalg_obj.matprod12(d, bmat(:, 1:npt));
            vlag(npt + 1:npt + n) = linalg_obj.matprod21(bmat, [reshape(wcheck, [], 1); reshape(d, [], 1)]);
            % The following line is equivalent to the above one, but handles WCHECK and D separately.
            % %vlag(npt + 1:npt + n) = matprod(bmat(:, 1:npt), wcheck) + matprod(bmat(:, npt + 1:npt + n), d)

            % The following line sets VLAG(KREF) to the correct value.
            vlag(kref) = vlag(kref) + consts_obj.ONE;

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(vlag) == npt + n, "SIZE(VLAG) == NPT + N", srname);
                tol = max(consts_obj.TEN ^ max(-8, -consts_obj.MAXPOW10), min(0.1, consts_obj.TEN ^ min(12, consts_obj.MAXPOW10) * consts_obj.EPS * double(npt + n)));
                debug_obj.wassert(abs(sum(vlag(1:npt), 'all') - consts_obj.ONE) / double(npt) <= tol || floor(-log10(eps(class(0.0)))) < floor(-log10(eps(class(0.0)))), "SUM(VLAG(1:NPT)) == 1", srname);
            end

        end
        function beta = calbeta(obj, kref, bmat, d, xpt, zmat, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This function calculates BETA for a given step D with respect to XREF = XPT(:, KREF). This
            % subroutine is usually invoked with KREF = KOPT, which correspond to the current best interpolation
            % point as well as the center of the trust region. See (4.12) and (4.26) of the NEWUOA paper.
            %--------------------------------------------------------------------------------------------------%

            % Common modules
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();
            linalg_obj = linalg_mod();


            % Inputs

            % BMAT(N, NPT + N)
            % D(N)
            % XPT(N, NPT)
            % ZMAT(NPT, NPT - N - 1)
            % Absent in BOBYQA, being equivalent to IDZ = 1

            % Outputs
            beta = NaN;

            % Local variables
            srname = "CALBETA";
            idz_loc = NaN;
            n = NaN;
            npt = NaN;
            dsq = NaN;
            dvlag = NaN;
            dxref = NaN;
            vlag = NaN(size(xpt, 1) + size(xpt, 2), 1);
            wcheck = NaN(size(zmat, 1), 1);
            wmv = NaN(size(xpt, 1) + size(xpt, 2), 1);
            wvlag = NaN;
            xref = NaN(size(xpt, 1), 1);
            xrefsq = NaN;

            % Sizes
            n = size(xpt, 1);
            npt = size(xpt, 2);

            % Read IDZ, which is absent from BOBYQA, being equivalent to IDZ = 1.
            idz_loc = 1;
            ipObj = inputParser();
            addParameter(ipObj, 'idz', NaN);
            parse(ipObj, varargin{:});
            idz = ipObj.Results.idz;
            if ~ismember('idz', ipObj.UsingDefaults)
                idz_loc = idz;
            end

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(n >= 1 && npt >= n + 2, "N >= 1, NPT >= N + 2", srname);
                debug_obj.assert(idz_loc >= 1 && idz_loc <= size(zmat, 2) + 1, "1 <= IDZ <= SIZE(ZMAT, 2) + 1", srname);
                debug_obj.assert(kref >= 1 && kref <= npt, "1 <= KREF <= NPT", srname);
                debug_obj.assert(size(bmat, 1) == n && size(bmat, 2) == npt + n, "SIZE(BMAT)==[N, NPT+N]", srname);
                debug_obj.assert(linalg_obj.issymmetric(bmat(:, npt + 1:npt + n)), "BMAT(:, NPT+1:NPT+N) is symmetric", srname);
                debug_obj.assert(size(zmat, 1) == npt && size(zmat, 2) == npt - n - 1, "SIZE(ZMAT) == [NPT, NPT - N - 1]", srname);
                debug_obj.assert(numel(d) == n && all(infnan_obj.is_finite(d), 'all'), "SIZE(D) == N, D is finite", srname);
                debug_obj.assert(all(infnan_obj.is_finite(xpt), 'all'), "XPT is finite", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            xref(:) = xpt(:, kref); % Read XREF.

            %--------------------------------------------------------------------------------------------------%
            % N.B.: When checking the NEWUOA paper, note that the paper takes KREF = KOPT and XREF = XOPT.
            %--------------------------------------------------------------------------------------------------%

            % Set WCHECK to the first NPT entries of (w-v) for w and v in (4.10) and (4.24) of the NEWUOA paper.
            wcheck(:) = linalg_obj.matprod12(d, xpt);
            wcheck(:) = wcheck .* (consts_obj.HALF * wcheck + linalg_obj.matprod12(xref, xpt));

            % WMV is the vector (w-v) for w and v in (4.10) and (4.24) of the NEWUOA paper.
            wmv(:) = [reshape(wcheck, [], 1); reshape(d, [], 1)];
            % The following two lines set VLAG to H*(w-v).
            vlag(1:npt) = obj.omega_mul(idz_loc, zmat, wcheck) + linalg_obj.matprod12(d, bmat(:, 1:npt));
            vlag(npt + 1:npt + n) = linalg_obj.matprod21(bmat, wmv);
            % The following line is equivalent to the above one, but handles WCHECK and D separately.
            % %VLAG(NPT + 1:NPT + N) = MATPROD(BMAT(:, 1:NPT), WCHECK) + MATPROD(BMAT(:, NPT + 1:NPT + N), D)

            % Set BETA = HALF*||XREF + D||^4 - (W-V)'*H*(W-V) - [XREF'*(X+XREF)]^2 + HALF*||XREF||^4. See
            % equations (4.10), (4.12), (4.24), and (4.26) of the NEWUOA paper.
            dxref = linalg_obj.inprod(d, xref);
            dsq = linalg_obj.inprod(d, d);
            xrefsq = linalg_obj.inprod(xref, xref);
            dvlag = linalg_obj.inprod(d, vlag(npt + 1:npt + n));
            wvlag = linalg_obj.inprod(wcheck, vlag(1:npt));
            beta = dxref ^ 2 + dsq * (xrefsq + dxref + dxref + consts_obj.HALF * dsq) - dvlag - wvlag;
            %---------------------------------------------------------------------------------------------------%
            % The last line is equivalent to either of the following lines, but performs better numerically.
            % %BETA = DXREF**2 + DSQ * (XREFSQ + DXREF + DXREF + HALF * DSQ) - INPROD(VLAG, WMV)  ! not good
            % %BETA = DXREF**2 + DSQ * (XREFSQ + DXREF + DXREF + HALF * DSQ) - WVLAG - DVLAG  ! bad
            %---------------------------------------------------------------------------------------------------%

            % N.B.:
            % 1. Mathematically, the following two quantities are equal:
            % DXREF**2 + DSQ * (XREFSQ + DXREF + DXREF + HALF * DSQ) ,
            % HALF * (INPROD(X, X)**2 + INPROD(XREF, XREF)**2) - INPROD(X, XREF)**2 with X = XREF + D.
            % However, the first (by Powell) is a better numerical scheme. According to the first formulation,
            % this quantity is in the order of ||D||^2*||XREF||^2 if ||XREF|| >> ||D||, which is normally the case.
            % However, each term in the second formulation has an order of ||XREF||^4. Thus much cancellation
            % will occur in the second formulation. In addition, the first formulation contracts the rounding
            % error in (XREFSQ + DXREF + DXREF + HALF * DSQ) by a factor of ||D||^2, which is typically small.
            % 2. We can evaluate INPROD(VLAG, WMV) as INPROD(VLAG(1:NPT), WCHECK) + INPROD(VLAG(NPT+1:NPT+N),D)
            % if it is desirable to handle WCHECK and D separately due to their significantly different magnitudes.

            % The following line sets VLAG(KREF) to the correct value if we intend to output VLAG.
            % %VLAG(KREF) = VLAG(KREF) + ONE

            %====================%
            %  Calculation ends  %
            %====================%

        end
        function den = calden(obj, kref, bmat, d, xpt, zmat, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This function calculates DEN for a given step D with respect to XREF = XPT(:, KREF). DEN is an
            % array of length NPT, and DEN(K) is the value of SIGMA in (4.12) of the NEWUOA paper if XPT(:, K)
            % is replaced with XPT(:, KREF)+D. This value appears as a DENominator in the updating formula of
            % the matrix H as detailed in (4.11) of the NEWUOA paper. This subroutine is usually invoked with
            % KREF = KOPT, which correspond to the current best interpolation point as well as the center of the
            % trust region.
            %--------------------------------------------------------------------------------------------------%

            % Common modules
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            infnan_obj = infnan_mod();
            linalg_obj = linalg_mod();


            % Inputs

            % BMAT(N, NPT + N)
            % D(N)
            % XPT(N, NPT)
            % ZMAT(NPT, NPT - N - 1)
            % Absent in BOBYQA, being equivalent to IDZ = 1

            % Outputs
            den = NaN(size(xpt, 2), 1);

            % Local variables
            srname = "CALDEN";
            idz_loc = NaN;
            n = NaN;
            npt = NaN;
            beta = NaN;
            hdiag = NaN(size(xpt, 2), 1);
            vlag = NaN(size(xpt, 1) + size(xpt, 2), 1);

            % Sizes
            n = size(xpt, 1);
            npt = size(xpt, 2);

            % Read IDZ, which is absent from BOBYQA, being equivalent to IDZ = 1.
            idz_loc = 1;
            ipObj = inputParser();
            addParameter(ipObj, 'idz', NaN);
            parse(ipObj, varargin{:});
            idz = ipObj.Results.idz;
            if ~ismember('idz', ipObj.UsingDefaults)
                idz_loc = idz;
            end

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(n >= 1 && npt >= n + 2, "N >= 1, NPT >= N + 2", srname);
                debug_obj.assert(idz_loc >= 1 && idz_loc <= size(zmat, 2) + 1, "1 <= IDZ <= SIZE(ZMAT, 2) + 1", srname);
                debug_obj.assert(kref >= 1 && kref <= npt, "1 <= KREF <= NPT", srname);
                debug_obj.assert(size(bmat, 1) == n && size(bmat, 2) == npt + n, "SIZE(BMAT)==[N, NPT+N]", srname);
                debug_obj.assert(linalg_obj.issymmetric(bmat(:, npt + 1:npt + n)), "BMAT(:, NPT+1:NPT+N) is symmetric", srname);
                debug_obj.assert(size(zmat, 1) == npt && size(zmat, 2) == npt - n - 1, "SIZE(ZMAT) == [NPT, NPT - N - 1]", srname);
                debug_obj.assert(numel(d) == n && all(infnan_obj.is_finite(d), 'all'), "SIZE(D) == N, D is finite", srname);
                debug_obj.assert(all(infnan_obj.is_finite(xpt), 'all'), "XPT is finite", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            hdiag(:) = -sum(fortran.power(zmat(:, 1:idz_loc - 1), 2), 2) + sum(fortran.power(zmat(:, idz_loc:size(zmat, 2)), 2), 2);
            vlag(:) = obj.calvlag_lfqint(kref, bmat, d, xpt, zmat, 'idz', idz_loc);
            beta = obj.calbeta(kref, bmat, d, xpt, zmat, 'idz', idz_loc);
            den(:) = hdiag * beta + fortran.power(vlag(1:npt), 2);

            %====================%
            %  Calculation ends  %
            %====================%

        end
        function vlag = calvlag_qint(~, pl, d, xref, kref)
            %--------------------------------------------------------------------------------------------------%
            % This function evaluates VLAG = [LFUNC_1(XREF+D), ..., LFUNC_NPT(XREF+D)] for a quadratic
            % interpolation problem, where LFUNC_K is the K-the Lagrange function, and XREF is the KREF-th
            % interpolation node. This subroutine is usually invoked with KREF = KOPT, which correspond to the
            % current best interpolation point as well as the center of the trust region.
            % The coefficients of LFUNC_K are provided by PL(:, K) so that LFUNC_K(Y) = <Y, G> + <HESSIAN*Y, Y>,
            % where G is PL(1:N, K), and HESSIAN is the symmetric matrix whose upper triangular part is stored
            % in PL(N+1:N*(N+3)/2, K) column by column. Note the following:
            % 1. For K /= KREF, LFUNC_K(XREF + D) = LFUNC_K(XREF + D) - LFUNC_K(XREF) as LFUNC_K(XREF) = 0.
            % 2. When K = KREF, LFUNC_K(XREF + D) = LFUNC_K(XREF + D) - LFUNC_K(XREF) + 1 as LFUNC_K(XREF) = 1.
            % Therefore, the function first calculates VLAG(K) = QUADINC_GHV(PL(:, K), D, XREF) for each K, and
            % then increase VLAG(KREF) by 1.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            linalg_obj = linalg_mod();
            % Inputs



            % Outputs
            vlag = NaN(size(pl, 2), 1);
            % Local variables
            srname = "CALVLAG_QINT";
            ih = NaN;
            n = NaN;
            j = NaN;
            s = NaN(numel(xref), 1);
            w = NaN(size(pl, 1), 1);

            % Sizes
            n = fix(numel(xref));

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(d) == n, "SIZE(D) = N", srname);
                debug_obj.assert(size(pl, 2) == (n + 1) * (n + 2) / 2 && size(pl, 1) == size(pl, 2) - 1, "SIZE(PL) = [N*(N+3)/2, (N+1)*(N+2)/2]", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            s(:) = xref + d;

            w(1:n) = d;
            for j = 1:n
                ih = n + (j - 1) * j / 2;
                w(ih + 1:ih + j) = d(1:j) * s(j) + d(j) * xref(1:j);
                w(ih + j) = consts_obj.HALF * w(ih + j);
            end

            vlag(:) = linalg_obj.matprod12(w, pl); % VLAG(K) = QUADINC_GHV(PL(:, K), D, XREF)
            vlag(kref) = vlag(kref) + consts_obj.ONE;

            %====================%
            % Calculation ends   %
            %====================%

        end
        function ij = setij(~, n, npt, varargin)
            %--------------------------------------------------------------------------------------------------%
            % Set IJ to a 2-by-(NPT-2*N-1) integer array so that IJ(:, K) = [P(K + 2*N + 1), Q(K + 2*N + 1)],
            % with P and Q defined in (2.4) of the BOBYQA paper as well as Section 3 of the NEWUOA paper.
            % If NPT <= 2*N + 1, then IJ is empty. Assume that NPT >= 2*N + 2. Then SIZE(IJ) = [2, NPT-2*N-1].
            % IJ contains integers between 1 and N. When NPT = (N+1)*(N+2)/2, the columns of IJ correspond to
            % a permutation of {{I, J} : 1 <= I /= J <= N}; when NPT < (N+1)*(N+2)/2, they correspond to the
            % first NPT - 2*N - 1 elements of such a permutation. The permutation is enumerated first in the
            % ascending order of |I - J| and then in the ascending order of MAX{I, J}. We do not distinguish
            % between {I, J} and {J, I}, which represent the same set. If we want to ensure IJ(1, :) > IJ(2, :),
            % then we can specify SORTING_DIRECTION = 'descend'.
            %
            % This function is used in the initialization of NEWUOA, BOBYQA, and LINCOA.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            linalg_obj = linalg_mod();
            string_obj = string_mod();

            % Inputs



            % Outputs
            ij = NaN(2, max(0, npt - 2 * n - 1));
            % Local variables
            srname = "SETIJ";
            k = NaN;
            ell = NaN(max(0, npt - 2 * n - 1), 1);

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(n >= 1 && npt >= n + 2, "N >= 1, NPT >= N + 2", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            ell(:) = fix(reshape((n:npt - n - 2), [], 1) ./ n); % The ell below (2.4) of the BOBYQA paper.
            ij(1, :) = reshape((n:npt - n - 2), [], 1) - n * ell + 1;
            ij(2, :) = mod(ij(1, :) + ell - 1, n) + 1; % MODULO(K-1, N) + 1 = K-N for K in [N+1, 2N]
            ipObj = inputParser();
            addParameter(ipObj, 'sorting_direction', "");
            parse(ipObj, varargin{:});
            sorting_direction = ipObj.Results.sorting_direction;
            if ~ismember('sorting_direction', ipObj.UsingDefaults)
                ij(:, :) = linalg_obj.sort_i2(ij, 'dim', 1, 'direction', sorting_direction); % SORTING_DIRECTION is 'DESCEND' of 'ASCEND'

            end
            %%MATLAB: (N.B.: Fortran MODULO == MATLAB `mod`, Fortran MOD == MATLAB `rem`)
            %%ell = floor((n : npt-n-2) / n);
            %%ij(1, :) = (n : npt-n-2) - n*ell + 1;
            %%ij(2, :) = mod(ij(1, :) + ell - 1, n) + 1;  % mod(k-1,n) + 1 = k-n for k in [n+1,2n]
            %%if nargin >= 3
            %%    ij = sort(ij, 2, sorting_direction)  % `sorting_direction` is 'descend' of 'ascend'
            %%end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(size(ij, 1) == 2 && size(ij, 2) == max(0, npt - 2 * n - 1), "SIZE(IJ) == [2, NPT - 2*N - 1]", srname);
                debug_obj.assert(all(ij >= 1 & ij <= n, 'all'), "1 <= IJ <= N", srname);
                if ismember('sorting_direction', ipObj.UsingDefaults)
                    debug_obj.assert(all(ij(1, :) ~= ij(2, :), 'all'), "IJ(1, :) /= IJ(2, :)", srname);
                else
                    if string_obj.lower(sorting_direction) == "descend"
                        debug_obj.assert(all(ij(1, :) > ij(2, :), 'all'), "IJ(1, :) > IJ(2, :)", srname);
                    elseif string_obj.lower(sorting_direction) == "ascend"
                        debug_obj.assert(all(ij(1, :) < ij(2, :), 'all'), "IJ(1, :) < IJ(2, :)", srname);
                    end
                end
            end
        end

    end
end
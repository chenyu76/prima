classdef shiftbase_mod
    %--------------------------------------------------------------------------------------------------%
    % This module contains a subroutine that shifts the base point from XBASE to XBASE + XPT. It is
    % used in NEWUOA, BOBYQA, and LINCOA.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the NEWUOA paper.
    %
    % Dedicated to late Professor M. J. D. Powell FRS (1936--2015).
    %
    % Started: July 2020
    %
    % Last Modified: Thu 14 Aug 2025 07:35:47 AM CST
    %--------------------------------------------------------------------------------------------------%

    methods
        function varargout = shiftbase(obj, varargin)
            if numel(varargin) >= 7 && numel(varargin) <= 8 && isvector(varargin{2}) && (~isvector(varargin{3}) && ~isscalar(varargin{3})) && (~isvector(varargin{4}) && ~isscalar(varargin{4}))
                [varargout{1:nargout}] = obj.shiftbase_lfqint(varargin{:});
            else
                [varargout{1:nargout}] = obj.shiftbase_qint(varargin{:});
            end
        end
        function [xbase, xpt, bmat, hq] = shiftbase_lfqint(~, kopt, xbase, xpt, zmat, bmat, pq, hq, varargin)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine shifts the base point from XBASE to XBASE + XOPT and updates BMAT and HQ
            % accordingly. PQ and ZMAT remain the same after the shifting. See Section 7 of the NEWUOA paper.
            % N.B.:
            % 1. In Powell's implementation of NEWUOA, the quadratic model is represented by [GQ, PQ, HQ], where
            % GQ is the gradient of the quadratic model at XBASE. In that case, GQ should be updated by
            % GQ = GQ + HESS_MUL(XOPT, XPT, PQ, HQ) in this subroutine, where PQ is the un-updated version.
            % However, Powell implemented BOBYQA and LINCOA without GQ but with GOPT, which is the gradient at
            % XBASE + XOPT. Note that GOPT remains unchanged when XBASE is shifted. In our implementation,
            % NEWUOA also uses GOPT instead of GQ.
            % 2. [IDZ, ZMAT] provides the factorization of Omega in (3.17) of the NEWUOA paper; in specific,
            % Omega = sum_{i=1}^{NPT-N-1} s_i*ZMAT(:,i)*ZMAT(:,i)^T, s_i = -1 if i < IDZ, and si = 1 if i >= IDZ.
            % In precise arithmetic, IDZ should be always 1; to cope with rounding errors, NEWUOA and LINCOA
            % allow IDZ = -1 (see (4.18)--(4.20) of the NEWUOA paper); in BOBYQA, IDZ is always 1, and the
            % rounding errors are handled by the RESCUE subroutine (Sec. 5 of the BOBYQA paper).
            %--------------------------------------------------------------------------------------------------%

            % Common modules
            consts_obj = prima_mat.common.consts_mod();
            debug_obj = prima_mat.common.debug_mod();
            infnan_obj = prima_mat.common.infnan_mod();
            linalg_obj = prima_mat.common.linalg_mod();


            % Inputs

            % PQ(NPT)
            % ZMAT(NPT, NPT - N - 1)
            % Absent in BOBYQA, being equivalent to IDZ = 1

            % In-outputs
            % BMAT(N, NPT + N)
            % HQ(N, N)
            % XBASE(N)
            % XPT(N, NPT)

            % Local variables
            srname = "SHIFTBASE_LFQINT";


            bymat = NaN(numel(xbase));
            %real(RP) :: htol

            sxpt = NaN(size(xpt, 2), 1);
            v = NaN(numel(xbase), 1);
            vxopt = NaN(numel(xbase));
            xopt = NaN(numel(xbase), 1);

            xptxav = NaN(size(xpt, 1), size(xpt, 2));
            ymat = NaN(size(xpt, 1), size(xpt, 2));
            yzmat = NaN(numel(xbase), size(zmat, 2));
            yzmat_c = NaN(numel(xbase), size(zmat, 2));

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
                debug_obj.assert(numel(xbase) == n && all(infnan_obj.is_finite(xbase), 'all'), "SIZE(XBASE) == N, XBASE is finite", srname);
                debug_obj.assert(all(infnan_obj.is_finite(xpt), 'all'), "XPT is finite", srname);
                debug_obj.assert(kopt >= 1 && kopt <= npt, "1 <= KOPT <= NPT", srname);
                debug_obj.assert(idz_loc >= 1 && idz_loc <= size(zmat, 2) + 1, "1 <= IDZ <= SIZE(ZMAT, 2) + 1", srname);
                debug_obj.assert(size(bmat, 1) == n && size(bmat, 2) == npt + n, "SIZE(BMAT)==[N, NPT+N]", srname);
                debug_obj.assert(linalg_obj.issymmetric(bmat(:, npt + 1:npt + n)), "BMAT(:, NPT+1:NPT+N) is symmetric", srname);
                debug_obj.assert(size(zmat, 1) == npt && size(zmat, 2) == npt - n - 1, "SIZE(ZMAT) == [NPT, NPT - N - 1]", srname);
                debug_obj.assert(numel(pq) == npt, "SIZE(PQ) = NPT", srname);
                debug_obj.assert(size(hq, 1) == n && linalg_obj.issymmetric(hq), "HQ is an NxN symmetric matrix", srname);
                % The following test cannot be passed.
                %htol = max(TEN**max(-10, -MAXPOW10), min(1.0E-1_RP, TEN**min(10, MAXPOW10) * EPS)) ! Tolerance for error in H
                %call assert(errh(idz_loc, bmat, zmat, xpt) <= htol, 'H = W^{-1} in (3.12) of the NEWUOA paper', srname)

            end

            %====================%
            % Calculation starts %
            %====================%

            % Read XOPT.
            xopt(:) = xpt(:, kopt);
            xoptsq = linalg_obj.inprod(xopt, xopt);

            % Update BMAT. See (7.11)--(7.12) of the NEWUOA paper and the elaborations around.
            % XPTXAV corresponds to XPT - XAV in the NEWUOA paper, with XAV = (X0 + XOPT)/2.
            xptxav(:, :) = xpt - consts_obj.HALF * fortran.spread(xopt, 'dim', 2, 'ncopies', npt);
            %%MATLAB: xptxav = xpt - xopt/2  % xopt should be a column! Implicit expansion
            %sxpt = matprod(xopt, xptxav)
            sxpt(:) = linalg_obj.matprod12(xopt, xpt) - consts_obj.HALF * xoptsq; % This one seems to work better numerically.

            % First, make the changes to BMAT that do not depend on ZMAT.
            qxoptq = consts_obj.QUART * xoptsq;
            for k = 1:npt
                ymat(:, k) = sxpt(k) * xptxav(:, k) + qxoptq * xopt;
            end
            %%MATLAB: ymat = xptxav .* sxpt + qxoptq * xopt  % sxpt should be a row, xopt should be a column
            %ymat(:, kopt) = HALF * xoptsq * xopt ! This makes no difference according to a test on 20220406
            bymat(:, :) = linalg_obj.matprod22(bmat(:, 1:npt), ymat.'); % BMAT(:, 1:NPT) is not updated yet.
            bmat(:, npt + 1:npt + n) = bmat(:, npt + 1:npt + n) + (bymat + bymat.');
            % Then the revisions of BMAT that depend on ZMAT are calculated.
            yzmat(:, :) = linalg_obj.matprod22(ymat, zmat);
            yzmat_c(:, :) = yzmat;
            yzmat_c(:, 1:idz_loc - 1) = -yzmat(:, 1:idz_loc - 1); % IDZ_LOC is usually small. So this assignment is cheap.
            bmat(:, npt + 1:npt + n) = bmat(:, npt + 1:npt + n) + linalg_obj.matprod22(yzmat, yzmat_c.');
            bmat(:, 1:npt) = bmat(:, 1:npt) + linalg_obj.matprod22(yzmat_c, zmat.');

            % Update the quadratic model. Note that PQ remains unchanged. For HQ, see (7.14) of the NEWUOA paper.
            %v = matprod(xptxav, pq)  ! Vector V in (7.14) of the NEWUOA paper
            v(:) = linalg_obj.matprod21(xpt, pq) - consts_obj.HALF * sum(pq, 'all') * xopt; % This one seems to work better numerically.
            vxopt(:, :) = linalg_obj.outprod(v, xopt); %%MATLAB: vxopt = v * xopt';  % v and xopt should be both columns
            hq(:, :) = (vxopt + vxopt.') + hq; %call r2update(hq, ONE, xopt, v)
            %call symmetrize(hq)  ! Do this if the update above does not ensure symmetry.

            % The following instructions complete the shift of XBASE.
            xbase(:) = xbase + xopt;
            xpt(:, :) = xpt - xopt;
            xpt(:, kopt) = consts_obj.ZERO;
            %%MATLAB: xpt = xpt - xopt; xpt(:, kopt) = 0;  % xopt should be a column! Implicit expansion

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(xbase) == n && all(infnan_obj.is_finite(xbase), 'all'), "SIZE(XBASE) == N, XBASE is finite", srname);
                debug_obj.assert(size(xpt, 1) == n && size(xpt, 2) == npt, "SIZE(XPT) == [N, NPT]", srname);
                debug_obj.assert(all(abs(xpt(:, kopt)) <= 0, 'all'), "XPT(:, KOPT) == 0", srname);
                debug_obj.assert(all(infnan_obj.is_finite(xpt), 'all'), "XPT is finite", srname);
                debug_obj.assert(size(bmat, 1) == n && size(bmat, 2) == npt + n, "SIZE(BMAT)==[N, NPT+N]", srname);
                debug_obj.assert(linalg_obj.issymmetric(bmat(:, npt + 1:npt + n)), "BMAT(:, NPT+1:NPT+N) is symmetric", srname);
                debug_obj.assert(size(hq, 1) == n && linalg_obj.issymmetric(hq), "HQ is an NxN symmetric matrix", srname);
                % The following test cannot be passed.
                %call assert(errh(idz_loc, bmat, zmat, xpt) <= htol, 'H = W^{-1} in (3.12) of the NEWUOA paper', srname)

            end

        end
        function [pl, pq, xbase, xpt] = shiftbase_qint(~, kopt, pl, pq, xbase, xpt)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine shifts the base point from XBASE to XBASE + XOPT, and make the corresponding
            % changes to the gradients of the Lagrange functions and the quadratic model. See the discussion
            % below (40) of the UOBYQA paper.
            %--------------------------------------------------------------------------------------------------%

            % Common modules
            consts_obj = prima_mat.common.consts_mod();
            debug_obj = prima_mat.common.debug_mod();
            infnan_obj = prima_mat.common.infnan_mod();
            linalg_obj = prima_mat.common.linalg_mod();


            % Inputs


            % In-outputs
            % XBASE(N)
            % XPT(N, NPT)
            % PL(NPT-1, NPT)
            % PQ(NPT-1)

            % Local variables
            srname = "SHIFTBASE_QINT";


            xopt = NaN(numel(xbase), 1);

            % Sizes
            n = size(xpt, 1);
            npt = size(xpt, 2);

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(npt == (n + 1) * (n + 2) / 2, "NPT = (N+1)(N+2)/2", srname);
                debug_obj.assert(kopt >= 1 && kopt <= npt, "1 <= KOPT <= NPT", srname);
                debug_obj.assert(size(pl, 1) == npt - 1 && size(pl, 2) == npt, "SIZE(PL) == [NPT-1, NPT]", srname);
                debug_obj.assert(numel(pq) == npt - 1, "SIZE(PQ) == NPT-1", srname);
                debug_obj.assert(numel(xbase) == n && all(infnan_obj.is_finite(xbase), 'all'), "SIZE(XOPT) == N, XOPT is finite", srname);
                debug_obj.assert(size(xpt, 1) == n && size(xpt, 2) == npt && all(infnan_obj.is_finite(xpt), 'all'), "SIZE(XPT) == [N, NPT], XPT is finite", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            % Shift the base point from XBASE to XBASE + XOPT.
            xopt(:) = xpt(:, kopt);
            xbase(:) = xbase + xopt;
            xpt(:, :) = xpt - xopt;
            xpt(:, kopt) = consts_obj.ZERO;

            % Update the gradient of the model
            pq(1:n) = pq(1:n) + linalg_obj.smat_mul_vec(pq(n + 1:npt - 1), xopt);

            % Update the gradient of the Lagrange functions.
            for k = 1:npt
                pl(1:n, k) = pl(1:n, k) + linalg_obj.smat_mul_vec(pl(n + 1:npt - 1, k), xopt);
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(size(pl, 1) == npt - 1 && size(pl, 2) == npt, "SIZE(PL) == [NPT-1, NPT]", srname);
                debug_obj.assert(numel(pq) == npt - 1, "SIZE(PQ) == NPT-1", srname);
                debug_obj.assert(numel(xbase) == n && all(infnan_obj.is_finite(xbase), 'all'), "SIZE(XBASE) == N, XBASE is finite", srname);
                debug_obj.assert(size(xpt, 1) == n && size(xpt, 2) == npt && all(infnan_obj.is_finite(xpt), 'all'), "SIZE(XPT) == [N, NPT], XPT is finite", srname);
            end

        end

    end
end
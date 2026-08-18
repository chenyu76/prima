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



            % PQ(NPT)
            % ZMAT(NPT, NPT - N - 1)
            % Absent in BOBYQA, being equivalent to IDZ = 1


            % BMAT(N, NPT + N)
            % HQ(N, N)
            % XBASE(N)
            % XPT(N, NPT)



            bymat = NaN(numel(xbase));
            %real(RP) :: htol

            sxpt = NaN(size(xpt, 2), 1);
            v = NaN(numel(xbase), 1);
            vxopt = NaN(numel(xbase));
            xopt = NaN(numel(xbase), 1);

            xptxav = NaN(size(xpt, 1), size(xpt, 2));
            ymat = NaN(size(xpt, 1), size(xpt, 2));
            yzmat = NaN(numel(xbase), size(zmat, 2));


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


            %====================%
            % Calculation starts %
            %====================%

            % Read XOPT.
            xopt(:) = xpt(:, kopt);
            xoptsq = sum(xopt .* xopt, 'all');

            % Update BMAT. See (7.11)--(7.12) of the NEWUOA paper and the elaborations around.
            % XPTXAV corresponds to XPT - XAV in the NEWUOA paper, with XAV = (X0 + XOPT)/2.
            xptxav(:, :) = xpt - 0.5 * xopt;
            %%MATLAB: xptxav = xpt - xopt/2  % xopt should be a column! Implicit expansion
            %sxpt = matprod(xopt, xptxav)
            sxpt(:) = xpt.' * xopt - 0.5 * xoptsq; % This one seems to work better numerically.

            % First, make the changes to BMAT that do not depend on ZMAT.
            qxoptq = 0.25 * xoptsq;
            for k = 1:npt
                ymat(:, k) = sxpt(k) * xptxav(:, k) + qxoptq * xopt;
            end
            %%MATLAB: ymat = xptxav .* sxpt + qxoptq * xopt  % sxpt should be a row, xopt should be a column
            %ymat(:, kopt) = HALF * xoptsq * xopt ! This makes no difference according to a test on 20220406
            bymat(:, :) = bmat(:, 1:npt) * ymat.'; % BMAT(:, 1:NPT) is not updated yet.
            bmat(:, npt + 1:npt + n) = bmat(:, npt + 1:npt + n) + (bymat + bymat.');
            % Then the revisions of BMAT that depend on ZMAT are calculated.
            yzmat(:, :) = ymat * zmat;
            yzmat_c = yzmat;
            yzmat_c(:, 1:idz_loc - 1) = -yzmat(:, 1:idz_loc - 1); % IDZ_LOC is usually small. So this assignment is cheap.
            bmat(:, npt + 1:npt + n) = bmat(:, npt + 1:npt + n) + yzmat * yzmat_c.';
            bmat(:, 1:npt) = bmat(:, 1:npt) + yzmat_c * zmat.';

            % Update the quadratic model. Note that PQ remains unchanged. For HQ, see (7.14) of the NEWUOA paper.
            %v = matprod(xptxav, pq)  ! Vector V in (7.14) of the NEWUOA paper
            v(:) = xpt * pq - 0.5 * sum(pq, 'all') * xopt; % This one seems to work better numerically.
            vxopt(:, :) = v * xopt.'; %%MATLAB: vxopt = v * xopt';  % v and xopt should be both columns
            hq(:, :) = (vxopt + vxopt.') + hq; %call r2update(hq, ONE, xopt, v)
            %call symmetrize(hq)  ! Do this if the update above does not ensure symmetry.

            % The following instructions complete the shift of XBASE.
            xbase(:) = xbase + xopt;
            xpt(:, :) = xpt - xopt;
            xpt(:, kopt) = 0.0;
            %%MATLAB: xpt = xpt - xopt; xpt(:, kopt) = 0;  % xopt should be a column! Implicit expansion

            %====================%
            %  Calculation ends  %
            %====================%



        end
        function [pl, pq, xbase, xpt] = shiftbase_qint(~, kopt, pl, pq, xbase, xpt)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine shifts the base point from XBASE to XBASE + XOPT, and make the corresponding
            % changes to the gradients of the Lagrange functions and the quadratic model. See the discussion
            % below (40) of the UOBYQA paper.
            %--------------------------------------------------------------------------------------------------%



            linalg_obj = prima_mat.common.linalg_mod();


            % XBASE(N)
            % XPT(N, NPT)
            % PL(NPT-1, NPT)
            % PQ(NPT-1)



            xopt = NaN(numel(xbase), 1);


            n = size(xpt, 1);
            npt = size(xpt, 2);


            %====================%
            % Calculation starts %
            %====================%

            % Shift the base point from XBASE to XBASE + XOPT.
            xopt(:) = xpt(:, kopt);
            xbase(:) = xbase + xopt;
            xpt(:, :) = xpt - xopt;
            xpt(:, kopt) = 0.0;

            % Update the gradient of the model
            pq(1:n) = pq(1:n) + linalg_obj.smat_mul_vec(pq(n + 1:npt - 1), xopt);

            % Update the gradient of the Lagrange functions.
            for k = 1:npt
                pl(1:n, k) = pl(1:n, k) + linalg_obj.smat_mul_vec(pl(n + 1:npt - 1, k), xopt);
            end

            %====================%
            %  Calculation ends  %
            %====================%



        end

    end
end
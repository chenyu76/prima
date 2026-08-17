classdef rescue_mod
    %--------------------------------------------------------------------------------------------------%
    % This module contains the RESCUE subroutine described in Section 5 of the BOBYQA paper.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the BOBYQA paper.
    %
    % N.B.:
    % 1. According to a test on 20220425, the invocations of RESCUE is rare --- it is never invoked
    % on CUTEst unconstrained or bound constrained problems with at most 50 variables unless heavy noise
    % is imposed on the function evaluation.
    % 2. Zaikun (20230321): According to a test on 20230321 on problems of at most 200 variables, it
    % affects (not necessarily worsen) the performance of BOBYQA quite marginally if RESCUE is
    % completely disabled. Therefore, in the first implementation of an algorithm based on the
    % derivative-free PSB, it seems safe to ignore RESCUE. It is similar for the IDZ technique of
    % NEWUOA/LINCOA.
    %
    % Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
    %
    % Started: February 2022
    %
    % Last Modified: Thu 14 Aug 2025 07:34:53 AM CST
    %--------------------------------------------------------------------------------------------------%

    methods
        function [kopt, nf, fhist, fval, gopt, hq, pq, sl, su, xbase, xhist, xpt, bmat, zmat, info] = rescue(obj, calfun, solver, iprint, maxfun, delta, ftarget, xl, xu, kopt, nf, fhist, fval, gopt, hq, pq, sl, su, xbase, xhist, xpt, bmat, zmat)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine implements "the method of RESCUE" introduced in Section 5 of BOBYQA paper. The
            % purpose of this subroutine is to replace a few interpolation points by new points in order to
            % improve the geometry of the interpolation set and the conditioning of the interpolation system.
            % This is done in the following way.
            %
            % 1. Define a set of "provisional interpolation points" XPT_PROV around the current XOPT. Similar to
            % the construction of the initial interpolation set, XPT_PROV is obtained by perturbing XOPT subject
            % to the bound constraints along one or two coordinate directions, the latter taking place only if
            % NPT >= 2*N+2. See (5.4)--(5.5) of the BOBYQA paper for details. After defining XPT_PROV, set BMAT
            % and ZMAT to represent the H matrix defined in (2.7) of the BOBYQA paper corresponding to XPT_PROV.
            % N.B.: In the code, XPT_PROV is not formed explicitly, but represented implicitly by XPT, PTSID,
            % and PTSAUX.
            % 2. For each "original interpolation point" in XPT, check whether it can replace a point in XPT_PROV
            % without damaging the geometry of XPT_PROV, which is indicated by the denominator SIGMA in the
            % updating formula of H due to the replacement (see (4.9) of the BOBYQA paper). If yes, update
            % XPT_PROV by the performing such an replacement. Continue doing this until all the original points
            % get reinstated in XPT_PROV, or we cannot find any original point that can safely replace
            % a provisional point. Note the following.
            % 2.1. Suppose that the KORIG-th original point is going to replace the KPROV-th provisional point.
            % Then we first exchange the KORIG-th and KPROV-th provisional points, and then replace the new
            % KORIG-th provisional point by the KORIG-th original point. BMAT and ZMAT are updated accordingly.
            % In this way, FVAL(KORIG) is consistent with XPT_PROV(:, KORIG), so that we need not update FVAL.
            % 2.1. The KOPT-th original point (i.e., XOPT) is always reinstated in XPT_PROV. This is done by
            % replacing XPT_PROV(:, 1) with XOPT. This is the first replacement to perform.
            % 2.2. After XOPT is reinstated, the original points are ranked according the scores saved in SCORE.
            % SCORE is initialized to the squares of the distances from between the original point and XOPT.
            % The KORIG is set to the index of the point with the smallest positive score. If XPT(:, KORIG)
            % cannot replace any provisional point safely, then set SCORE(KORIG) to -SCORE(KORIG) - SCOREINC;
            % otherwise, set SCORE(KORIG) to 0 and all the other scores to their absolute values. In this way,
            % an original point that fails to replace any provisional point during the current attempt will be
            % skipped until another original point succeeds in doing so; moreover, the failing point will get
            % a lower priority in later attempts. An original point that successfully replaces a provisional
            % point will not be tried again due to the zero score.
            % 2.3. Once a provisional point is replaced, it will be marked (by setting the corresponding entry
            % of PTSID to zero) so that it will not be replaced again.
            % 3. When the above procedure finishes, normally most original points are reinstated in XPT_PROV, so
            % that XPT_PROV differs from XPT only at very few positions. Set XPT to XPT_PROV, update FVAL at the
            % new interpolation points by evaluating F, and then quadratic interpolant [GQ, PQ, HQ] accordingly.
            %
            % At the end of the subroutine, the elements of BMAT and ZMAT are set in a well-conditioned way to
            % the values that are appropriate for the new interpolation points. The elements of GOPT, HQ and PQ
            % are also revised to the values that are appropriate to the final quadratic model.
            %
            % The arguments NF, KOPT, XL, XU, IPRINT, MAXFUN, XBASE, XPT, FVAL GOPT, HQ, PQ, BMAT, ZMAT, SL
            % and, SU have the same meanings as the corresponding arguments of BOBYQB on the entry to RESCUE.
            % DELTA is the current trust region radius.
            % PTSAUX is a 2-by-N real array. For J = 1, 2, ..., N, PTSAUX(1, J) and PTSAUX(2, J) specify the two
            %   positions of provisional interpolation points when a nonzero step is taken along e_J (the J-th
            %   coordinate direction) through XBASE + XOPT, as specified below. Usually these steps have length
            %   DELTA, but other lengths are chosen if necessary in order to satisfy the bound constraints.
            % PTSID is an integer array of length NPT. Its components denote provisional new positions of the
            %   interpolation points. The K-th point is a candidate for change if and only if PTSID(K) is
            %   nonzero. In this case let IP and IQ be the integer parts of PTSID(K) and (PTSID(K)-IP)*(N+1). If
            %   IP and IQ are both positive, the step from XBASE + XOPT to the new K-th interpolation point is
            %   PTSAUX(1, IP)*e_IP + PTSAUX(1, IQ)*e_IQ. Otherwise the step is either PTSAUX(1, IP)*e_IP or
            %   PTSAUX(2, IQ)*e_IQ in the cases IQ=0 or  IP=0, respectively.
            %--------------------------------------------------------------------------------------------------%

            % Common modules
            checkexit_obj = prima_mat.common.checkexit_mod();
            consts_obj = prima_mat.common.consts_mod();
            debug_obj = prima_mat.common.debug_mod();
            evaluate_obj = prima_mat.common.evaluate_mod();
            history_obj = prima_mat.common.history_mod();
            infnan_obj = prima_mat.common.infnan_mod();
            infos_obj = prima_mat.common.infos_mod();
            linalg_obj = prima_mat.common.linalg_mod();
            message_obj = prima_mat.common.message_mod();

            powalg_obj = prima_mat.common.powalg_mod();
            string_obj = prima_mat.common.string_mod();
            xinbd_obj = prima_mat.common.xinbd_mod();


            % Inputs



            % XL(N)
            % XU(N)

            % In-outputs


            % FHIST(MAXFHIST)
            % FVAL(NPT)
            % GOPT(N)
            % HQ(N, N)
            % PQ(NPT)
            % SL(N)
            % SU(N)
            % XBASE(N)
            % XHIST(N, MAXXHIST)
            % XPT(N, NPT)

            % Outputs

            %  BMAT(N, NPT + N)
            % ZMAT(NPT, NPT-N-1)

            % Local variables
            srname = "RESCUE";
            ij = NaN(2, max(0, size(xpt, 2) - 2 * size(xpt, 1) - 1));
            ip = NaN;
            iq = NaN;


            k = NaN;

            korig = NaN;
            kprov = NaN;


            subinfo = NaN;
            mask = false(size(xpt, 1), 1);
            beta = NaN;
            bsum = NaN;
            den = NaN(size(xpt, 2), 1);
            f = NaN;

            hcol = NaN(size(bmat, 2), 1);
            hdiag = NaN(size(xpt, 2), 1);
            moderr = NaN;
            pqinc = NaN(size(xpt, 2), 1);
            ptsaux = NaN(2, size(xpt, 1));
            ptsid = NaN(size(xpt, 2), 1);
            score = NaN(size(xpt, 2), 1);


            v = NaN(size(xpt, 1), 1);
            vlag = NaN(size(xpt, 1) + size(xpt, 2), 1);
            vquad = NaN;
            wmv = NaN(size(xpt, 1) + size(xpt, 2), 1);
            x = NaN(size(xpt, 1), 1);
            xnew = NaN(size(xpt, 1), 1);
            xopt = NaN(size(xpt, 1), 1);
            xp = NaN;
            xq = NaN;
            xxpt = NaN(size(xpt, 2), 1);

            n = size(xpt, 1);
            npt = size(xpt, 2);
            maxxhist = size(xhist, 2);
            maxfhist = numel(fhist);
            maxhist = max(maxxhist, maxfhist);

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(abs(iprint) <= 3, "IPRINT is 0, 1, -1, 2, -2, 3, or -3", srname);
                debug_obj.assert(n >= 1, "N >= 1", srname);
                debug_obj.assert(npt >= n + 2, "NPT >= N+2", srname);
                debug_obj.assert(maxfun >= npt + 1, "MAXFUN >= NPT+1", srname);
                debug_obj.assert(kopt >= 1 && kopt <= npt, "1 <= KOPT <= NPT", srname);
                debug_obj.assert(delta > 0, "DELTA > 0", srname);
                debug_obj.assert(numel(fval) == npt && ~any(infnan_obj.is_nan_sp(fval) | infnan_obj.is_posinf(fval), 'all'), "SIZE(FVAL) == NPT and FVAL is not NaN/+Inf", srname);
                debug_obj.assert(~any(fval < fval(kopt), 'all'), "FVAL(KOPT) is the smallest in FVAL", srname);
                debug_obj.assert(maxfhist * (maxfhist - maxhist) == 0, "SIZE(FHIST) == 0 or MAXHIST", srname);
                debug_obj.assert(numel(xl) == n && numel(xu) == n, "SIZE(XL) == N == SIZE(XU)", srname);
                debug_obj.assert(numel(sl) == n && all(sl <= 0, 'all'), "SIZE(SL) == N, SL <= 0", srname);
                debug_obj.assert(numel(su) == n && all(su >= 0, 'all'), "SIZE(SU) == N, SU >= 0", srname);
                debug_obj.assert(numel(gopt) == n, "SIZE(GOPT) == N", srname);
                debug_obj.assert(size(hq, 1) == n && linalg_obj.issymmetric(hq), "HQ is n-by-n and symmetric", srname);
                debug_obj.assert(numel(pq) == npt, "SIZE(PQ) == NPT", srname);
                debug_obj.assert(numel(xbase) == n && all(infnan_obj.is_finite(xbase), 'all'), "SIZE(XBASE) == N, XBASE is finite", srname);
                debug_obj.assert(all(xbase >= xl & xbase <= xu, 'all'), "XL <= XBASE <= XU", srname);
                debug_obj.assert(size(xhist, 1) == n && maxxhist * (maxxhist - maxhist) == 0, "SIZE(XHIST, 1) == N, SIZE(XHIST, 2) == 0 or MAXHIST", srname);
                debug_obj.assert(all(infnan_obj.is_finite(xhist(:, 1:min(nf, maxxhist))), 'all'), "XHIST is finite", srname);
                for k = 1:min(nf, maxxhist)
                    debug_obj.assert(all(xhist(:, k) >= xl, 'all') && all(xhist(:, k) <= xu, 'all'), "XL <= XHIST <= XU", srname);
                end
                debug_obj.assert(all(infnan_obj.is_finite(xpt), 'all'), "XPT is finite", srname);
                debug_obj.assert(all(xpt >= reshape(sl, [], 1), 'all') && all(xpt <= reshape(su, [], 1), 'all'), "SL <= XPT <= SU", srname);
                debug_obj.assert(size(bmat, 1) == n && size(bmat, 2) == npt + n, "SIZE(BMAT) == [N, NPT+N]", srname);
                debug_obj.assert(size(zmat, 1) == npt && size(zmat, 2) == npt - n - 1, "SIZE(ZMAT) == [NPT, NPT-N-1]", srname);
                debug_obj.assert(maxhist >= 0 && maxhist <= maxfun, "0 <= MAXHIST <= MAXFUN", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            info = infos_obj.INFO_DFT;

            % Do nothing if NF already reaches it upper bound.
            % To please Fortran compilers, set BMAT and ZMAT before returning, though they will not be used.
            if nf >= maxfun
                bmat = repmat(consts_obj.ZERO, size(bmat));
                zmat = repmat(consts_obj.ZERO, size(zmat));
                info = infos_obj.MAXFUN_REACHED;
                return
            end

            % Shift the interpolation points so that XOPT becomes the origin.
            xopt(:) = xpt(:, kopt);
            sl(:) = min(sl - xopt, consts_obj.ZERO);
            su(:) = max(su - xopt, consts_obj.ZERO);
            xbase(:) = min(max(xl, xbase + xopt), xu);
            xpt(:, :) = xpt - reshape(xopt, [], 1);
            xpt(:, kopt) = consts_obj.ZERO;

            % Update HQ so that HQ and PQ define the second derivatives of the model after XBASE has been
            % shifted to the trust region centre.
            v(:) = linalg_obj.matprod21(xpt, pq) + consts_obj.HALF * sum(pq, 'all') * xopt;
            hq = linalg_obj.r2_sym(hq, consts_obj.ONE, xopt, v);

            % Set the elements of PTSAUX.
            ptsaux(1, :) = min(delta, su);
            ptsaux(2, :) = max(-delta, sl);
            mask(:) = (ptsaux(1, :) + ptsaux(2, :) < 0);
            ptsaux([1, 2], linalg_obj.trueloc(mask)) = ptsaux([2, 1], linalg_obj.trueloc(mask));
            mask(:) = (abs(ptsaux(2, :)) < consts_obj.HALF * abs(ptsaux(1, :)));
            ptsaux(2, linalg_obj.trueloc(mask)) = consts_obj.HALF * ptsaux(1, linalg_obj.trueloc(mask));

            % Set the identifiers of the artificial interpolation points that are along a coordinate direction
            % from XOPT, and set the corresponding nonzero elements of BMAT and ZMAT.
            sfrac = consts_obj.HALF / double(n + 1);
            ptsid(1) = sfrac;
            bmat = repmat(consts_obj.ZERO, size(bmat));
            zmat = repmat(consts_obj.ZERO, size(zmat));
            for k = 1:n
                ptsid(k + 1) = double(k) + sfrac;
                if k <= npt - n - 1
                    ptsid(k + n + 1) = double(k) / double(n + 1) + sfrac;
                    temp = consts_obj.ONE / (ptsaux(1, k) - ptsaux(2, k));
                    bmat(k, k + 1) = -temp + consts_obj.ONE / ptsaux(1, k);
                    bmat(k, k + n + 1) = temp + consts_obj.ONE / ptsaux(2, k);
                    bmat(k, 1) = -bmat(k, k + 1) - bmat(k, k + n + 1);
                    zmat(1, k) = sqrt(consts_obj.TWO) / abs(ptsaux(1, k) * ptsaux(2, k));
                    zmat(k + 1, k) = zmat(1, k) * ptsaux(2, k) * temp;
                    zmat(k + n + 1, k) = -zmat(1, k) * ptsaux(1, k) * temp;
                else
                    bmat(k, 1) = -consts_obj.ONE / ptsaux(1, k);
                    bmat(k, k + 1) = consts_obj.ONE / ptsaux(1, k);
                    bmat(k, k + npt) = -consts_obj.HALF * ptsaux(1, k) ^ 2;
                end
            end

            % Set any remaining identifiers with their nonzero elements of ZMAT.
            ij(:, :) = powalg_obj.setij(n, npt);
            for k = 2 * n + 2:npt
                ip = ij(1, k - 2 * n - 1);
                iq = ij(2, k - 2 * n - 1);
                ptsid(k) = double(ip) + double(iq) / double(n + 1) + sfrac;
                temp = consts_obj.ONE / (ptsaux(1, ip) * ptsaux(1, iq));
                zmat([1, k], k - n - 1) = temp;
                zmat([ip + 1, iq + 1], k - n - 1) = -temp;
            end

            % Update BMAT, ZMAT, ans PTSID so that the 1st and the KOPT-th provisional points are exchanged.
            % After the exchanging, the KOPT-th point in the provisional set becomes the zero vector, which is
            % exactly the KOPT-th original point (after the shift of XBASE at the beginning of the subroutine).
            if kopt ~= 1
                bmat(:, [1, kopt]) = bmat(:, [kopt, 1]);
                zmat([1, kopt], :) = zmat([kopt, 1], :);
            end
            ptsid(1) = ptsid(kopt);
            ptsid(kopt) = consts_obj.ZERO;

            % The squares of the distances from XOPT to the other interpolation points are set at SCORE, which
            % will be used to define the index KORIG in the loop below.  Increments of SCOREINC may be added
            % later to these scores to balance the consideration of the choice of point that is going to become
            % current. Note that, in Powell's BOBYQA code, the initial scores are the squares of the distances,
            % but there is no square in the BOBYQA paper (see the paragraph between (5.9) and (5.10) of the
            % BOBYQA paper). The latter seem to work better in a test on 20221125.
            %score = sum(xpt**2, dim=1)  ! Powell's BOBYQA code
            score(:) = sqrt(sum(xpt .^ 2, 1)); % Powell's BOBYQA paper
            % In theory, SCORE(KOPT) = 0. Make sure this so that KOPT will be skipped when we choose KORIG below.
            score(kopt) = consts_obj.ZERO;
            scoreinc = max(score, [], 'all');

            % NPROV is the number of provisional points that has not yet been replaced with original points.
            nprov = npt - 1;

            % Even without an upper bound for the loop counter, the following loop runs for at most NPT^2 times:
            % for each value of NPROV, we need at most NPT loops to find an original point that can safely
            % replace a provisional point; if such a pair of origin and provisional points are found, then NPROV
            % will de reduced by 1; otherwise, SCORE will become all zero or negative, and the loop will exit.
            % Originally, it is a WHILE loop, but we change it to a DO loop to avoid infinite cycling.
            % N.B.: Overflow will occur in NPT^2 if NPT > 180 and IK = 16. The following is a workaround, which
            % is **not needed in Python/MATLAB/Julia/R. In MATLAB, we can just take maxiter = npt^2**.
            maxiter = fix(min(10 ^ min(floor(log10(double(intmax('int64')))), floor(log10(double(intmax('int64'))))), npt ^ 2)); %%MATLAB: maxiter = npt^2;
            for iter = 1:maxiter
                % %DO WHILE (ANY(SCORE > 0) .AND. NPROV > 1)   ! WHILE version.
                % %IF (ALL(SCORE <= 0) .AND. NPROV <= 0) THEN ! Powell's code. May not take any provisional point.
                % %IF (ALL(SCORE <= 0) .AND. NPROV <= 2) THEN  ! Retain at least two provisional points.
                if all(score <= 0, 'all') || nprov <= 1
                    % Retain at least one provisional point.
                    break
                end

                % Pick the index KORIG of an original point that has not yet replaced one of the provisional
                % points, giving attention to the closeness to XOPT and to previous tries with KORIG.
                korig = fortran.minloc(score, 'mask', (score > 0), 'dim', 1);

                % Calculate VLAG and BETA for the required updating of the H matrix if XPT(:, KORIG) is
                % reinstated in the set of interpolation points, which means to replace a point in the
                % following provisional interpolation set XPT_PROV defined in (5.4)--(5.5) of the BOBYQA paper.
                % 1. XPT_PROV(:, KOPT) = 0;
                % 2. For each K /= KOPT, if PTSID(K) == 0, then XPT_PROV(:, K) = XPT(:, K); if PTSID(K) > 0,
                % then XPT_PROV(:, K) has nonzeros only at IP (if IP > 0), IQ (if IQ > 0) positions, where IP
                % and IQ are the P(J) and Q(J) defined in and below (2.4) of the BOBYQA paper.

                % First, form the (W - V) vector for XPT(:, KORIG).
                % In the code below, WMV = W - V = w(XNEW) - w(XOPT) without the (NPT+1)th entry, where
                % XNEW = XPT(:, KORIG), XOPT = XPT_PROV = 0, and w(.) is defined by (6.3) of the NEWUOA paper
                % (as well as (4.10) of the BOBYQA paper). Since XOPT = 0, we see from (6.3) that w(XOPT)(K) = 0
                % for all K except that w(XOPT)(NPT+1) = 1. Thus WMV is the same at w(XNEW) without the
                % (NPT+1)-th entry. Therefore, WMV= [HALF*MATPROD(XNEW, XPT_PROV)**2, XNEW].
                for k = 1:npt
                    if k == kopt
                        wmv(k) = consts_obj.ZERO;
                    elseif ptsid(k) <= 0
                        % Indeed, PTSID >= 0. So PTSID(K) <= 0 means PTSID(K) = 0.
                        wmv(k) = linalg_obj.inprod(xpt(:, korig), xpt(:, k));
                    else
                        ip = floor(ptsid(k)); % IP = 0 if 0 < PTSID(K) < 1.
                        iq = floor(double(n + 1) * ptsid(k) - double((n + 1) * ip));
                        debug_obj.assert(ip >= 0 && ip <= npt && iq >= 0 && iq <= npt, "0 <= IP, IQ <= NPT", srname);
                        if ip > 0 && iq > 0
                            wmv(k) = xpt(ip, korig) * ptsaux(1, ip) + xpt(iq, korig) * ptsaux(1, iq);
                        elseif ip > 0
                            wmv(k) = xpt(ip, korig) * ptsaux(1, ip);
                        elseif iq > 0
                            wmv(k) = xpt(iq, korig) * ptsaux(2, iq);
                        else
                            wmv(k) = consts_obj.ZERO;
                        end
                    end
                    wmv(k) = consts_obj.HALF * wmv(k) * wmv(k);
                end
                wmv(npt + 1:npt + n) = xpt(:, korig);

                % Now calculate VLAG = H*WMV + e_KOPT according to (4.26) of the NEWUOA paper except VLAG(KOPT).
                vlag(1:npt) = linalg_obj.matprod21(zmat, linalg_obj.matprod12(wmv(1:npt), zmat)) + linalg_obj.matprod12(wmv(npt + 1:npt + n), bmat(:, 1:npt));
                vlag(npt + 1:npt + n) = linalg_obj.matprod21(bmat, wmv(1:npt + n));

                % Now calculate BETA. According to (4.12) of the NEWUOA paper (also (4.10) of the BOBYQA paper),
                % BETA = HALF*||XNEW - XOPT||^4 - WMV'*H*WMV. To calculate WMX'*H*WMV, note that
                % WMV'*H*WMV = WMV' * [Z*Z', B2^T; B1, B2] * WMV with Z = ZMAT, B1 = BMAT(:, 1:NPT), and
                % B2 = BMAT(:, NPT+1:NPT+N). Denoting W1 = WMV(1:NPT) and W2 = WMV(NPT+1:NPT+N), we then have
                % WMV'*H*WMV = ||W1'*Z||^2 + 2*W1'*B1*W2 + W1'*B2*W2 = ||W1'*Z||^2 + W1'(B1*W2 + [B1, B2]*WMV).
                bsum = linalg_obj.inprod(wmv(1:n), linalg_obj.matprod21(bmat(:, 1:npt), wmv(1:npt)) + linalg_obj.matprod21(bmat, wmv));
                beta = consts_obj.HALF * sum(xpt(:, korig) .^ 2, 'all') ^ 2 - sum(linalg_obj.matprod12(wmv(1:npt), zmat) .^ 2, 'all') - bsum;

                % Finally, set VLAG(KOPT) to the correct value.
                vlag(kopt) = vlag(kopt) + consts_obj.ONE;

                % For all K with PTSID(K) > 0, calculate the denominator DEN(K) = SIGMA in the updating formula
                % of H for XPT(:, KORIG) to replace XPT_PROV(:, K).
                den(:) = consts_obj.ZERO;
                hdiag(linalg_obj.trueloc(ptsid > 0)) = sum(zmat(linalg_obj.trueloc(ptsid > 0), :) .^ 2, 2);
                den(linalg_obj.trueloc(ptsid > 0)) = hdiag(linalg_obj.trueloc(ptsid > 0)) * beta + vlag(linalg_obj.trueloc(ptsid > 0)) .^ 2;

                % Attempt setting KPROV to the index of the provisional point to be replaced with the KORIG-th
                % original interpolation point. We choose KPROV by maximizing DEN(KPROV), which will be the
                % denominator SIGMA in the updating formula (4.9). In order to avoid a small denominator, we
                % consider it proper to replace the KPROV-th provisional point with the KORIG-th original point
                % only if DEN(KPROV) = MAXVAL(DEN) > C*MAXVAL(VLAG(1:NPT)**2), where C is a relatively small
                % positive constant --- C = 1 is achievable if the rounding errors were not severe; Powell took
                % C = 0.01, which prefers strongly the original point to the provisional (new) point, as the
                % latter necessitate new function evaluations. If this inequality is not achievable for the
                % current KORIG, then we will update SCORE(KORIG) to a negative value and continue the loop with
                % the next KORIG, which is set to MINLOC(SCORE, MASK=(SCORE > 0)) at the beginning of the loop,
                % skipping the original interpolation points with a nonpositive score. When a KORIG rendering
                % the aforesaid inequality is found, SCORE(KORIG) will be set to zero, and all the scores will
                % be reset to their absolute values, so that future attempts will try the original points that
                % have not succeeded in replacing a provisional point. The update of SCORE reflects an adaptive
                % ranking of the original points: points that are closer to XOPT have higher priority, and a
                % point will be ranked lower if it fails to fulfill MAXVAL(DEN) > C*MAXVAL(VLAG(1:NPT)**2).
                % Even if KORIG cannot satisfy this condition for now, it may validate the inequality in future
                % attempts, as BMAT and ZMAT will be updated.
                if ~(infnan_obj.is_finite(sum(abs(vlag), 'all')) && any(den > 5.0e-2 * max(vlag(1:npt) .^ 2, [], 'all'), 'all'))
                    % The above condition works a bit better than Powell's version below due to the factor 0.05.
                    % %IF (.NOT. (ANY(DEN > 1.0E-2_RP * MAXVAL(VLAG(1:NPT)**2)))) THEN  ! Powell' code
                    score(korig) = -score(korig) - scoreinc;
                    continue
                end
                kprov = fortran.maxloc(den, 'mask', (~infnan_obj.is_nan_sp(den)), 'dim', 1);
                %%MATLAB: [~, kprov] = max(den, [], 'omitnan');

                % Update BMAT, ZMAT, VLAG, and PTSID to exchange the KPROV-th and KORIG-th provisional points.
                % After the exchanging, the KORIG-th original point will replace the KORIG-th provisional point.
                if kprov ~= korig
                    bmat(:, [kprov, korig]) = bmat(:, [korig, kprov]);
                    zmat([kprov, korig], :) = zmat([korig, kprov], :);
                    vlag([kprov, korig]) = vlag([korig, kprov]);
                end
                ptsid(kprov) = ptsid(korig);

                % Set PTSID(KORIG) = 0 so that the KORIG-th provisional point (after the exchanging) will be
                % skipped in the later loops.
                ptsid(korig) = consts_obj.ZERO;
                % Set SCORE(KORIG) = 0 so that the KORIG-th original point will be skipped in later loops.
                score(korig) = consts_obj.ZERO;
                % Reset SCORE to ABS(SCORE) so that all the original points with a nonzero score will be checked
                % in later loops.
                score(:) = abs(score);

                % Update the BMAT and ZMAT matrices so that the KORIG-th original point replaces the KORIG-th
                % provisional point.
                [bmat, zmat] = obj.updateh_rsc(korig, beta, vlag, bmat, zmat);

                % NPROV is the number of provisional points that has not yet been replaced with original points.
                nprov = nprov - 1;
            end

            % All the final positions of the interpolation points have been chosen although any changes have not
            % been included yet in XPT. Also the final BMAT and ZMAT matrices are complete, but, apart from the
            % shift of XBASE, the updating of the quadratic model remains to be done. The following cycle
            % through the new interpolation points begins by putting the new point in XPT(:, KPT) and by setting
            % PQ(KPT) to zero. A return occurs if MAXFUN prohibits another value of F or when all the new
            % interpolation points are included in the model.
            kbase = kopt;
            fbase = fval(kopt);
            if nprov > 0
                for kpt = 1:npt
                    if ptsid(kpt) <= 0
                        continue
                    end

                    % Absorb PQ(KPT)*XPT(:, KPT)*XPT(:, KPT)^T into the explicit part of the Hessian of the
                    % quadratic model. Implement R1UPDATE properly so that it ensures HQ is symmetric.
                    hq = linalg_obj.r1_sym(hq, pq(kpt), xpt(:, kpt));
                    pq(kpt) = consts_obj.ZERO;

                    ip = floor(ptsid(kpt));
                    iq = floor(double(n + 1) * ptsid(kpt) - double((n + 1) * ip));

                    % Update XPT(:, KPT) to the new point. It contains at most two nonzeros XP and XQ at the IP
                    % and IQ entries.
                    xp = consts_obj.ZERO;
                    xq = consts_obj.ZERO;
                    xnew(:) = consts_obj.ZERO;
                    if ip > 0 && iq > 0
                        xp = ptsaux(1, ip);
                        xnew(ip) = xp;
                        xq = ptsaux(1, iq);
                        xnew(iq) = xq;
                    elseif ip > 0
                        % IP > 0, IQ == 0
                        xp = ptsaux(1, ip);
                        xnew(ip) = xp;
                    elseif iq > 0
                        % IP == 0, IQ > 0
                        xq = ptsaux(2, iq);
                        xnew(iq) = xq;
                    end

                    % Zaikun 20240314: Skip the new point if it is too close to XPT(:, KPT), the point to replace.
                    % Indeed, it may even happen that XNEW == XPT(:, KPT), which did occur when RP = REAL16 (half
                    % precision) and led to an infinite cycling, because RESCUE did not make any change to XPT,
                    % and later the algorithm decided to call RESCUE again with the same data. This was fixed by
                    % the skipping, and by terminating the algorithm if RESCUE is requested for two times
                    % without any new function evaluations in between, which was the behavior of Powell's code.
                    % Skipping an XNEW that is close but not identical to XPT(:, KPT) will cause discrepancy
                    % between [BMAT, ZMAT] and XPT, since the former has been updated, but it is not severe as
                    % the difference between XNEW and XPT(:, KPT) is tiny.
                    if sum(abs(xnew - xpt(:, kpt)), 'all') <= 1.0e-2 * delta || ~infnan_obj.is_finite(sum(abs(xnew), 'all'))
                        continue
                    end
                    xpt(:, kpt) = xnew;

                    % Calculate F at the new interpolation point, and set MODERR to the factor that is going to
                    % multiply the KPT-th Lagrange function when the model is updated to provide interpolation
                    % to the new function value.
                    x(:) = xinbd_obj.xinbd(xbase, xpt(:, kpt), xl, xu, sl, su); % In precise arithmetic, X = XBASE + XPT(:, KPT).
                    f = evaluate_obj.evaluatef(calfun, x);
                    nf = nf + 1;

                    % Print a message about the function evaluation according to IPRINT.
                    message_obj.fmsg(solver, "Rescue", iprint, nf, delta, f, x);
                    % Save X, F into the history.
                    [xhist, fhist] = history_obj.savehist(nf, x, xhist, f, fhist);

                    % Update FVAL and KOPT.
                    fval(kpt) = f;
                    if f < fval(kopt)
                        kopt = kpt;
                    end

                    % Check whether to exit
                    subinfo = checkexit_obj.checkexit_unc(maxfun, nf, f, ftarget, x);
                    if subinfo ~= infos_obj.INFO_DFT
                        info = subinfo;
                        break
                    end

                    % Set VQUAD to the value of the current model at the new XPT(:, KPT), which has at most two
                    % nonzeros XP and XQ at the IP and IQ entries respectively.
                    vquad = fbase;
                    if ip > 0 && iq > 0
                        vquad = vquad + xp * (gopt(ip) + consts_obj.HALF * xp * hq(ip, ip));
                        vquad = vquad + xq * (gopt(iq) + consts_obj.HALF * xq * hq(iq, iq));
                        vquad = vquad + xp * xq * hq(ip, iq);
                        xxpt(:) = xp * xpt(ip, :) + xq * xpt(iq, :);
                    elseif ip > 0
                        % IP > 0, IQ == 0
                        vquad = vquad + xp * (gopt(ip) + consts_obj.HALF * xp * hq(ip, ip));
                        xxpt(:) = xp * xpt(ip, :);
                    elseif iq > 0
                        % IP == 0, IQ > 0
                        vquad = vquad + xq * (gopt(iq) + consts_obj.HALF * xq * hq(iq, iq));
                        xxpt(:) = xq * xpt(iq, :);
                    end
                    vquad = vquad + consts_obj.HALF * linalg_obj.inprod(xxpt, pq .* xxpt);
                    % N.B.: INPROD(XXPT, PQ * XXPT) = INPROD(X, HESS_MUL(X, XPT, PQ))

                    % Update the quadratic model.
                    moderr = f - vquad;
                    gopt(:) = gopt + moderr * bmat(:, kpt);
                    pqinc(:) = moderr * linalg_obj.matprod21(zmat, zmat(kpt, :));
                    pq(linalg_obj.trueloc(ptsid <= 0)) = pq(linalg_obj.trueloc(ptsid <= 0)) + pqinc(linalg_obj.trueloc(ptsid <= 0));
                    for k = 1:npt
                        if ptsid(k) <= 0
                            continue
                        end
                        ip = floor(ptsid(k));
                        iq = floor(double(n + 1) * ptsid(k) - double((n + 1) * ip));
                        if ip > 0 && iq > 0
                            hq(ip, ip) = hq(ip, ip) + pqinc(k) * ptsaux(1, ip) ^ 2;
                            hq(iq, iq) = hq(iq, iq) + pqinc(k) * ptsaux(1, iq) ^ 2;
                            hq(ip, iq) = hq(ip, iq) + pqinc(k) * ptsaux(1, ip) * ptsaux(1, iq);
                            hq(iq, ip) = hq(ip, iq);
                        elseif ip > 0
                            % IP > 0, IQ == 0
                            hq(ip, ip) = hq(ip, ip) + pqinc(k) * ptsaux(1, ip) ^ 2;
                        elseif iq > 0
                            % IP == 0, IP > 0
                            hq(iq, iq) = hq(iq, iq) + pqinc(k) * ptsaux(2, iq) ^ 2;
                        end
                    end
                    ptsid(kpt) = consts_obj.ZERO;
                end
            end

            % Update GOPT if necessary.
            if kopt ~= kbase
                gopt(:) = gopt + powalg_obj.hess_mul(xpt(:, kopt), xpt, pq, 'hq', hq);
            end

            %--------------------------------------------------------------------------------------------------%
            % Zaikun 20221123: What if we rebuild the model? It seems to worsen the performance of BOBYQA. Why?
            % %hq = ZERO
            % %pq = omega_mul(1_IK, zmat, fval - fval(kopt))
            % %gopt = matprod(bmat(:, 1:npt), fval - fval(kopt)) + hess_mul(xpt(:, kopt), xpt, pq)
            %--------------------------------------------------------------------------------------------------%

            %--------------------------------------------------------------------------------------------------%
            % Zaikun 20221123: Shouldn't we correct the models using the new [BMAT, ZMAT]?!
            % In this way, we do not even need the quadratic model received by RESCUE is an interpolant.
            % %real(RP) :: qval(size(xpt, 2))
            % %qval = [(quadinc(xpt(:, k) - xpt(:, kopt), xpt, gopt, pq, hq), k=1, npt)]
            % %pq = pq + omega_mul(1_IK, zmat, fval - qval - fval(kopt))
            % %gopt = gopt + matprod(bmat(:, 1:npt), fval - qval - fval(kopt)) + hess_mul(xpt(:, kopt), xpt, pq)
            %--------------------------------------------------------------------------------------------------%

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(kopt >= 1 && kopt <= npt, "1 <= KOPT <= NPT", srname);
                debug_obj.assert(numel(fhist) == maxfhist, "SIZE(FHIST) == MAXFHIST", srname);
                debug_obj.assert(numel(fval) == npt && ~any(infnan_obj.is_nan_sp(fval) | infnan_obj.is_posinf(fval), 'all'), "SIZE(FVAL) == NPT and FVAL is not NaN/+Inf", srname);
                debug_obj.assert(~any(fval < fval(kopt), 'all'), "FVAL(KOPT) is the smallest in FVAL", srname);
                debug_obj.assert(numel(sl) == n && numel(su) == n, "SIZE(SL) == N == SIZE(SU)", srname);
                debug_obj.assert(numel(gopt) == n, "SIZE(GOPT) == N", srname);
                debug_obj.assert(size(hq, 1) == n && linalg_obj.issymmetric(hq), "HQ is n-by-n and symmetric", srname);
                debug_obj.assert(numel(pq) == npt, "SIZE(PQ) == NPT", srname);
                debug_obj.assert(numel(xbase) == n && all(infnan_obj.is_finite(xbase), 'all'), "SIZE(XBASE) == N, XBASE is finite", srname);
                debug_obj.assert(all(xbase >= xl & xbase <= xu, 'all'), "XL <= XBASE <= XU", srname);
                debug_obj.assert(size(xhist, 1) == n && maxxhist * (maxxhist - maxhist) == 0, "SIZE(XHIST, 1) == N, SIZE(XHIST, 2) == 0 or MAXHIST", srname);
                debug_obj.assert(~any(infnan_obj.is_nan_sp(xhist(:, 1:min(nf, maxxhist))), 'all'), "XHIST does not contain NaN", srname);
                % The last calculated X can be Inf (finite + finite can be Inf numerically).
                for k = 1:min(nf, maxxhist)
                    debug_obj.assert(all(xhist(:, k) >= xl, 'all') && all(xhist(:, k) <= xu, 'all'), "XL <= XHIST <= XU", srname);
                end
                debug_obj.assert(size(xpt, 1) == n && size(xpt, 2) == npt, "SIZE(XPT) == [N, NPT]", srname);
                debug_obj.assert(all(infnan_obj.is_finite(xpt), 'all'), "XPT is finite", srname);
                debug_obj.assert(all(xpt >= reshape(sl, [], 1), 'all') && all(xpt <= reshape(su, [], 1), 'all'), "SL <= XPT <= SU", srname);
                debug_obj.assert(size(bmat, 1) == n && size(bmat, 2) == npt + n, "SIZE(BMAT) == [N, NPT+N]", srname);
                debug_obj.assert(linalg_obj.issymmetric(bmat(:, npt + 1:npt + n)), "BMAT(:, NPT+1:NPT+N) is symmetric", srname);
                debug_obj.assert(size(zmat, 1) == npt && size(zmat, 2) == npt - n - 1, "SIZE(ZMAT) == [NPT, NPT-N-1]", srname);

                for j = 1:npt
                    hcol(1:npt) = linalg_obj.matprod21(zmat, zmat(j, :));
                    hcol(npt + 1:npt + n) = bmat(:, j);
                    debug_obj.assert(floor(-log10(eps(class(0.0)))) < floor(-log10(eps(class(0.0)))) || sum(abs(hcol), 'all') > 0, "Column " + string_obj.int2str(j) + " of H is nonzero", srname);
                end
            end

        end
        function [bmat, zmat, info] = updateh_rsc(~, knew, beta, vlag_in, bmat, zmat, varargin)
            % %%% N.B.: UPDATEH_RSC is only used by RESCUE.
            %--------------------------------------------------------------------------------------------------%
            % This subroutine updates arrays BMAT and ZMAT in order to replace the interpolation point
            % XPT(:, KNEW) by XNEW = XPT(:, KOPT) + D. See Section 4 of the BOBYQA paper. [BMAT, ZMAT] describes
            % the matrix H in the BOBYQA paper (eq. 2.7), which is the inverse of the coefficient matrix of the
            % KKT system for the least-Frobenius norm interpolation problem: ZMAT holds a factorization of the
            % leading NPT*NPT submatrix OMEGA of H, the factorization being OMEGA = ZMAT*ZMAT^T; BMAT holds the
            % last N ROWs of H except for the (NPT+1)th column. Note that the (NPT + 1)th row and (NPT + 1)th
            % column of H are not stored as they are unnecessary for the calculation.
            %--------------------------------------------------------------------------------------------------%

            % Common modules
            consts_obj = prima_mat.common.consts_mod();
            debug_obj = prima_mat.common.debug_mod();
            infnan_obj = prima_mat.common.infnan_mod();
            infos_obj = prima_mat.common.infos_mod();
            linalg_obj = prima_mat.common.linalg_mod();
            string_obj = prima_mat.common.string_mod();

            % Inputs


            % VLAG(NPT + N)

            % In-outputs
            % BMAT(N, NPT + N)
            % ZMAT(NPT, NPT-N-1)

            % Outputs


            % Local variables
            srname = "UPDATEH_RSC";


            hcol = NaN(size(bmat, 2), 1);


            v1 = NaN(size(bmat, 1), 1);
            v2 = NaN(size(bmat, 1), 1);
            vlag = NaN(numel(vlag_in), 1);


            % Sizes.
            n = size(bmat, 1);
            npt = size(bmat, 2) - size(bmat, 1);

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(n >= 1, "N >= 1", srname);
                debug_obj.assert(npt >= n + 2, "NPT >= N+2", srname);
                debug_obj.assert(knew >= 1 && knew <= npt, "1 <= KNEW <= NPT", srname);
                debug_obj.assert(size(bmat, 1) == n && size(bmat, 2) == npt + n, "SIZE(BMAT)==[N, NPT+N]", srname);
                debug_obj.assert(linalg_obj.issymmetric(bmat(:, npt + 1:npt + n)), "BMAT(:, NPT+1:NPT+N) is symmetric", srname);
                debug_obj.assert(size(zmat, 1) == npt && size(zmat, 2) == npt - n - 1, "SIZE(ZMAT) == [NPT, NPT-N-1]", srname);
                debug_obj.assert(numel(vlag_in) == npt + n, "SIZE(VLAG) == NPT + N", srname);

                for j = 1:npt
                    hcol(1:npt) = linalg_obj.matprod21(zmat, zmat(j, :));
                    hcol(npt + 1:npt + n) = bmat(:, j);
                    debug_obj.assert(floor(-log10(eps(class(0.0)))) < floor(-log10(eps(class(0.0)))) || sum(abs(hcol), 'all') > 0, "Column " + string_obj.int2str(j) + " of H is nonzero", srname);
                end

                % The following is too expensive to check.
                %tol = 1.0E-2_RP
                %call wassert(errh(bmat, zmat, xpt) <= tol .or. precision(0.0_RP) < precision(0.0D0), &
                %    & 'H = W^{-1} in (2.7) of the BOBYQA paper', srname)

            end

            %====================%
            % Calculation starts %
            %====================%

            ipObj = inputParser();
            addParameter(ipObj, 'info', NaN);
            parse(ipObj, varargin{:});
            info = ipObj.Results.info;
            if nargout >= 3
                info = infos_obj.INFO_DFT;
            end

            % We must not do anything if KNEW is 0. This can only happen sometimes after a trust-region step.
            if knew <= 0
                % KNEW < 0 is impossible if the input is correct.
                return
            end

            % Read VLAG, and calculate parameters for the updating formula (4.9) and (4.14) of the BOBYQA paper.
            vlag(:) = vlag_in;
            tau = vlag(knew);
            % In theory, DENOM can also be calculated after ZMAT is rotated below. However, this worsened the
            % performance of BOBYQA in a test on 20220413.
            denom = sum(zmat(knew, :) .^ 2, 'all') * beta + tau ^ 2;

            % Quite rarely, due to rounding errors, VLAG or BETA may not be finite, or DENOM may not be
            % positive. In such cases, [BMAT, ZMAT] would be destroyed by the update, and hence we would rather
            % not update them at all. Or should we simply terminate the algorithm?
            if ~(infnan_obj.is_finite(sum(abs(vlag), 'all') + abs(beta)) && denom > 0)
                if nargout >= 3
                    info = infos_obj.DAMAGING_ROUNDING;
                end
                return
            end

            % After the following line, VLAG = H*w - e_KNEW in the NEWUOA paper (where t = KNEW).
            vlag(knew) = vlag(knew) - consts_obj.ONE;

            % Apply Givens rotations to put zeros in the KNEW-th row of ZMAT. After this, ZMAT(KNEW, :) contains
            % only one nonzero at ZMAT(KNEW, 1). Entries of ZMAT are treated as 0 if the moduli are quite small.
            for j = 2:npt - n - 1
                if abs(zmat(knew, j)) > 1.0e-20 * max(abs(zmat), [], 'all')
                    % This threshold is by Powell
                    grot = linalg_obj.planerot(zmat(knew, [1, j]));
                    zmat(:, [1, j]) = linalg_obj.matprod22(zmat(:, [1, j]), grot.');
                end
                zmat(knew, j) = consts_obj.ZERO;
            end

            % Put the KNEW-th column of the unupdated H (except for the (NPT+1)th entry) into HCOL.
            hcol(1:npt) = zmat(knew, 1) * zmat(:, 1);
            hcol(npt + 1:npt + n) = bmat(:, knew);

            % Complete the updating of ZMAT. See (4.14) of the BOBYQA paper.
            sqrtdn = sqrt(denom);
            zknew1 = zmat(knew, 1) / sqrtdn;
            zmat(:, 1) = (tau / sqrtdn) * zmat(:, 1) - zknew1 * vlag(1:npt);
            zmat(knew, 1) = zknew1; % Because TAU = VLAG(KNEW) + 1. Powell's code does not have this.

            % Finally, update the matrix BMAT. It implements the last N rows of (4.9) in the BOBYQA paper.
            alpha = hcol(knew);
            v1(:) = (alpha * vlag(npt + 1:npt + n) - tau * hcol(npt + 1:npt + n)) ./ denom;
            v2(:) = (-beta * hcol(npt + 1:npt + n) - tau * vlag(npt + 1:npt + n)) ./ denom;
            bmat(:, :) = bmat + linalg_obj.outprod(v1, vlag) + linalg_obj.outprod(v2, hcol); %call r2update(bmat, ONE, v1, vlag, ONE, v2, hcol)
            % N.B.: The use of OUTPROD is expensive memory-wise, but it is not our concern in this implementation.
            % Numerically, the update above does not guarantee BMAT(:, NPT+1 : NPT+N) to be symmetric.
            A_slice = linalg_obj.symmetrize(bmat(:, npt + 1:npt + n)); bmat(:, npt + 1:npt + n) = A_slice;

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(size(bmat, 1) == n && size(bmat, 2) == npt + n, "SIZE(BMAT)==[N, NPT+N]", srname);
                debug_obj.assert(linalg_obj.issymmetric(bmat(:, npt + 1:npt + n)), "BMAT(:, NPT+1:NPT+N) is symmetric", srname);
                debug_obj.assert(size(zmat, 1) == npt && size(zmat, 2) == npt - n - 1, "SIZE(ZMAT) == [NPT, NPT-N-1]", srname);

                for j = 1:npt
                    hcol(1:npt) = linalg_obj.matprod21(zmat, zmat(j, :));
                    hcol(npt + 1:npt + n) = bmat(:, j);
                    debug_obj.assert(floor(-log10(eps(class(0.0)))) < floor(-log10(eps(class(0.0)))) || sum(abs(hcol), 'all') > 0, "Column " + string_obj.int2str(j) + " of H is nonzero", srname);
                end

                % The following is too expensive to check.
                % %if (n * npt <= 50) then
                % %    xpt_test = xpt
                % %    xpt_test(:, knew) = xpt(:, kopt) + d
                % %    call assert(errh(bmat, zmat, xpt_test) <= tol .or. precision(0.0_RP) < precision(0.0D0), &
                % %        & 'H = W^{-1} in (2.7) of the BOBYQA paper', srname)
                % %end if

            end
        end

    end
end
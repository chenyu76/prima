classdef trustregion_cobyla_mod
    %--------------------------------------------------------------------------------------------------%
    % This module provides subroutines concerning the trust-region calculations of COBYLA.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code and the COBYLA paper.
    %
    % Dedicated to the late Professor M. J. D. Powell FRS (1936--2015).
    %
    % Started: June 2021
    %
    % Last Modified: Saturday, March 16, 2024 AM03:37:33
    %--------------------------------------------------------------------------------------------------%

    methods
        function d = trstlp(obj, A, b, delta, g)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine calculates an N-component vector D by the following two stages. In the first
            % stage, D is set to the shortest vector that minimizes the greatest violation of the constraints
            %       A^T * D <= B,  K = 1, 2, 3, ..., M,
            % subject to the Euclidean length of D being at most DELTA. If its length is strictly less than
            % DELTA, then the second stage uses the resultant freedom in D to minimize the objective function
            %       G^T * D
            % subject to no increase in any greatest constraint violation.
            %
            % It is possible but rare that a degeneracy may prevent D from attaining the target length DELTA.
            %
            % CVIOL is the largest constraint violation of the current D: MAXVAL([A^T*D-B, ZERO]).
            % ICON is the index of a most violated constraint if CVIOL is positive.
            %
            % NACT is the number of constraints in the active set and IACT(1), ...,IACT(NACT) are their indices,
            % while the remainder of IACT contains a permutation of the remaining constraint indices.
            % N.B.: NACT <= min(M, N). Obviously, NACT <= M. In addition, The constraints in IACT(1, ..., NACT)
            % have linearly independent gradients (see the comments above the instructions that delete a
            % constraint from the active set to make room for the new active constraint with index IACT(ICON));
            % it can also be seen from the update of NACT: starting from 0, NACT is incremented only if NACT < N.
            %
            % Further, Z is an orthogonal matrix whose first NACT columns can be regarded as the result of
            % Gram-Schmidt applied to the active constraint gradients. For J = 1, 2, ..., NACT, the number
            % ZDOTA(J) is the scalar product of the J-th column of Z with the gradient of the J-th active
            % constraint. D is the current vector of variables and here the residuals of the active constraints
            % should be zero. Further, the active constraints have nonnegative Lagrange multipliers that are
            % held at the beginning of VMULTC. The remainder of this vector holds the residuals of the inactive
            % constraints at D, the ordering of the components of VMULTC being in agreement with the permutation
            % of the indices of the constraints that is in IACT. All these residuals are nonnegative, which is
            % achieved by the shift CVIOL that makes the least residual zero.
            %
            % N.B.:
            % 0. In Powell's implementation, the constraints are A^T * D >= B. In other words, the A and B in
            % our implementation are the negative of those in Powell's implementation.
            % 1. The algorithm was NOT documented in the COBYLA paper. A note should be written to introduce it!
            % 2. As a major part of the algorithm (see TRSTLP_SUB), the code maintains and updates the QR
            % factorization of A(:, IACT(1:NACT)), i.e., the gradients of all the active (linear) constraints.
            % The matrix Z is indeed Q, and the vector ZDOTA is the diagonal of R. The factorization is updated
            % by Givens rotations when an index is added in or removed from IACT.
            % 3. There are probably better algorithms available for the trust-region linear programming problem.
            %--------------------------------------------------------------------------------------------------%

            % Common modules
            consts_obj = prima_mat.common.consts_mod();
            debug_obj = prima_mat.common.debug_mod();
            infnan_obj = prima_mat.common.infnan_mod();
            linalg_obj = prima_mat.common.linalg_mod();

            % Inputs
            % A(N, M)
            % B(M)

            % G(N)

            % Outputs
            d = NaN(size(A, 1), 1); % D(N)

            % Local variables
            srname = "TRSTLP";

            iact = NaN(numel(b) + 1, 1);


            nact = NaN;
            A_aug = NaN(size(A, 1), size(A, 2) + 1);
            b_aug = NaN(numel(b) + 1, 1);

            vmultc = NaN(numel(b) + 1, 1);
            z = NaN(numel(d));

            % Sizes
            m = size(A, 2);
            n = size(A, 1);

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(n >= 1, "N >= 1", srname);
                debug_obj.assert(m >= 0, "M >= 0", srname);
                debug_obj.assert(numel(g) == n, "SIZE(G) == N", srname);
                debug_obj.assert(numel(d) == n, "SIZE(D) == N", srname);
                debug_obj.assert(numel(b) == m, "SIZE(B) == M", srname);
                debug_obj.assert(delta > 0, "DELTA > 0", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            % Form A_aug and B_aug. This allows the gradient of the objective function to be regarded as the
            % gradient of a constraint in the second stage.
            A_aug(:, :) = reshape([reshape(A, 1, []), reshape(g, 1, [])], [n, m + 1]); %%MATLAB: A_aug = [A, g];
            b_aug(:) = [reshape(b, [], 1); consts_obj.ZERO]; %%MATLAB: b_aug = [b; 0];

            % Scale the problem if A_aug contains large values. Otherwise, floating point exceptions may occur.
            % Note that the trust-region step is scale invariant.
            % N.B.: It is faster and safer to scale by multiplying a reciprocal than by division. See
            % https://fortran-lang.discourse.group/t/ifort-ifort-2021-8-0-1-0e-37-1-0e-38-0/
            for i = 1:m + 1                % Note that SIZE(A, 2) = SIZE(B) = M + 1 /= M.
                if max(abs(A_aug(:, i)), [], 'all') > 1.0e12
                    modscal = max(consts_obj.TWO * consts_obj.REALMIN, consts_obj.ONE / max(abs(A_aug(:, i)), [], 'all')); % MAX: avoid underflow.
                    A_aug(:, i) = A_aug(:, i) * modscal;
                    b_aug(i) = b_aug(i) * modscal;
                end
            end

            % Stage 1: minimize the l_infinity constraint violation of the linearized constraints.
            [iact_slice, nact, d, vmultc_slice, z] = obj.trstlp_sub(iact(1:m), nact, 1, A_aug(:, 1:m), b_aug(1:m), delta, d, vmultc(1:m), z); iact(1:m) = iact_slice; vmultc(1:m) = vmultc_slice;

            % Stage 2: minimize the linearized objective without increasing the l_infinity constraint violation.
            [iact, nact, d, vmultc, z] = obj.trstlp_sub(iact, nact, 2, A_aug, b_aug, delta, d, vmultc, z);

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(d) == n, "SIZE(D) == N", srname);
                debug_obj.assert(all(infnan_obj.is_finite(d), 'all'), "D is finite", srname);
                % Due to rounding, it may happen that ||D|| > DELTA, but ||D|| > 2*DELTA is highly improbable.
                debug_obj.assert(linalg_obj.p_norm(d) <= consts_obj.TWO * delta, "||D|| <= 2*DELTA", srname);
            end
        end
        function [iact, nact, d, vmultc, z] = trstlp_sub(~, iact, nact, stage, A, b, delta, d, vmultc, z)
            %--------------------------------------------------------------------------------------------------%
            % This subroutine does the real calculations for TRSTLP, both stage 1 and stage 2.
            % Major differences between stage 1 and stage 2:
            % 1. Initialization. Stage 2 inherits the values of some variables from stage 1, so they are
            % initialized in stage 1 but not in stage 2.
            % 2. CVIOL. CVIOL is updated after at iteration in stage 1, while it remains a constant in stage 2.
            % 3. SDIRN. See the definition of SDIRN in the code for details.
            % 4. OPTNEW. The two stages have different objectives, so OPTNEW is updated differently.
            % 5. STEP. STEP <= CVIOL in stage 1.
            %--------------------------------------------------------------------------------------------------%

            % Common modules
            consts_obj = prima_mat.common.consts_mod();
            debug_obj = prima_mat.common.debug_mod();
            infnan_obj = prima_mat.common.infnan_mod();
            linalg_obj = prima_mat.common.linalg_mod();
            powalg_obj = prima_mat.common.powalg_mod();

            % Inputs

            % A(N, MCON)
            % B(M)


            % In-outputs
            % IACT(MCON)

            % D(N)
            % VMULTC(MCON)
            % Z(N, N)

            % Local variables
            srname = "TRSTLP_SUB";


            nactsav = NaN;


            %real(RP) :: cvold
            cvsabs = NaN(numel(b), 1);
            cvshift = NaN(numel(b), 1);
            dd = NaN;
            dnew = NaN(numel(d), 1);
            dold = NaN(numel(d), 1);
            frac = NaN;
            fracmult = NaN(numel(vmultc), 1);


            sd = NaN;
            sdirn = NaN(numel(d), 1);
            sqrtd = NaN;
            ss = NaN;
            step = NaN;
            vmultd = NaN(numel(vmultc), 1);
            zdasav = NaN(size(z, 2), 1);
            zdota = NaN(size(z, 2), 1);

            % Sizes
            mcon = size(A, 2);
            n = size(A, 1);

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(n >= 1, "N >= 1", srname);
                debug_obj.assert(stage == 1 || stage == 2, "STAGE == 1 or 2", srname);
                debug_obj.assert((mcon >= 0 && stage == 1) || (mcon >= 1 && stage == 2), "MCON >= 1 in stage 1 and MCON >= 0 in stage 2", srname);
                debug_obj.assert(numel(b) == mcon, "SIZE(B) == MCON", srname);
                debug_obj.assert(numel(iact) == mcon, "SIZE(IACT) == MCON", srname);
                debug_obj.assert(numel(vmultc) == mcon, "SIZE(VMULTC) == MCON", srname);
                debug_obj.assert(numel(d) == n, "SIZE(D) == N", srname);
                debug_obj.assert(size(z, 1) == n && size(z, 2) == n, "SIZE(Z) == [N, N]", srname);
                debug_obj.assert(delta > 0, "DELTA > 0", srname);
                if stage == 2
                    debug_obj.assert(all(infnan_obj.is_finite(d), 'all') && linalg_obj.p_norm(d) <= consts_obj.TWO * delta, "D is finite and ||D|| <= 2*DELTA at the beginning of stage 2", srname);
                    debug_obj.assert((nact >= 0 && nact <= min(mcon, n)), "0 <= NACT <= MIN(MCON, N) at the beginning of stage 2", srname);
                    debug_obj.assert(floor(-log10(eps(class(0.0)))) < floor(-log10(eps(class(0.0)))) || all(vmultc(1:mcon - 1) >= 0, 'all'), "VMULTC >= 0 at the beginning of stage 2", srname);
                    % N.B.: Stage 1 defines only VMULTC(1:M); VMULTC(M+1) is undefined!

                end
            end

            %====================%
            % Calculation starts %
            %====================%

            % Initialization according to STAGE.
            if stage == 1
                iact(:) = linalg_obj.linspace_i(1, mcon, mcon); %%MATLAB: iact = (1:mcon);  % Row vector
                % N.B.: 1. The MATLAB version of LINSPACE returns a row vector. Take a transpose if needed.
                % 2. In MATLAB, linspace(1, mcon, mcon) can also be written as (1:mcon).
                nact = 0;
                d(:) = consts_obj.ZERO;
                cviol = linalg_obj.maximum1([consts_obj.ZERO; reshape(-b, [], 1)]);
                vmultc(:) = cviol + b;
                z(:, :) = linalg_obj.eye1(n);
                if mcon == 0 || cviol <= 0
                    % Check whether a quick return is possible. Make sure the In-outputs have been initialized.
                    return
                end

                if all(infnan_obj.is_nan_sp(b), 'all')
                    return
                else
                    icon = fix(fortran.maxloc(-b, 'mask', (~infnan_obj.is_nan_sp(b)), 'dim', 1));
                    %%MATLAB: [~, icon] = max(b, [], 'omitnan');
                end
                m = mcon;
                sdirn(:) = consts_obj.ZERO;
            else
                if linalg_obj.inprod(d, d) >= delta ^ 2
                    % Check whether a quick return is possible.
                    return
                end

                iact(mcon) = mcon;
                vmultc(mcon) = consts_obj.ZERO;
                m = mcon - 1;
                icon = mcon;

                % In Powell's code, stage 2 uses the ZDOTA and CVIOL calculated by stage 1. Here we re-calculate
                % them so that they need not be passed from stage 1 to 2, and hence the coupling is reduced.
                cviol = linalg_obj.maximum1([consts_obj.ZERO; reshape(linalg_obj.matprod12(d, A(:, 1:m)) - b(1:m), [], 1)]);
            end
            zdota(1:nact) = reshape(arrayfun(@(k) linalg_obj.inprod(z(:, k), A(:, iact(k))), 1:nact), [], 1);
            %%MATLAB: zdota(1:nact) = sum(z(:, 1:nact) .* A(:, iact(1:nact)), 1);  % Row vector

            % More initialization.
            optold = consts_obj.REALMAX;
            nactold = nact;
            nfail = 0;

            %----------------------------------------------------------------------------------------------%
            % Zaikun 20211011: VMULTD is computed from scratch at each iteration, but VMULTC is inherited.
            %----------------------------------------------------------------------------------------------%

            % Powell's code can encounter infinite cycling, which did happen when testing the following CUTEst
            % problems: DANWOODLS, GAUSS1LS, GAUSS2LS, GAUSS3LS, KOEBHELB, TAX13322, TAXR13322. Indeed, in all
            % these cases, Inf/NaN appear in D due to extremely large values in A (up to 10^219). To resolve
            % this, we set the maximal number of iterations to MAXITER, and terminate if Inf/NaN occurs in D.
            % The formulation of MAXITER below contains a precaution against overflow. In MATLAB/Python/Julia/R,
            % we can write maxiter = min(10000, 100*max(m, n))
            maxiter = fix(min(10 ^ min(4, floor(log10(double(intmax('int64'))))), 100 * fix(max(m, n))));
            for iter = 1:maxiter
                if consts_obj.DEBUGGING
                    debug_obj.assert(floor(-log10(eps(class(0.0)))) < floor(-log10(eps(class(0.0)))) || all(vmultc >= 0, 'all'), "VMULTC >= 0", srname);
                end
                if stage == 1
                    optnew = cviol;
                else
                    optnew = linalg_obj.inprod(d, A(:, mcon));
                end

                % End the current stage of the calculation if 3 consecutive iterations have either failed to
                % reduce the best calculated value of the objective function or to increase the number of active
                % constraints since the best value was calculated. This strategy prevents cycling, but there is
                % a remote possibility that it will cause premature termination.
                if optnew < optold || nact > nactold
                    nactold = nact;
                    nfail = 0;
                else
                    nfail = nfail + 1;
                end
                optold = min(optold, optnew);
                if nfail == 3
                    break
                end

                % If ICON exceeds NACT, then we add the constraint with index IACT(ICON) to the active set.
                if icon > nact
                    zdasav(1:nact) = zdota(1:nact);
                    nactsav = nact;
                    [z, zdota, nact] = powalg_obj.qradd_Rdiag(A(:, iact(icon)), z, zdota, nact); % QRADD may update NACT to NACT + 1.
                    % Indeed, it suffices to pass ZDOTA(1:MIN(N, NACT+1)) to QRADD as follows.
                    % %call qradd(A(:, iact(icon)), z, zdota(1:min(n, nact + 1_IK)), nact)

                    if nact == nactsav + 1
                        % N.B.: It is problematic to index arrays using [NACT, ICON] when NACT == ICON.
                        % Zaikun 20211012: Why should VMULTC(NACT) = 0?
                        if nact ~= icon
                            vmultc([icon, nact]) = [vmultc(nact), consts_obj.ZERO];
                            iact([icon, nact]) = iact([nact, icon]);
                        else
                            vmultc(nact) = consts_obj.ZERO;
                        end
                    else
                        % Zaikun 20211011:
                        % 1. VMULTD is calculated from scratch for the first time (out of 2) in one iteration.
                        % 2. NOTE that IACT has not been updated to replace IACT(NACT) with IACT(ICON). Thus
                        % A(:, IACT(1:NACT)) is the UNUPDATED version before QRADD (Z(:, 1:NACT) remains the
                        % same before and after QRADD). Therefore, if we supply ZDOTA to LSQR (as Rdiag) as
                        % Powell did, we should use the UNUPDATED version, namely ZDASAV.
                        vmultd(1:nact) = linalg_obj.lsqr_Rdiag(A(:, iact(1:nact)), A(:, iact(icon)), 'Q', z(:, 1:nact), 'Rdiag', zdasav(1:nact));
                        if ~any(vmultd(1:nact) > 0 & iact(1:nact) <= m, 'all')
                            % N.B.: This can be triggered by NACT == 0 (among other possibilities)! This is
                            % important, because NACT will be used as an index in the sequel.
                            break
                        end
                        % VMULTD(NACT+1:MCON) is not used, but we have to initialize it in Fortran, or compilers
                        % complain about the WHERE construct below (another solution: restrict WHERE to 1:NACT).
                        vmultd(nact + 1:mcon) = -consts_obj.ONE; % SIZE(VMULTD) = MCON

                        % Revise the Lagrange multipliers. The revision is not applicable to VMULTC(NACT + 1:M).
                        fracmult(:) = consts_obj.REALMAX;
                        fracmult(vmultd > 0 & iact <= m) = vmultc(vmultd > 0 & iact <= m) ./ vmultd(vmultd > 0 & iact <= m);
                        %%MATLAB: mask = (vmultd > 0 & iact <= m); fracmult(mask) = vmultc(mask) / vmultd(mask);
                        % Only the places with VMULTD > 0 and IACT <= M is relevant blow, if any.
                        frac = min(fracmult(1:nact), [], 'all'); % FRACMULT(NACT+1:MCON) may contain garbage.
                        vmultc(1:nact) = max(consts_obj.ZERO, vmultc(1:nact) - frac * vmultd(1:nact));

                        % Reorder the active constraints so that the one to be replaced is at the end of the list.
                        % Exit if the new value of ZDOTA(NACT) is not acceptable. Powell's condition for the
                        % following IF: .NOT. ABS(ZDOTA(NACT)) > 0. Note that it is different from
                        % 'ABS(ZDOTA(NACT) <= 0)', as ZDOTA(NACT) can be NaN.
                        % N.B.: We cannot arrive here with NACT == 0, which should have triggered an exit above.
                        if infnan_obj.is_nan_sp(zdota(nact)) || abs(zdota(nact)) <= consts_obj.EPS ^ 2
                            break
                        end
                        vmultc([icon, nact]) = [consts_obj.ZERO, frac]; % VMULTC([ICON, NACT]) is valid as ICON > NACT.
                        iact([icon, nact]) = iact([nact, icon]);
                    end

                    % In stage 2, ensure that the objective continues to be treated as the last active constraint.
                    % Zaikun 20211011, 20211111: Is it guaranteed for stage 2 that IACT(NACT-1) = MCON when
                    % IACT(NACT) /= MCON??? If not, then how does the following procedure ensure that MCON is
                    % the last of IACT(1:NACT)?
                    if stage == 2 && iact(nact) ~= mcon
                        if nact <= 1
                            % We must exit, as NACT-1 is used as an index below. Powell's code does not have this.
                            break
                        end
                        [z, Rdiag_slice] = powalg_obj.qrexc_Rdiag(A(:, iact(1:nact)), z, zdota(1:nact), nact - 1); zdota(1:nact) = Rdiag_slice;
                        % Indeed, it suffices to pass Z(:, 1:NACT) to QREXC as follows.
                        % %call qrexc(A(:, iact(1:nact)), z(:, 1:nact), zdota(1:nact), nact - 1_IK)
                        iact([nact - 1, nact]) = iact([nact, nact - 1]);
                        vmultc([nact - 1, nact]) = vmultc([nact, nact - 1]);
                    end
                    % Zaikun 20211117: It turns out that the last few lines do not guarantee IACT(NACT) == N in
                    % stage 2; the following test cannot be passed. IS THIS A BUG?!
                    % %call assert(iact(nact) == mcon .or. stage == 1, 'IACT(NACT) == MCON in stage 2', srname)

                    % Powell's code does not have the following. It avoids subsequent floating point exceptions.
                    %------------------------------------------------------------------------------------------%
                    if infnan_obj.is_nan_sp(zdota(nact)) || abs(zdota(nact)) <= consts_obj.EPS ^ 2
                        break
                    end
                    %------------------------------------------------------------------------------------------%

                    % Set SDIRN to the direction of the next change to the current vector of variables.
                    % Usually during stage 1 the vector SDIRN gives a search direction that reduces all the
                    % active constraint violations by one simultaneously.
                    if stage == 1
                        sdirn(:) = sdirn - ((linalg_obj.inprod(sdirn, A(:, iact(nact))) + consts_obj.ONE) / zdota(nact)) * z(:, nact);
                    else
                        sdirn(:) = -(consts_obj.ONE / zdota(nact)) * z(:, nact);
                        % SDIRN = Z(:, NACT)/(A(:,IACT(NACT))^T*Z(:, NACT))
                        % SDIRN^T*A(:, IACT(NACT)) = 1, SDIRN is orthogonal to A(:, IACT(1:NACT-1)) and is
                        % parallel to Z(:, NACT).
                    end
                else                    % ICON <= NACT
                    % Delete the constraint with the index IACT(ICON) from the active set, which is done by
                    % reordering IACT(ICONT:NACT) into [IACT(ICON+1:NACT), IACT(ICON)] by pairwise exchanges
                    % and then reduce NACT to NACT - 1. In theory, ICON > 0.
                    debug_obj.validate(icon > 0, "ICON > 0", srname);
                    [z, Rdiag_slice] = powalg_obj.qrexc_Rdiag(A(:, iact(1:nact)), z, zdota(1:nact), icon); zdota(1:nact) = Rdiag_slice; % QREXC does nothing if ICON==NACT.
                    % Indeed, it suffices to pass Z(:, 1:NACT) to QREXC as follows.
                    % %call qrexc(A(:, iact(1:nact)), z(:, 1:nact), zdota(1:nact), icon)
                    iact(icon:nact) = [reshape(iact(icon + 1:nact), [], 1); iact(icon)];
                    vmultc(icon:nact) = [reshape(vmultc(icon + 1:nact), [], 1); vmultc(icon)];
                    nact = nact - 1;

                    % Powell's code does not have the following. It avoids subsequent exceptions.
                    %------------------------------------------------------------------------------------------%
                    % Zaikun 20221212: In theory, NACT > 0 in stage 2, as the objective function should always
                    % be considered as an "active constraint" --- more precisely, IACT(NACT) = MCON. However,
                    % looking at the code, I cannot see why in stage 2 NACT must be positive after the reduction
                    % above. It did happen in stage 1 that NACT became 0 after the reduction --- this is
                    % extremely rare, and it was never observed until 20221212, after almost one year of
                    % random tests. Maybe NACT is theoretically positive even in stage 1?
                    if stage == 2 && nact <= 0
                        break % If this case ever occurs, we have to exit, as NACT is used as an index below.

                    end
                    if nact > 0
                        if infnan_obj.is_nan_sp(zdota(nact)) || abs(zdota(nact)) <= consts_obj.EPS ^ 2
                            break
                        end
                    end
                    %------------------------------------------------------------------------------------------%

                    % Set SDIRN to the direction of the next change to the current vector of variables.
                    if stage == 1
                        sdirn(:) = sdirn - linalg_obj.inprod(sdirn, z(:, nact + 1)) * z(:, nact + 1);
                        % SDIRN is orthogonal to Z(:, NACT+1)

                    else
                        sdirn(:) = -(consts_obj.ONE / zdota(nact)) * z(:, nact);
                        % SDIRN = Z(:, NACT)/(A(:,IACT(NACT))^T*Z(:, NACT))
                        % SDIRN^T*A(:, IACT(NACT)) = 1, SDIRN is orthogonal to A(:, IACT(1:NACT-1)) and is
                        % parallel to Z(:, NACT).
                    end
                end

                % Calculate the step to the trust region boundary or take the step that reduces CVIOL to 0.
                %----------------------------------------------------------------------------------------------%
                % The following calculation of STEP is adopted from NEWUOA/BOBYQA/LINCOA. It seems to improve
                % the performance of COBYLA. We also found that removing the precaution about underflows is
                % beneficial to the overall performance of COBYLA --- the underflows are harmless anyway.
                dd = delta ^ 2 - linalg_obj.inprod(d, d);
                ss = linalg_obj.inprod(sdirn, sdirn);
                sd = linalg_obj.inprod(sdirn, d);
                if dd <= 0 || ss <= consts_obj.EPS * delta ^ 2 || infnan_obj.is_nan_sp(sd)
                    break
                end
                % SQRTD: square root of a discriminant. The MAXVAL avoids SQRTD < ABS(SD) due to underflow.
                sqrtd = max([sqrt(ss * dd + sd ^ 2), abs(sd), sqrt(ss * dd)], [], 'all');
                if sd > 0
                    step = dd / (sqrtd + sd);
                else
                    step = (sqrtd - sd) / ss;
                end
                % STEP < 0 should not happen. STEP can be 0 or NaN when, e.g., SD or SS becomes Inf.
                if step <= 0 || ~infnan_obj.is_finite(step)
                    break
                end
                % Powell's approach and comments are as follows.
                %----------------------------------------------------------------%
                % The two statements below that include the factor EPS prevent
                % some harmless underflows that occurred in a test calculation
                % (Zaikun: here, EPS is the machine epsilon; Powell's original
                % code used 1.0E-6, and Powell's code was written in SINGLE
                % PRECISION). Further, we skip the step if it could be zero within
                % a reasonable tolerance for computer rounding errors.
                %
                % %dd = delta**2 - sum(d**2, mask=(abs(d) >= EPS * delta))
                % %ss = inprod(sdirn, sdirn)
                % %if (dd  <= 0) then
                % %    exit
                % %end if
                % %sd = inprod(sdirn, d)
                % %if (abs(sd) >= EPS * sqrt(ss * dd)) then
                % %    step = dd / (sqrt(ss * dd + sd**2) + sd)
                % %else
                % %    step = dd / (sqrt(ss * dd) + sd)
                % %end if
                %----------------------------------------------------------------%
                %----------------------------------------------------------------------------------------------%

                if stage == 1
                    if linalg_obj.isminor0(cviol, step)
                        break
                    end
                    step = min(step, cviol);
                end

                % Set DNEW to the new variables if STEP is the steplength, and reduce CVIOL to the corresponding
                % maximum residual if stage 1 is being done.
                dnew(:) = d + step * sdirn;
                if stage == 1
                    %cvold = cviol
                    cviol = linalg_obj.maximum1([consts_obj.ZERO; reshape(linalg_obj.matprod12(dnew, A(:, iact(1:nact))) - b(iact(1:nact)), [], 1)]);
                    % N.B.: CVIOL will be used when calculating VMULTD(NACT+1 : MCON).

                end

                % Zaikun 20211011:
                % 1. VMULTD is computed from scratch for the second (out of 2) time in one iteration.
                % 2. VMULTD(1:NACT) and VMULTD(NACT+1:MCON) are calculated separately with no coupling.
                % 3. VMULTD will be calculated from scratch again in the next iteration.
                % Set VMULTD to the VMULTC vector that would occur if D became DNEW. A device is included to
                % force VMULTD(K)=ZERO if deviations from this value can be attributed to computer rounding
                % errors. First calculate the new Lagrange multipliers.
                vmultd(1:nact) = -linalg_obj.lsqr_Rdiag(A(:, iact(1:nact)), dnew, 'Q', z(:, 1:nact), 'Rdiag', zdota(1:nact));
                if stage == 2
                    vmultd(nact) = max(consts_obj.ZERO, vmultd(nact)); % This seems never activated.

                end
                % Complete VMULTD by finding the new constraint residuals. (Powell wrote "Complete VMULTC ...")
                cvshift(:) = cviol - (linalg_obj.matprod12(dnew, A(:, iact)) - b(iact)); % Only CVSHIFT(nact+1:mcon) is needed.
                cvsabs(:) = linalg_obj.matprod12(abs(dnew), abs(A(:, iact))) + abs(b(iact)) + cviol;
                cvshift(linalg_obj.trueloc(linalg_obj.isminor1(cvshift, cvsabs))) = consts_obj.ZERO;
                %%MATLAB: cvshift(isminor(cvshift, cvsabs)) = 0;
                vmultd(nact + 1:mcon) = cvshift(nact + 1:mcon);

                % Calculate the fraction of the step from D to DNEW that will be taken.
                fracmult(:) = consts_obj.REALMAX;
                fracmult(vmultd < 0) = vmultc(vmultd < 0) ./ (vmultc(vmultd < 0) - vmultd(vmultd < 0));
                %%MATLAB: mask = (vmultd < 0); fracmult(mask) = vmultc(mask) / (vmultc(mask) - vmultd(mask));
                % Only the places with VMULTD < 0 is relevant below, if any.
                icon = fix(fortran.minloc([consts_obj.ONE; reshape(fracmult, [], 1)], 'dim', 1) - 1);
                frac = min([consts_obj.ONE; reshape(fracmult, [], 1)], [], 'all');
                %%MATLAB: [frac, icon] = min([1, fracmult]); icon = icon - 1

                % Update D, VMULTC and CVIOL.
                dold(:) = d;
                d(:) = (consts_obj.ONE - frac) * d + frac * dnew;
                vmultc(:) = max(consts_obj.ZERO, (consts_obj.ONE - frac) * vmultc + frac * vmultd);
                % Exit in case of Inf/NaN in D or VMULTC.
                if ~(infnan_obj.is_finite(sum(abs(d), 'all')) && infnan_obj.is_finite(sum(abs(vmultc), 'all')))
                    d(:) = dold; % Should we restore also IACT, NACT, VMULTC, and Z?
                    break
                end

                if stage == 1
                    %cviol = (ONE - frac) * cvold + frac * cviol  ! Powell's version
                    % In theory, CVIOL = MAXVAL([MATPROD(D, A) - B, ZERO]), yet the CVIOL updated as above
                    % can be quite different from this value if A has huge entries (e.g., > 1E20).
                    cviol = linalg_obj.maximum1([consts_obj.ZERO; reshape(linalg_obj.matprod12(d, A) - b, [], 1)]);
                end

                if icon < 1 || icon > mcon
                    % In Powell's code, the condition is ICON == 0. Indeed, ICON < 0 cannot hold unless
                    % FRACMULT contains only NaN, which should not happen; ICON > MCON should never occur.
                    break
                end
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(numel(iact) == mcon, "SIZE(IACT) == MCON", srname);
                debug_obj.assert(numel(vmultc) == mcon, "SIZE(VMULTC) == MCON", srname);
                debug_obj.assert(floor(-log10(eps(class(0.0)))) < floor(-log10(eps(class(0.0)))) || all(vmultc >= 0, 'all'), "VMULTC >= 0", srname);
                debug_obj.assert(numel(d) == n, "SIZE(D) == N", srname);
                debug_obj.assert(all(infnan_obj.is_finite(d), 'all'), "D is finite", srname);
                debug_obj.assert(linalg_obj.p_norm(d) <= consts_obj.TWO * delta, "||D|| <= 2*DELTA", srname);
                debug_obj.assert(size(z, 1) == n && size(z, 2) == n, "SIZE(Z) == [N, N]", srname);
                debug_obj.assert(nact >= 0 && nact <= min(mcon, n), "0 <= NACT <= MIN(MCON, N)", srname);
            end

        end
        function delta = trrad(~, delta_in, dnorm, eta1, eta2, gamma1, gamma2, ratio)
            %--------------------------------------------------------------------------------------------------%
            % This function updates the trust region radius according to RATIO and DNORM.
            %--------------------------------------------------------------------------------------------------%

            % Generic module
            consts_obj = prima_mat.common.consts_mod();
            infnan_obj = prima_mat.common.infnan_mod();
            debug_obj = prima_mat.common.debug_mod();


            % Input
            % Current trust-region radius
            % Norm of current trust-region step
            % Ratio threshold for contraction
            % Ratio threshold for expansion
            % Contraction factor
            % Expansion factor
            % Reduction ratio

            % Outputs
            delta = NaN;

            % Local variables
            srname = "TRRAD";

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(delta_in >= dnorm && dnorm > 0, "DELTA_IN >= DNORM > 0", srname);
                debug_obj.assert(eta1 >= 0 && eta1 <= eta2 && eta2 < 1, "0 <= ETA1 <= ETA2 < 1", srname);
                debug_obj.assert(eta1 >= 0 && eta1 <= eta2 && eta2 < 1, "0 <= ETA1 <= ETA2 < 1", srname);
                debug_obj.assert(gamma1 > 0 && gamma1 < 1 && gamma2 > 1, "0 < GAMMA1 < 1 < GAMMA2", srname);
                % By the definition of RATIO in ratio.f90, RATIO cannot be NaN unless the actual reduction is
                % NaN, which should NOT happen due to the moderated extreme barrier.
                debug_obj.assert(~infnan_obj.is_nan_sp(ratio), "RATIO is not NaN", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            if ratio <= eta1
                delta = gamma1 * dnorm; % Powell's UOBYQA/NEWUOA.
                %delta = gamma1 * delta_in  ! Powell's COBYLA/LINCOA.
                %delta = min(gamma1 * delta_in, dnorm)  ! Powell's BOBYQA.

            elseif ratio <= eta2
                delta = max(gamma1 * delta_in, dnorm); % Powell's UOBYQA/NEWUOA/BOBYQA/LINCOA

            else
                delta = max(gamma1 * delta_in, gamma2 * dnorm); % Powell's NEWUOA/BOBYQA.
                %delta = max(delta_in, gamma2 * dnorm)  ! Modified version. Works well for UOBYQA.
                % For noise-free CUTEst problems of <= 100 variables, Powell's version works slightly better
                % than the modified one.
                %delta = max(delta_in, 1.25_RP * dnorm, dnorm + rho)  ! Powell's UOBYQA.
                %delta = min(max(gamma1 * delta_in, gamma2 * dnorm), sqrt(gamma2) * delta_in)  ! Powell's LINCOA.
            end

            % For noisy problems, the following may work better.
            % %if (ratio <= eta1) then
            % %    delta = gamma1 * dnorm
            % %elseif (ratio <= eta2) then  ! Ensure DELTA >= DELTA_IN
            % %    delta = delta_in
            % %else  ! Ensure DELTA > DELTA_IN with a constant factor
            % %    delta = max(delta_in * (1.0_RP + gamma2) / 2.0_RP, gamma2 * dnorm)
            % %end if

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(delta > 0, "DELTA > 0", srname);
            end

        end

    end
end
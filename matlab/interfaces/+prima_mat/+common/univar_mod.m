classdef univar_mod
    %--------------------------------------------------------------------------------------------------%
    % This module implements functions that approximately optimizes univariate functions. They are used
    % in NEWUOA (TRSAPP, BIGLAG, and BIGDEN) and BOBYQA (TRSBOX).
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net) based on Powell's code.
    %
    % Started: January 2021
    %
    % Last Modified: Sunday, April 24, 2022 PM02:15:25
    %
    % N.B.:
    % 0. Here, the optimization is performed using only function values by sampling the objective
    % function on a grid. Indeed, in NEWUOA and BOBYQA, the objective function optimized here is either
    % a trigonometric function or a rational function. So derivatives can be used, but a simple grid
    % search suffices because highly precise solutions are not necessary.
    % 1. Both CIRCLE_MIN and CIRCLE_MAXABS require an input GRID_SIZE, the size of the grid used in
    % the search. Powell chose GRID_SIZE = 50 in NEWUOA. MAGICALLY, this number works the best for
    % NEWUOA in tests on CUTest problems. Larger (e.g., 60, 100) or smaller (e.g., 20, 40) values will
    % worsen the performance of NEWUOA. Why?
    %--------------------------------------------------------------------------------------------------%

    methods
        function angle = circle_min(~, fun, args, grid_size)
            %--------------------------------------------------------------------------------------------------%
            % This function seeks an approximate minimizer of a 2*PI-periodic function FUN(X, ARGS), where the
            % scalar X is the decision variable, and ARGS is a given vector of parameters. It evaluates the
            % function at an evenly distributed "grid" on [0, 2*PI], the number of grid points being GRID_SIZE.
            % Then it takes the grid point with the least value of FUN, and improves the point by a step that
            % minimizes the quadratic that interpolates FUN on this point and its two nearest neighbours.
            % The objective function FUN can represent the parametrization of a function defined on the circle,
            % which explains the name of this function.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = prima_mat.common.consts_mod();
            debug_obj = prima_mat.common.debug_mod();
            infnan_obj = prima_mat.common.infnan_mod();
            linalg_obj = prima_mat.common.linalg_mod();

            % Inputs



            % Outputs
            angle = NaN;

            % Local variables
            srname = "CIRCLE_MIN";


            agrid = NaN(grid_size + 1, 1);


            fgrid = NaN(grid_size, 1);


            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(grid_size >= 3, "GRID_SIZE >= 3", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            agrid(:) = linalg_obj.linspace_r(consts_obj.ZERO, consts_obj.TWO * consts_obj.PI, grid_size + 1); % Size: GRID_SIZE+1; the last entry will be unused
            fgrid(:) = reshape(cell2mat(arrayfun(@(k) fun(agrid(k), args), (1:grid_size), "UniformOutput", false)), [], 1);
            %%MATLAB: fgrid = arrayfun(@(angle) fun(angle, args), agrid(1:grid_size));  % Same shape as `agrid`

            if all(infnan_obj.is_nan(fgrid), 'all')
                angle = consts_obj.ZERO;
                return
            end

            kopt = fix(fortran.minloc(fgrid, 'mask', (~infnan_obj.is_nan(fgrid)), 'dim', 1));
            fopt = fgrid(kopt);
            %%MATLAB: [fopt, kopt] = min(fgrid, [], 'omitnan');
            fprev = fgrid(mod(kopt - 2, grid_size) + 1); % Corresponds to KOPT - 1
            fnext = fgrid(mod(kopt, grid_size) + 1); % Corresponds to KOPT + 1

            step = consts_obj.ZERO;
            if abs(fprev - fnext) > 0
                fprev = fprev - fopt;
                fnext = fnext - fopt;
                step = consts_obj.HALF * (fprev - fnext) / (fprev + fnext);
            end

            if infnan_obj.is_finite(step) && abs(step) > 0
                unit_angle = (consts_obj.TWO * consts_obj.PI) / double(grid_size);
                angle = (double(kopt - 1) + step) * unit_angle;
                % 1. AGRID(KOPT) = (KOPT-1) * UNIT_ANGLE. 2. ANGLE may not be in [0, 2*PI].

            else
                angle = agrid(kopt);
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(infnan_obj.is_finite(angle), "ANGLE is finite", srname);
            end

        end
        function angle = circle_maxabs(~, fun, args, grid_size)
            %--------------------------------------------------------------------------------------------------%
            % This function seeks an approximate maximizer of the absolute value of a 2*PI-periodic function
            % FUN(X, ARGS), where the scalar X is the decision variable, and ARGS is a given vector of parameters.
            % It evaluates the function at an evenly distributed "grid" on [0, 2*PI], the number of grid points
            % being GRID_SIZE. Then it takes the grid point with the largest value of |FUN|, and improves the
            % point by a step that maximizes the absolute value of the quadratic that interpolates FUN on this
            % point and its two nearest neighbours. The objective function FUN can represent the parametrization
            % of a function defined on the circle, which explains the name of this function.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = prima_mat.common.consts_mod();
            debug_obj = prima_mat.common.debug_mod();
            infnan_obj = prima_mat.common.infnan_mod();
            linalg_obj = prima_mat.common.linalg_mod();

            % Inputs



            % Outputs
            angle = NaN;

            % Local variables
            srname = "CIRCLE_MAXABS";


            agrid = NaN(grid_size + 1, 1);


            fgrid = NaN(grid_size, 1);


            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(grid_size >= 3, "GRID_SIZE >= 3", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            agrid(:) = linalg_obj.linspace_r(consts_obj.ZERO, consts_obj.TWO * consts_obj.PI, grid_size + 1); % Size: GRID_SIZE+1; the last entry is not used
            fgrid(:) = reshape(cell2mat(arrayfun(@(k) fun(agrid(k), args), (1:grid_size), "UniformOutput", false)), [], 1);
            %%MATLAB: fgrid = arrayfun(@(angle) fun(angle, args), agrid(1:grid_size));  % Same shape as `agrid`

            if all(infnan_obj.is_nan(fgrid), 'all')
                angle = consts_obj.ZERO;
                return
            end

            kopt = fix(fortran.maxloc(abs(fgrid), 'mask', (~infnan_obj.is_nan(fgrid)), 'dim', 1));
            %%MATLAB: [~, kopt] = max(abs(fgrid), [], 'omitnan');
            fopt = fgrid(kopt);
            fprev = fgrid(mod(kopt - 2, grid_size) + 1); % Corresponds to KOPT - 1
            fnext = fgrid(mod(kopt, grid_size) + 1); % Corresponds to KOPT + 1

            step = consts_obj.ZERO;
            if abs(fprev - fnext) > 0
                fprev = fprev - fopt;
                fnext = fnext - fopt;
                step = consts_obj.HALF * (fprev - fnext) / (fprev + fnext);
            end

            if infnan_obj.is_finite(step) && abs(step) > 0
                unit_angle = (consts_obj.TWO * consts_obj.PI) / double(grid_size);
                angle = (double(kopt - 1) + step) * unit_angle;
                % 1. AGRID(KOPT) = (KOPT-1) * UNIT_ANGLE. 2. ANGLE may not be in [0, 2*PI].

            else
                angle = agrid(kopt);
            end

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(infnan_obj.is_finite(angle), "ANGLE is finite", srname);
            end

        end
        function x = interval_max(~, fun, lb, ub, args, grid_size)
            %--------------------------------------------------------------------------------------------------%
            % This function seeks an approximate maximizer of a function F(X, ARGS) for X in [LB, UB], where
            % ARGS is a vector of parameters. It evaluates the function at an evenly distributed "grid" on
            % [LB, UB], the number of grid points being GRID_SIZE. Then it takes the grid point with the largest
            % value of FUN, and improves the point by a step that maximizes the quadratic that interpolates
            % FUN on this point and its two nearest neighbours unless the point is LB or UB.
            %--------------------------------------------------------------------------------------------------%
            consts_obj = prima_mat.common.consts_mod();
            debug_obj = prima_mat.common.debug_mod();
            infnan_obj = prima_mat.common.infnan_mod();
            linalg_obj = prima_mat.common.linalg_mod();

            % Inputs



            % Outputs
            x = NaN;

            % Local variables
            srname = "INTERVAL_MAX";


            fgrid = NaN(grid_size, 1);


            xgrid = NaN(grid_size, 1);

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(infnan_obj.is_finite(lb) && infnan_obj.is_finite(ub) && lb <= ub, "LB <= UB and they are finite", srname);
                debug_obj.assert(grid_size >= 3, "GRID_SIZE >= 3", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            if ub <= lb
                x = lb;
                return
            end

            xgrid(:) = linalg_obj.linspace_r(lb, ub, grid_size);
            fgrid(:) = reshape(cell2mat(arrayfun(@(k) fun(xgrid(k), args), (1:grid_size), "UniformOutput", false)), [], 1);
            %%MATLAB: fgrid = arrayfun(@(x) fun(x, args), xgrid(1:grid_size));  % Same shape as `xgrid`

            if all(infnan_obj.is_nan(fgrid), 'all')
                x = lb;
                return
            end

            kopt = fix(fortran.maxloc(fgrid, 'mask', (~infnan_obj.is_nan(fgrid)), 'dim', 1));
            fopt = fgrid(kopt);
            %%MATLAB: [fopt, kopt] = min(fgrid, [], 'omitnan');

            if kopt == 1
                x = lb;
            elseif kopt == grid_size
                x = ub;
            else
                fprev = fgrid(kopt - 1);
                fnext = fgrid(kopt + 1);
                step = consts_obj.ZERO;
                if abs(fprev - fnext) > 0
                    step = consts_obj.HALF * ((fnext - fprev) / (fopt + fopt - fprev - fnext));
                end
                if infnan_obj.is_finite(step) && abs(step) > 0
                    x = lb + (ub - lb) * (double(kopt - 1) + step) / double(grid_size - 1);
                    % N.B.: 1. XGRID(KOPT) = LB + (UB-LB)*(KOPT - 1)/(GRID_SIZE -1)
                    % 2. XGRID(KOPT-1) <= X <= XGRID(KOPT+1), as X maximizes the quadratic interpolant.

                else
                    x = xgrid(kopt);
                end
            end

            %====================%
            %  Calculation ends  %
            %====================%

            % Postconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(lb <= x && x <= ub, "LB <= X <= UB", srname);
            end

        end

    end
end
classdef xinbd_mod


    methods
        function x = xinbd(~, xbase, step, xl, xu, sl, su)
            %--------------------------------------------------------------------------------------------------%
            % This function sets X to XBASE + STEP, paying careful attention to the following bounds.
            % 1. XBASE is a point between XL and XU (guaranteed);
            % 2. STEP is a step between SL and SU (may be with rounding errors);
            % 3. SL = XL - XBASE, SU = XU - XBASE;
            % 4. X should be between XL and XU.
            %--------------------------------------------------------------------------------------------------%
            % Common modules
            consts_obj = consts_mod();
            debug_obj = debug_mod();
            linalg_obj = linalg_mod();
            infnan_obj = infnan_mod();


            % Inputs



            % Outputs
            x = NaN(numel(xbase), 1);

            % Local variables
            srname = "XINBD";
            n = NaN;
            s = NaN(numel(xbase), 1);

            % Sizes
            n = fix(numel(xbase));

            % Preconditions
            if consts_obj.DEBUGGING
                debug_obj.assert(all(infnan_obj.is_finite(xbase), 'all'), "SIZE(XBASE) == N, XBASE is finite", srname);
                debug_obj.assert(numel(xl) == n && numel(xu) == n, "SIZE(XL) == N == SIZE(XU)", srname);
                debug_obj.assert(all(xbase >= xl & xbase <= xu, 'all'), "XL <= XBASE <= XU", srname);
                debug_obj.assert(numel(sl) == n && numel(su) == n, "SIZE(SL) == N == SIZE(SU)", srname);
                debug_obj.assert(all(step + 100.0 * consts_obj.EPS_custom * max(consts_obj.ONE, abs(step)) >= sl & step - 100.0 * consts_obj.EPS_custom * max(consts_obj.ONE, abs(step)) <= su, 'all'), "SL <= STEP <= SU", srname);
            end

            %====================%
            % Calculation starts %
            %====================%

            s(:) = max(sl, min(su, step));
            x(:) = max(xl, min(xu, xbase + s));
            x(linalg_obj.trueloc(s <= sl)) = xl(linalg_obj.trueloc(s <= sl));
            x(linalg_obj.trueloc(s >= su)) = xu(linalg_obj.trueloc(s >= su));

            %====================%
            %  Calculation ends  %
            %====================%

            if consts_obj.DEBUGGING
                debug_obj.assert(numel(x) == n && all(x >= xl & x <= xu, 'all'), "SIZE(X) == N, XL <= X <= XU", srname);
                debug_obj.assert(all(x <= xl | step > sl, 'all'), "X == XL if STEP <= SL", srname);
                debug_obj.assert(all(x >= xu | step < su, 'all'), "X == XU if STEP >= SU", srname);
            end

        end

    end
end
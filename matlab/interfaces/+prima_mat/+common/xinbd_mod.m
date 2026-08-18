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



            % Inputs



            % Outputs
            x = NaN(numel(xbase), 1);

            % Local variables


            s = NaN(numel(xbase), 1);

            % Sizes


            % Preconditions


            %====================%
            % Calculation starts %
            %====================%

            s(:) = max(sl, min(su, step));
            x(:) = max(xl, min(xu, xbase + s));
            x(s <= sl) = xl(s <= sl);
            x(s >= su) = xu(s >= su);

            %====================%
            %  Calculation ends  %
            %====================%



        end

    end
end
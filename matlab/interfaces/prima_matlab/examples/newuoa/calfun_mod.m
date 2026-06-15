%--------------------------------------------------------------------------------------------------%
% This is an example to illustrate the usage of the solver.
%
% The objective function is trivial. This is intentional, as the focus is how to use the API.
%--------------------------------------------------------------------------------------------------%
%------------------------- THE MODULE THAT IMPLEMENTS CALFUN, CALLBACK_FCN ------------------------%
classdef calfun_mod

    properties
        RP;
        IK;
    end

    methods
        function obj = calfun_mod()
            obj.RP = class(0.0);
            obj.IK = class(0);
        end
        % Objective function
        function f = calfun(~, x)

            % Inputs


            % Outputs


            f = fortran.power((x(1) - 5.0), 2) + fortran.power((x(2) - 4.0), 2);

        end
        % Callback function
        function terminate = callback_fcn(~, x, f, nf, tr, varargin)


            ipObj = inputParser();
            addParameter(ipObj, 'cstrv', NaN);
            addParameter(ipObj, 'nlconstr', NaN);
            addParameter(ipObj, 'terminate', false);
            parse(ipObj, varargin{:});
            cstrv = ipObj.Results.cstrv;
            nlconstr = ipObj.Results.nlconstr;
            terminate = ipObj.Results.terminate;
            % Suppress compiler warning about unused variable
            % Suppress compiler warning about unused variable

            fprintf(1, 'Best point so far: x = [%6.4f, %6.4f], f = %6.3f, nf = %d, tr = %d\n', x(1), x(2), f, nf, tr);

            terminate = false;

        end

    end
end
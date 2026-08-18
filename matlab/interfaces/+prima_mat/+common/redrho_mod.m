classdef redrho_mod
    %--------------------------------------------------------------------------------------------------%
    % This module provides a function that calculates RHO when it needs to be reduced.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net).
    %
    % Started: September 2021
    %
    % Last Modified: Monday, November 06, 2023 PM07:45:58
    %--------------------------------------------------------------------------------------------------%

    methods
        function rho = redrho(~, rho_in, rhoend)
            %--------------------------------------------------------------------------------------------------%
            % This function calculates RHO when it needs to be reduced.
            % The scheme is shared by UOBYQA, NEWUOA, BOBYQA, LINCOA. For COBYLA, Powell's code reduces RHO by
            % `RHO = HALF * RHO; IF (RHO <= 1.5_RP * RHOEND) RHO = RHOEND`, as specified in (11) of the COBYLA
            % paper. However, this scheme seems to work better, especially after we introduce DELTA.
            %--------------------------------------------------------------------------------------------------%



            rho = NaN;


            %====================%
            % Calculation starts %
            %====================%

            rho_ratio = rho_in / rhoend;

            if rho_ratio > 250.0
                rho = 0.1 * rho_in;
            elseif rho_ratio <= 16.0
                rho = rhoend;
            else
                rho = sqrt(rho_ratio) * rhoend; %rho = sqrt(rho_in * rhoend)
            end

            %====================%
            %  Calculation ends  %
            %====================%



        end

    end
end
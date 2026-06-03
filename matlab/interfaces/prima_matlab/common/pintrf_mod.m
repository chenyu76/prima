classdef pintrf_mod
    %--------------------------------------------------------------------------------------------------%
    % This is a module specifying the abstract interfaces OBJ, OBJCON, and CALLBACK. OBJ evaluates the
    % objective function for unconstrained, bound constrained, and linearly constrained problems; OBJCON
    % evaluates the objective and constraint functions for nonlinearly constrained problems; CALLBACK
    % is a callback function that is called after each iteration of the solvers to report the progress
    % and optionally request termination.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net).
    %
    % Started: July 2020.
    %
    % Last Modified: Friday, December 22, 2023 PM01:23:44
    %--------------------------------------------------------------------------------------------------%

    methods
        function obj = pintrf_mod()
            %%%%%% Users must provide the implementation of OBJ or OBJCON. !!!!!!

        end

    end
end
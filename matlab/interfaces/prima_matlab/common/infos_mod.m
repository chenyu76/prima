classdef infos_mod
    %--------------------------------------------------------------------------------------------------%
    % This is a module defining exit flags.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net).
    %
    % Started: July 2020.
    %
    % Last Modified: Sunday, May 21, 2023 PM03:00:41
    %--------------------------------------------------------------------------------------------------%
    properties
        consts_obj;
        INFO_DFT;
        SMALL_TR_RADIUS;
        FTARGET_ACHIEVED;
        TRSUBP_FAILED;
        MAXFUN_REACHED;
        MAXTR_REACHED;
        NAN_INF_X;
        NAN_INF_F;
        NAN_INF_MODEL;
        NO_SPACE_BETWEEN_BOUNDS;
        DAMAGING_ROUNDING;
        ZERO_LINEAR_CONSTRAINT;
        CALLBACK_TERMINATE;
        INVALID_INPUT;
        ASSERTION_FAILS;
        VALIDATION_FAILS;
        MEMORY_ALLOCATION_FAILS;
    end

    methods
        function obj = infos_mod()
            obj.consts_obj = consts_mod();




















            obj.INFO_DFT = 0;
            obj.SMALL_TR_RADIUS = 0;
            obj.FTARGET_ACHIEVED = 1;
            obj.TRSUBP_FAILED = 2;
            obj.MAXFUN_REACHED = 3;
            obj.MAXTR_REACHED = 20;
            obj.NAN_INF_X = -1;
            obj.NAN_INF_F = -2;
            obj.NAN_INF_MODEL = -3;
            obj.NO_SPACE_BETWEEN_BOUNDS = 6;
            obj.DAMAGING_ROUNDING = 7;
            obj.ZERO_LINEAR_CONSTRAINT = 8;
            obj.CALLBACK_TERMINATE = 30;

            % Stop-codes.
            % The following codes are used by ERROR STOP as stop-codes, which should be default integers.
            obj.INVALID_INPUT = 100;
            obj.ASSERTION_FAILS = 101;
            obj.VALIDATION_FAILS = 102;
            obj.MEMORY_ALLOCATION_FAILS = 103;
        end

    end
end
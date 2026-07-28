classdef memory_mod
    %--------------------------------------------------------------------------------------------------%
    % This module provides subroutines concerning memory management.
    %
    % In particular, the intrinsic ALLOCATE is wrapped into the procedure SAFEALLOC, which may be a
    % controversial practice. We choose to do this because it has helped us a couple of times to locate
    % bugs or problems in our code or even in compilers (e.g., Absoft). See the below for discussions:
    % https://fortran-lang.discourse.group/t/best-practice-of-allocating-memory-in-fortran
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net).
    %
    % Started: July 2020
    %
    % Last Modified: Wednesday, February 28, 2024 AM12:20:23
    %--------------------------------------------------------------------------------------------------%

    methods
        function varargout = cstyle_sizeof(obj, varargin)
            if numel(varargin) == 1
                [varargout{1:nargout}] = obj.size_of_sp(varargin{:});
            elseif numel(varargin) == 1
                [varargout{1:nargout}] = obj.size_of_dp(varargin{:});
            else
                [varargout{1:nargout}] = obj.size_of_qp(varargin{:});
            end
        end
        function varargout = safealloc(obj, varargin)
            if numel(varargin) == 2 && ischar(varargin{1}) && isscalar(varargin{1})
                [varargout{1:nargout}] = obj.alloc_character(varargin{:});
            elseif numel(varargin) == 2 && islogical(varargin{1}) && isvector(varargin{1})
                [varargout{1:nargout}] = obj.alloc_lvector(varargin{2});
            elseif numel(varargin) == 2 && isinteger(varargin{1}) && isvector(varargin{1})
                [varargout{1:nargout}] = obj.alloc_ivector(varargin{2});
            elseif numel(varargin) == 2 && isfloat(varargin{1}) && isvector(varargin{1})
                [varargout{1:nargout}] = obj.alloc_rvector_sp(varargin{2});
            elseif numel(varargin) == 2 && isfloat(varargin{1}) && isvector(varargin{1})
                [varargout{1:nargout}] = obj.alloc_rvector_dp(varargin{2});
            elseif numel(varargin) == 2 && isfloat(varargin{1}) && isvector(varargin{1})
                [varargout{1:nargout}] = obj.alloc_rvector_qp(varargin{2});
            elseif numel(varargin) == 3 && isinteger(varargin{1}) && (~isvector(varargin{1}) && ~isscalar(varargin{1}))
                [varargout{1:nargout}] = obj.alloc_imatrix(varargin{2}, varargin{3});
            elseif numel(varargin) == 3 && isfloat(varargin{1}) && (~isvector(varargin{1}) && ~isscalar(varargin{1}))
                [varargout{1:nargout}] = obj.alloc_rmatrix_sp(varargin{2}, varargin{3});
            elseif numel(varargin) == 3 && isfloat(varargin{1}) && (~isvector(varargin{1}) && ~isscalar(varargin{1}))
                [varargout{1:nargout}] = obj.alloc_rmatrix_dp(varargin{2}, varargin{3});
            else
                [varargout{1:nargout}] = obj.alloc_rmatrix_qp(varargin{2}, varargin{3});
            end
        end
        function y = size_of_sp(~, x)
            %--------------------------------------------------------------------------------------------------%
            % Return the storage size of X in Bytes, X being a REAL(SP) scalar.
            %--------------------------------------------------------------------------------------------------%

            % Inputs

            % Outputs
            y = NaN;

            % We prefer STORAGE_SIZE to C_SIZEOF, because the former is intrinsic while the later requires the
            % intrinsic module ISO_C_BINDING.
            y = fix(whos('x').bytes / numel(x) * 8 / 8); % Y = INT(C_SIZEOF(X), KIND(Y))
        end
        function y = size_of_dp(~, x)
            %--------------------------------------------------------------------------------------------------%
            % Return the storage size of X in Bytes, X being a REAL(DP) scalar.
            %--------------------------------------------------------------------------------------------------%

            % Inputs

            % Outputs
            y = NaN;

            y = fix(whos('x').bytes / numel(x) * 8 / 8);
        end
        function y = size_of_qp(~, x)
            %--------------------------------------------------------------------------------------------------%
            % Return the storage size of X in Bytes, X being a REAL(QP) scalar.
            %--------------------------------------------------------------------------------------------------%

            % Inputs

            % Outputs
            y = NaN;

            y = fix(whos('x').bytes / numel(x) * 8 / 8);
        end
        function x = alloc_rvector_sp(~, n)
            %--------------------------------------------------------------------------------------------------%
            % Allocate space for an allocatable REAL(SP) vector X, whose size is N after allocation.
            %--------------------------------------------------------------------------------------------------%

            debug_obj = prima_mat.common.debug_mod();

            % Inputs


            % Outputs


            % Local variables

            srname = "ALLOC_RVECTOR_SP";

            % Preconditions (checked even not debugging)
            debug_obj.validate(n >= 0, "N >= 0", srname);

            % According to the Fortran 2003 standard, when a procedure is invoked, any allocated ALLOCATABLE
            % object that is an actual argument associated with an INTENT(OUT) ALLOCATABLE dummy argument is
            % deallocated. So the following line is unnecessary since F2003 as X is INTENT(OUT):
            % %if (allocated(x)) deallocate (x)
            % Allocate memory for X. Initialize X to a compiler-independent strange value.
            x = NaN(n, 1); alloc_status = 0;
            x(:) = -realmax; % Costly if X is of a large size.
            % N.B.: Do not write ALLOCATE (X(1:N), STAT=ALLOC_STATUS, SOURCE=-HUGE(X)), because
            % 1. It is invalid to put X in the SOURCE specifier when it is being allocated;
            % 2. Absoft does not support the SOURCE keyword as of 2022.

            % Postconditions (checked even not debugging)
            debug_obj.validate(alloc_status == 0, "Memory allocation succeeds (ALLOC_STATUS == 0)", srname);
            debug_obj.validate(exist('x', 'var'), "X is allocated", srname);
            debug_obj.validate(numel(x) == n, "SIZE(X) == N", srname);
            debug_obj.validate(size(x, 1) == n, "LBOUND(X, 1) == 1, UBOUND(X, 1) == N", srname);
        end
        function x = alloc_rmatrix_sp(~, m, n)
            %--------------------------------------------------------------------------------------------------%
            % Allocate space for an allocatable REAL(SP) matrix X, whose size is (M, N) after allocation.
            %--------------------------------------------------------------------------------------------------%

            debug_obj = prima_mat.common.debug_mod();

            % Inputs


            % Outputs


            % Local variables

            srname = "ALLOC_RMATRIX_SP";

            % Preconditions (checked even not debugging)
            debug_obj.validate(m >= 0 && n >= 0, "M >= 0, N >= 0", srname);

            %if (allocated(x)) deallocate (x)  ! Unnecessary in F03 since X is INTENT(OUT)
            % Allocate memory for X. Initialize X to a compiler-independent strange value.
            x = NaN(m, n); alloc_status = 0; % Absoft does not support the SOURCE keyword as of 2022.
            x = repmat(-realmax, size(x)); % Costly if X is of a large size.

            % Postconditions (checked even not debugging)
            debug_obj.validate(alloc_status == 0, "Memory allocation succeeds (ALLOC_STATUS == 0)", srname);
            debug_obj.validate(exist('x', 'var'), "X is allocated", srname);
            debug_obj.validate(size(x, 1) == m && size(x, 2) == n, "SIZE(X) == [M, N]", srname);
            debug_obj.validate(size(x, 1) == m, "LBOUND(X, 1) == 1, UBOUND(X, 1) == M", srname);
            debug_obj.validate(size(x, 2) == n, "LBOUND(X, 2) == 1, UBOUND(X, 2) == N", srname);
        end
        function x = alloc_rvector_dp(~, n)
            %--------------------------------------------------------------------------------------------------%
            % Allocate space for an allocatable REAL(DP) vector X, whose size is N after allocation.
            %--------------------------------------------------------------------------------------------------%

            debug_obj = prima_mat.common.debug_mod();

            % Inputs


            % Outputs


            % Local variables

            srname = "ALLOC_RVECTOR_DP";

            % Preconditions (checked even not debugging)
            debug_obj.validate(n >= 0, "N >= 0", srname);

            % According to the Fortran 2003 standard, when a procedure is invoked, any allocated ALLOCATABLE
            % object that is an actual argument associated with an INTENT(OUT) ALLOCATABLE dummy argument is
            % deallocated. So the following line is unnecessary since F2003 as X is INTENT(OUT):
            % %if (allocated(x)) deallocate (x)
            % Allocate memory for X. Initialize X to a compiler-independent strange value.
            x = NaN(n, 1); alloc_status = 0; % Absoft does not support the SOURCE keyword as of 2022.
            x(:) = -realmax; % Costly if X is of a large size.

            % Postconditions (checked even not debugging)
            debug_obj.validate(alloc_status == 0, "Memory allocation succeeds (ALLOC_STATUS == 0)", srname);
            debug_obj.validate(exist('x', 'var'), "X is allocated", srname);
            debug_obj.validate(numel(x) == n, "SIZE(X) == N", srname);
            debug_obj.validate(size(x, 1) == n, "LBOUND(X, 1) == 1, UBOUND(X, 1) == N", srname);
        end
        function x = alloc_rmatrix_dp(~, m, n)
            %--------------------------------------------------------------------------------------------------%
            % Allocate space for an allocatable REAL(DP) matrix X, whose size is (M, N) after allocation.
            %--------------------------------------------------------------------------------------------------%

            debug_obj = prima_mat.common.debug_mod();

            % Inputs


            % Outputs


            % Local variables

            srname = "ALLOC_RMATRIX_DP";

            % Preconditions (checked even not debugging)
            debug_obj.validate(m >= 0 && n >= 0, "M >= 0, N >= 0", srname);

            %if (allocated(x)) deallocate (x)  ! Unnecessary in F03 since X is INTENT(OUT)
            % Allocate memory for X. Initialize X to a compiler-independent strange value.
            x = NaN(m, n); alloc_status = 0; % Absoft does not support the SOURCE keyword as of 2022.
            x = repmat(-realmax, size(x)); % Costly if X is of a large size.

            % Postconditions (checked even not debugging)
            debug_obj.validate(alloc_status == 0, "Memory allocation succeeds (ALLOC_STATUS == 0)", srname);
            debug_obj.validate(exist('x', 'var'), "X is allocated", srname);
            debug_obj.validate(size(x, 1) == m && size(x, 2) == n, "SIZE(X) == [M, N]", srname);
            debug_obj.validate(size(x, 1) == m, "LBOUND(X, 1) == 1, UBOUND(X, 1) == M", srname);
            debug_obj.validate(size(x, 2) == n, "LBOUND(X, 2) == 1, UBOUND(X, 2) == N", srname);
        end
        function x = alloc_rvector_qp(~, n)
            %--------------------------------------------------------------------------------------------------%
            % Allocate space for an allocatable REAL(QP) vector X, whose size is N after allocation.
            %--------------------------------------------------------------------------------------------------%

            debug_obj = prima_mat.common.debug_mod();

            % Inputs


            % Outputs


            % Local variables

            srname = "ALLOC_RVECTOR_QP";

            % Preconditions (checked even not debugging)
            debug_obj.validate(n >= 0, "N >= 0", srname);

            % According to the Fortran 2003 standard, when a procedure is invoked, any allocated ALLOCATABLE
            % object that is an actual argument associated with an INTENT(OUT) ALLOCATABLE dummy argument is
            % deallocated. So the following line is unnecessary since F2003 as X is INTENT(OUT):
            % %if (allocated(x)) deallocate (x)
            % Allocate memory for X. Initialize X to a compiler-independent strange value.
            x = NaN(n, 1); alloc_status = 0; % Absoft does not support the SOURCE keyword as of 2022.
            x(:) = -realmax; % Costly if X is of a large size.

            % Postconditions (checked even not debugging)
            debug_obj.validate(alloc_status == 0, "Memory allocation succeeds (ALLOC_STATUS == 0)", srname);
            debug_obj.validate(exist('x', 'var'), "X is allocated", srname);
            debug_obj.validate(numel(x) == n, "SIZE(X) == N", srname);
            debug_obj.validate(size(x, 1) == n, "LBOUND(X, 1) == 1, UBOUND(X, 1) == N", srname);
        end
        function x = alloc_rmatrix_qp(~, m, n)
            %--------------------------------------------------------------------------------------------------%
            % Allocate space for an allocatable REAL(QP) matrix X, whose size is (M, N) after allocation.
            %--------------------------------------------------------------------------------------------------%

            debug_obj = prima_mat.common.debug_mod();

            % Inputs


            % Outputs


            % Local variables

            srname = "ALLOC_RMATRIX_QP";

            % Preconditions (checked even not debugging)
            debug_obj.validate(m >= 0 && n >= 0, "M >= 0, N >= 0", srname);

            % %if (allocated(x)) deallocate (x)  ! Unnecessary in F03 since X is INTENT(OUT)
            % Allocate memory for X. Initialize X to a compiler-independent strange value.
            x = NaN(m, n); alloc_status = 0; % Absoft does not support the SOURCE keyword as of 2022.
            x = repmat(-realmax, size(x)); % Costly if X is of a large size.

            % Postconditions (checked even not debugging)
            debug_obj.validate(alloc_status == 0, "Memory allocation succeeds (ALLOC_STATUS == 0)", srname);
            debug_obj.validate(exist('x', 'var'), "X is allocated", srname);
            debug_obj.validate(size(x, 1) == m && size(x, 2) == n, "SIZE(X) == [M, N]", srname);
            debug_obj.validate(size(x, 1) == m, "LBOUND(X, 1) == 1, UBOUND(X, 1) == M", srname);
            debug_obj.validate(size(x, 2) == n, "LBOUND(X, 2) == 1, UBOUND(X, 2) == N", srname);
        end
        function x = alloc_lvector(~, n)
            %--------------------------------------------------------------------------------------------------%
            % Allocate space for an allocatable LOGICAL vector X, whose size is N after allocation.
            %--------------------------------------------------------------------------------------------------%

            debug_obj = prima_mat.common.debug_mod();

            % Inputs


            % Outputs


            % Local variables

            srname = "ALLOC_LVECTOR";

            % Preconditions (checked even not debugging)
            debug_obj.validate(n >= 0, "N >= 0", srname);

            % %if (allocated(x)) deallocate (x)  ! Unnecessary in F03 since X is INTENT(OUT)
            % Allocate memory for X. Initialize X to a compiler-independent value.
            x = false(n, 1); alloc_status = 0; % Absoft does not support the SOURCE keyword as of 2022.
            x(:) = false; % Costly if X is of a large size.

            % Postconditions (checked even not debugging)
            debug_obj.validate(alloc_status == 0, "Memory allocation succeeds (ALLOC_STATUS == 0)", srname);
            debug_obj.validate(exist('x', 'var'), "X is allocated", srname);
            debug_obj.validate(numel(x) == n, "SIZE(X) == N", srname);
            debug_obj.validate(size(x, 1) == n, "LBOUND(X, 1) == 1, UBOUND(X, 1) == N", srname);
        end
        function x = alloc_ivector(~, n)
            %--------------------------------------------------------------------------------------------------%
            % Allocate space for an allocatable INTEGER(IK) vector X, whose size is N after allocation.
            %--------------------------------------------------------------------------------------------------%

            debug_obj = prima_mat.common.debug_mod();

            % Inputs


            % Outputs


            % Local variables

            srname = "ALLOC_IVECTOR";

            % Preconditions (checked even not debugging)
            debug_obj.validate(n >= 0, "N >= 0", srname);

            % %if (allocated(x)) deallocate (x)  ! Unnecessary in F03 since X is INTENT(OUT)
            % Allocate memory for X. Initialize X to a compiler-independent strange value.
            x = NaN(n, 1); alloc_status = 0; % Absoft does not support the SOURCE keyword as of 2022.
            x(:) = -intmax('int32'); % Costly if X is of a large size.

            % Postconditions (checked even not debugging)
            debug_obj.validate(alloc_status == 0, "Memory allocation succeeds (ALLOC_STATUS == 0)", srname);
            debug_obj.validate(exist('x', 'var'), "X is allocated", srname);
            debug_obj.validate(numel(x) == n, "SIZE(X) == N", srname);
            debug_obj.validate(size(x, 1) == n, "LBOUND(X, 1) == 1, UBOUND(X, 1) == N", srname);
        end
        function x = alloc_imatrix(~, m, n)
            %--------------------------------------------------------------------------------------------------%
            % Allocate space for a INTEGER(IK) matrix X, whose size is (M, N) after allocation.
            %--------------------------------------------------------------------------------------------------%

            debug_obj = prima_mat.common.debug_mod();

            % Inputs


            % Outputs


            % Local variables

            srname = "ALLOC_IMATRIX";

            % Preconditions (checked even not debugging)
            debug_obj.validate(m >= 0 && n >= 0, "M >= 0, N >= 0", srname);

            % %if (allocated(x)) deallocate (x)  ! Unnecessary in F03 since X is INTENT(OUT)
            % Allocate memory for X. Initialize X to a compiler-independent strange value.
            x = NaN(m, n); alloc_status = 0; % Absoft does not support the SOURCE keyword as of 2022.
            x = repmat(-intmax('int32'), size(x)); % Costly if X is of a large size.

            % Postconditions (checked even not debugging)
            debug_obj.validate(alloc_status == 0, "Memory allocation succeeds (ALLOC_STATUS == 0)", srname);
            debug_obj.validate(exist('x', 'var'), "X is allocated", srname);
            debug_obj.validate(size(x, 1) == m && size(x, 2) == n, "SIZE(X) == [M, N]", srname);
            debug_obj.validate(size(x, 1) == m, "LBOUND(X, 1) == 1, UBOUND(X, 1) == M", srname);
            debug_obj.validate(size(x, 2) == n, "LBOUND(X, 2) == 1, UBOUND(X, 2) == N", srname);
        end
        function x = alloc_character(~, n)
            %--------------------------------------------------------------------------------------------------%
            % Allocate space for an allocatable character X, whose length is N after allocation.
            % N.B.: Here, we implement only the version with N being the default integer, even if IK = INT16. It
            % is unsafe to use INT16 as the length of a character variable. It may cause overflow in real2str,
            % as a double-precision vector of length ~3500 would be printed as a string longer than 65536.
            % On most modern platforms, the default integer kind is INT32, which is enough for printing
            % double-precision vectors of size ~ 10^8, being sufficient for this project.
            %--------------------------------------------------------------------------------------------------%
            debug_obj = prima_mat.common.debug_mod();

            % Inputs


            % Outputs


            % Local variables

            srname = "ALLOC_CHARACTER";

            % Preconditions (checked even not debugging)
            debug_obj.validate(n >= 0, "N >= 0", srname);

            % %if (allocated(x)) deallocate (x)  ! Unnecessary in F03 since X is INTENT(OUT)
            % Allocate memory for X. Initialize X to a compiler-independent value.
            x = ""; alloc_status = 0; % Absoft does not support the SOURCE keyword as of 2022.
            x = strjoin(repmat(" ", 1, n), ""); % Costly if X is of a large size.

            % Postconditions (checked even not debugging)
            debug_obj.validate(alloc_status == 0, "Memory allocation succeeds (ALLOC_STATUS == 0)", srname);
            debug_obj.validate(exist('x', 'var'), "X is allocated", srname);
            debug_obj.validate(strlength(x) == n, "LEN(X) == N", srname);
        end

    end
end
classdef consts_mod
    %--------------------------------------------------------------------------------------------------%
    % This is a module defining some constants.
    %
    % Coded by Zaikun ZHANG (www.zhangzk.net).
    %
    % Started: July 2020
    %
    % Last Modified: Mon 16 Feb 2026 03:39:56 PM CET
    %--------------------------------------------------------------------------------------------------%
    properties
        iso_fortran_env_obj;
        DEBUGGING;
        IK_DFT;
        RP_DFT;
        IK;
        RP;
        ZERO;
        ONE;
        TWO;
        HALF;
        QUART;
        TEN;
        TENTH;
        PI;
        EPS;
        REALMIN;
        REALMAX;
        MAXPOW10;
        TINYCV;
        FUNCMAX;
        CONSTRMAX;
        BOUNDMAX;
        SYMTOL_DFT;
        ORTHTOL_DFT;
        RHOBEG_DFT;
        RHOEND_DFT;
        FTARGET_DFT;
        CTOL_DFT;
        CWEIGHT_DFT;
        ETA1_DFT;
        ETA2_DFT;
        GAMMA1_DFT;
        GAMMA2_DFT;
        IPRINT_DFT;
        MAXFUN_DIM_DFT;
        MAXHISTMEM;
        MIN_MAXFILT;
        MAXFILT_DFT;
    end
    properties (Access = private)
        HALF_MAXPOW10;
        MHM;
    end

    methods
        function obj = consts_mod()
            %--------------------------------------------------------------------------------------------------%
            % Remarks:
            %
            % 1. REAL*4, REAL*8, INTEGER*4, INTEGER*8 are not Fortran standard expressions. Do not use them!
            %
            % 2. Never use KIND with a literal value, e.g., REAL(KIND = 8), because Fortran standards never
            % define what KIND = 8 means. There is NO guarantee that REAL(KIND = 8) will be legal, let alone
            % being double precision.
            %
            % 3. Fortran standard (as of F2018) specifies the following for types INTEGER and REAL.
            %
            %    - A processor shall provide ONE OR MORE representation methods that define sets of values for
            %    data of type integer; if the kind type parameter is not specified, the default kind value is
            %    KIND(0) and the type specified is DEFAULT INTEGER.
            %    - A processor shall provide TWO OR MORE approximation methods that define sets of values for
            %    data of type real; if the type keyword REAL is specified and the kind type parameter is not
            %    specified, the default kind value is KIND (0.0) and the type specified is DEFAULT REAL; If the
            %    type keyword DOUBLE PRECISION is specified, the kind value is KIND (0.0D0) and the type
            %    specified is DOUBLE PRECISION real; the decimal precision of the double precision real
            %    approximation method shall be greater than that of the default real method.
            %    - Default integer, default real, and default logical all occupy one storage unit. Double
            %    precision and complex occupy two storage units and double complex requires four storage units.
            %
            %    In other words, the standard only imposes that the following three types should be supported:
            %    - INTEGER(KIND(0)), i.e., default integer,
            %    - REAL(KIND(0.0)), i.e., default real (single-precision real),
            %    - REAL(KIND(0.0D0)), i.e., double-precision real.
            %
            %    Moreover, the following should be noted.
            %
            %    - Other types of INTEGER/REAL may not be available on all platforms (e.g., nvfortran 23.3 and
            %    flang 15.0.3 do not support REAL128).
            %    - The standard does not specify the range of the default integer. However, if the default real
            %    occupies 32 bits, which is normally the case, then the default integer occupies also 32 bits,
            %    and hence the range is probably [-2^32, 2^31-1], approximately [-2*10^9, 2*10^9].
            %    - The standard does not specify what the range and precision of the default real or the
            %    double-precision real, except that KIND(0.0D0) should have a greater precision than KIND(0.0)
            %    --- no requirement about the range.
            %
            %    Consequently, the following should be observed in all Fortran code.
            %
            %    - DO NOT use any kind parameter other than IK, IK_DFT, RP, RP_DFT, SP, or DP, unless you are
            %    sure that it is supported by your platform.
            %    - DO NOT make any assumption on the range of INTEGER, REAL, or REAL(0.0D0) unless you are sure.
            %    - Be cautious about OVERFLOW! In particular, for integers working as the lower/upper limit of
            %    arrays, overflow can lead to Segmentation Faults!
            %--------------------------------------------------------------------------------------------------%

            % Integer and real kinds. Unsupported kinds are negative.
            obj.iso_fortran_env_obj = fortran.iso_fortran_env();


            % Standard IO units

















            obj.DEBUGGING = (false); % Whether we are in debugging mode
            obj.IK_DFT = class(0); % Default integer kind
            obj.RP_DFT = class(0.0); % Default real kind

            % Define the integer kind to be used in the Fortran code.
            obj.IK = obj.IK_DFT;

            % Define the real kind to be used in the Fortran code.
            obj.RP = obj.iso_fortran_env_obj.real64;

            % Define some frequently used numbers.
            obj.ZERO = 0.0;
            obj.ONE = 1.0;
            obj.TWO = 2.0;
            obj.HALF = 0.5;
            obj.QUART = 0.25;
            obj.TEN = 10.0;
            obj.TENTH = 0.1;
            obj.PI = 3.141592653589793;

            % EPS is the machine epsilon, namely the smallest floating-point number such that 1.0 + EPS > 1.0.
            obj.EPS = eps(class(obj.ZERO));
            % REALMIN is the smallest positive normalized floating-point number, which is 2^(-1022) ~ 2.225E-308
            % for IEEE double precision. Taking double precision as an example, REALMIN in other languages:
            % MATLAB: realmin or realmin('double')
            % Python: numpy.finfo(numpy.float64).tiny
            % Julia: realmin(Float64)
            % R: double.xmin
            obj.REALMIN = realmin;
            % REALMAX is the largest positive floating-point number, which is 2^1023 * (2 - EPS) ~ 1.797E308
            % for IEEE double precision. Taking double precision as an example, REALMAX in other languages:
            % MATLAB: realmax or realmax('double')
            % Python: numpy.finfo(numpy.float64).max
            % Julia: realmax(Float64)
            % R: double.xmax
            obj.REALMAX = realmax;

            obj.MAXPOW10 = floor(log10(realmax(class(obj.ZERO))));
            obj.HALF_MAXPOW10 = floor(double(obj.MAXPOW10) / 2.0);

            % TINYCV is used in LINCOA. Powell set TINYCV = 1.0D-60. What about setting TINYCV = REALMIN?
            % N.B.: The `if` is a workaround for the following issues with LLVM flang 19.0.0 and nvfortran 24.3.0:
            % https://fortran-lang.discourse.group/t/flang-new-19-0-warning-overflow-on-power-with-integer-exponent/7801
            % https://forums.developer.nvidia.com/t/bug-of-nvfortran-24-3-0-fort1-terminated-by-signal-11/289026
            obj.TINYCV = obj.TEN ^ max(-60, -obj.MAXPOW10);
            % FUNCMAX is used in the moderated extreme barrier. All function values are projected to the
            % interval [-FUNCMAX, FUNCMAX] before passing to the solvers, and NaN is replaced with FUNCMAX.
            % CONSTRMAX plays a similar role for constraints.
            obj.FUNCMAX = obj.TEN ^ max(4, min(30, obj.HALF_MAXPOW10));
            obj.CONSTRMAX = obj.FUNCMAX;
            % Any bound with an absolute value at least BOUNDMAX is considered as no bound.
            obj.BOUNDMAX = obj.QUART * obj.REALMAX;

            % SYMTOL_DFT is the default tolerance for testing symmetry of matrices. It can be set to 0 if the
            % IEEE Standard for Floating-Point Arithmetic (IEEE 754) is respected, particularly if addition and
            % multiplication are commutative. However, as of 20220408, NAG nagfor does not ensure commutativity
            % for REAL128. Indeed, Fortran standards do not enforce IEEE 754, so compilers are not guaranteed to
            % respect it. Hence we set SYMTOL_DFT to a nonzero number when 1 is 1 (although we do not
            % intend to test symmetry in production, it may be tested if 0 is set to 1).
            % Update 20221226: When gfortran 12 is invoked with aggressive optimization options, it is buggy
            % with ALL() and ANY(). We set SYMTOL_DFT to REALMAX to signify this case and disable the check.
            % Update 20221229: ifx 2023.0.0 20221201 cannot ensure symmetry even up to 100*EPS if invoked
            % with aggressive optimization options and if the floating-point numbers are in single precision.
            % Update 20230307: ifx 2023.0.0 20221201 cannot ensure symmetry even up to 10*EPS if invoked with
            % -O3 and if the floating-point numbers are in single precision.
            % Update 20231002: HUAWEI BiSheng Compiler 2.1.0.B010 (flang) cannot ensure symmetry even up to
            % 1.0E2*EPS if invoked with -Ofast and if the floating-point numbers are in single precision.
            % This same is observed for arm-linux-compiler-22.1 on Kunpeng.
            % Update 20260129: AMD AOMP 22.0 cannot ensure symmetry up to TEN*EPS if invoked with -O3 -fast-math
            % and if the floating-point numbers are in single precision.
            %
            % Double or higher precision in released mode
            obj.SYMTOL_DFT = max(obj.TEN * obj.EPS, obj.TEN ^ max(-10, -obj.MAXPOW10));

            % ORTHTOL_DFT is the default tolerance for testing orthogonality of matrices.
            % In some cases, due to compiler bugs, we need to disable the test. We signify such cases by setting
            % ORTHTOL_DFT to REALMAX. For instance, NAG Fortran Compiler is buggy concerning half-precision
            % floating-point numbers before Release 7.2 Build 7201.
            obj.ORTHTOL_DFT = obj.REALMAX;



            % Some default values
            % RHOBEG: initial value of the trust region radius. Should be about one tenth of the greatest
            % expected change to a variable.
            obj.RHOBEG_DFT = obj.ONE;
            % RHOEND: final value of the trust region radius. Should indicate the accuracy required in the final
            % values of the variables.
            obj.RHOEND_DFT = obj.TEN ^ max(-6, -obj.MAXPOW10); % 1.0E-6
            % FTARGET: target value of the objective function. Solvers exit when finding a feasible point with
            % the objective function value no more than FTARGET.
            obj.FTARGET_DFT = -obj.REALMAX;
            % CTOL: tolerance for constraint violation. A point with constraint violation <= CTOL is considered feasible.
            obj.CTOL_DFT = sqrt(obj.EPS);
            % CWEIGHT: weight of constraint violation in the merit function used to select the output point.
            obj.CWEIGHT_DFT = obj.TEN ^ min(8, obj.MAXPOW10); % 1.0E8
            % ETA1: threshold of reduction ratio for shrinking the trust region radius.
            obj.ETA1_DFT = obj.TENTH;
            % ETA2: threshold of reduction ratio for expanding the trust region radius.
            obj.ETA2_DFT = 0.7;
            % GAMMA1: factor for shrinking the trust region radius.
            obj.GAMMA1_DFT = obj.HALF;
            % GAMMA2: factor for expanding the trust region radius.
            obj.GAMMA2_DFT = obj.TWO;
            % IPRINT: printing level. 0 means no printing.
            obj.IPRINT_DFT = 0;
            % MAXFUN_DIM_DFT*N is the maximal number of function evaluations.
            obj.MAXFUN_DIM_DFT = 500;

            % Maximal amount of memory (Byte) allowed for XHIST, FHIST, CONHIST, CHIST, and the filters.
            obj.MHM = 300 * 10 ^ 6;
            % Make sure that MAXHISTMEM does not exceed HUGE(0) to avoid overflow and memory errors.
            obj.MAXHISTMEM = min(obj.MHM, (intmax('int32') - 1) / 2);

            % Maximal length of the filter used in constrained solvers.
            obj.MIN_MAXFILT = 200; % Should be positive; < 200 is not recommended.
            obj.MAXFILT_DFT = 10 * obj.MIN_MAXFILT;
        end
        function varargout = INT16(obj, varargin)
            [varargout{1:nargout}] = obj.iso_fortran_env_obj.int16(varargin{:});
        end
        function varargout = INT32(obj, varargin)
            [varargout{1:nargout}] = obj.iso_fortran_env_obj.int32(varargin{:});
        end
        function varargout = INT64(obj, varargin)
            [varargout{1:nargout}] = obj.iso_fortran_env_obj.int64(varargin{:});
        end
        function varargout = DP(obj, varargin)
            [varargout{1:nargout}] = obj.iso_fortran_env_obj.real64(varargin{:});
        end
        function varargout = SP(obj, varargin)
            [varargout{1:nargout}] = obj.iso_fortran_env_obj.real32(varargin{:});
        end
        function varargout = QP(obj, varargin)
            [varargout{1:nargout}] = obj.iso_fortran_env_obj.real128(varargin{:});
        end
        function varargout = STDIN(obj, varargin)
            [varargout{1:nargout}] = obj.iso_fortran_env_obj.input_unit(varargin{:});
        end
        function varargout = STDOUT(obj, varargin)
            [varargout{1:nargout}] = obj.iso_fortran_env_obj.output_unit(varargin{:});
        end
        function varargout = STDERR(obj, varargin)
            [varargout{1:nargout}] = obj.iso_fortran_env_obj.error_unit(varargin{:});
        end

    end
end
function [a, D] = sppaoli(h, A, type, is_symmetric, iscorrect)
% Sparse Pauli decomposition of sparse matrix A (Error handling for type and is_symmetric only)
% Input parameters:
%   h           - Bandwidth parameter (required): numeric value (direct use), 'off' (internal calculation)
%   A           - Sparse matrix to be decomposed
%   type        - Output type: 'qua' (quaternary index), 'dec' (decimal index)
%   is_symmetric- Matrix symmetry: 'symmetric' (symmetric), 'asymmetric' (asymmetric)
%   iscorrect   - Validation switch: 'open' (enable validation), 'off' (disable validation)
% Output parameters:
%   a           - Pauli decomposition coefficient vector (a0, a1, ...)
%   D           - Pauli matrix index matrix: quaternary if type='qua', decimal if type='dec'
% Authors: Wu Feng, Li Chen, Wu Xuanlong, Guo Yansong, Zhu Li, Yang Yuxiang


%% 1. Required input check (h has no default value, must be provided)
if nargin < 4
    error('Error: At least 4 input arguments required: h, A, type, is_symmetric!');
end

%% 2. Check input matrix A: MUST be a full matrix (not sparse)
if issparse(A)
    error('Error: Input matrix A must be a FULL matrix, sparse matrix is not allowed!');
end


%% 3. Strict validation for parameter h (positive integer or 'off')
if ischar(h)
    if ~strcmp(h, 'off')
        error('Error: Parameter 1 (h) only accepts string ''off''!');
    end
elseif isnumeric(h)
    % Strict check: positive real integer (non-negative, non-fractional, non-complex)
    if ~isreal(h) || h <= 0 || mod(h, 1) ~= 0
        error('Error: Parameter 1 (h) must be a positive integer (e.g., 1,5,10) or ''off''!');
    end
else
    error('Error: Parameter 1 (h) only supports positive integer or string ''off''!');
end

%% 4. Error handling for parameter 3 (type) (original logic preserved)
valid_type = {'qua', 'dec', 'off'};
if ~ismember(type, valid_type)
    error('Error: Parameter 3 (type) must be ''qua'', ''dec'' or ''off''!');
end

%% 5. Error handling for parameter 4 (is_symmetric) (implemented following type parameter logic)
% Core validation: only strings 'symmetric' / 'asymmetric' allowed
valid_sym = {'symmetric', 'asymmetric'};
if ~ismember(is_symmetric, valid_sym)
    error('Error: Parameter 4 (is_symmetric) must be ''symmetric'' (symmetric) or ''asymmetric'' (asymmetric)!');
end

%% 6. Matrix decomposition branch (using valid is_symmetric parameter)
% Map type to integer: qua→0, dec→1, off→2 (consistent with C function)
type_num = 0;
if strcmp(type, 'dec')
    type_num = 1;
elseif strcmp(type, 'off')
    type_num = 2;
end

% Call modified C functions: paoli_s/paoli_is(h, A, type_num)
if strcmp(is_symmetric, 'symmetric')
    [a, D] = paoli_s(h, A, type_num);  % Decomposition for symmetric matrix
else
    [a, D] = paoli_is(h, A, type_num); % Decomposition for asymmetric matrix
end

%% 7. Original validation logic (unchanged, no additional parameter checks)
if nargin >= 5
    if strcmp(iscorrect, 'open')
        % Initialize Pauli matrix library
        PL = cell(4, 1);
        PL{1} = eye(2);
        PL{2} = [1 0; 0 -1];
        PL{3} = [0 1; 1 0];
        PL{4} = [0 -sqrt(-1); sqrt(-1) 0];
        
        N = size(A, 1);
        n = log2(N);
        Ki = sparse(N, N);
        Nd = length(a);
        
        % Reconstruct matrix and calculate error
        for i = 1:Nd
            if strcmp(type, 'qua')
                q = D(:, i) + 1;
            elseif strcmp(type, 'dec')
                q = dec2quat_vector(double(D(i)), n) + 1;
            end
            Gi = sparse(PL{q(1)});
            for k = 2:n
                Gi = kron(Gi, PL{q(k)});
            end
            Ki = Ki + Gi * a(i);
        end
        
        er = norm(A - Ki, 1);
        if er <= 1e-12
            disp(['Pauli decomposition correct, error (1-norm): ', num2str(er)]);
        else
            disp(['Pauli decomposition incorrect, error (1-norm): ', num2str(er)]);
        end
    elseif strcmp(iscorrect, 'off')
        disp('Pauli decomposition validation disabled');
    else
        error('Error: Parameter 5 (iscorrect) must be ''open'' or ''off''!');
    end
end
clear
clc

% 1. Define file path and load matrix
load_path = 'Path to store K1 matrix\K1.mat';% K1 matrix is stored in v7.3 format and must be a full matrix
data = load(load_path);  % Load MAT file into a structure
K1 = data.K1;            % Extract target matrix K1

% 2. New: Check matrix symmetry in main program (migrated from sppaoli function)
% Principle: Use Frobenius norm to check difference between matrix and its transpose; <1e-13 is considered symmetric
symmetry_norm = norm(K1 - K1', 'fro');  % Calculate Frobenius norm of matrix minus its transpose
if symmetry_norm < 1e-13
    symmetry_desc = 'symmetric';  % Condition true: symmetric
else
    symmetry_desc = 'asymmetric';% Condition false: asymmetric
end

% 3. Print result (string variable is defined, no syntax error)
fprintf('Matrix symmetry check result: %s (Frobenius norm difference: %.6e)\n', ...
        symmetry_desc, symmetry_norm);

% 4. Calculate bandwidth
% Ensure matrix is sparse (bandwidth function is more efficient for sparse matrices)
if ~issparse(K1)
    K1_sparse = sparse(K1);
else
    K1_sparse = K1;
end

% Call MATLAB built-in bandwidth function to compute upper and lower bandwidths
Lb_matlab = bandwidth(K1_sparse, 'lower');  % Lower bandwidth (below main diagonal)
Ub_matlab = bandwidth(K1_sparse, 'upper');  % Upper bandwidth (above main diagonal)

% Output bandwidth result (distinguish symmetric/non-symmetric)
if symmetry_norm < 1e-13
    h_matlab = Ub_matlab;  % For symmetric matrix, use upper bandwidth
else
    h_matlab = max(Ub_matlab, Lb_matlab);  % For non-symmetric matrix, use upper bandwidth by default
end


% 5. Initialize time statistics and loop parameters
h = h_matlab;  % Directly use bandwidth calculated by MATLAB
%h = 'off';       % Optional: keep 'off' to compute inside C function
total_time = 0;          % Total runtime accumulator
loop_count = 10;          % Number of loop runs
type_flag = 'off';       % output type 
is_symmetric = symmetry_desc; % symmetry
iscorrect_flag = 'off';  % Validation switch: 'open'/'off'

% 6. Run sparse Pauli decomposition in loop (pass symmetry parameter)
for i = 1:loop_count
    tic;  % Start timer
    
    % Core modification: call sppaoli with 5 parameters: bandwidth, matrix, type, symmetry, validation switch
    [a, DD] = sppaoli(h, K1, type_flag, is_symmetric, iscorrect_flag);
    
    elapsed_time = toc;  % Stop timer, get single runtime
    fprintf('Run %d time: %.6f seconds\n', i, elapsed_time);
    total_time = total_time + elapsed_time;  % Accumulate total time
end

% 7. Calculate and display average time
average_time = total_time / loop_count;
fprintf('------------------------\n');
fprintf('Average runtime: %.6f seconds\n', average_time);
function quatVec = dec2quat_vector(n_dec, n_digits)
    % dec2quat_vector Convert decimal number to fixed-length quaternary row vector
    % Input: 
    %   n_dec    - Non-negative decimal integer
    %   n_digits - Length of the quaternary vector
    % Output: 
    %   quatVec  - Corresponding fixed-length quaternary row vector, elements are integers 0-3
    
    % Input validation
    if n_dec < 0 || floor(n_dec) ~= n_dec
        error('First input must be a non-negative integer');
    end
    if n_digits < 1 || floor(n_digits) ~= n_digits
        error('Second input must be a positive integer');
    end
    
    % Calculate maximum value for n-digit quaternary
    max_value = 4^n_digits - 1;
    if n_dec > max_value
        error(['Exceeds range of ', num2str(n_digits), '-digit quaternary (max: ', num2str(max_value), ')']);
    end
    
    % Initialize output vector
    quatVec = zeros(1, n_digits);
    
    % Special case: quaternary representation of 0
    if n_dec == 0
        return;  % All-zero vector
    end
    
    % Convert to quaternary digits
    temp = n_dec;
    position = n_digits;  % Start from the rightmost (least significant bit)
    
    while temp > 0 && position >= 1
        remainder = mod(temp, 4);
        quatVec(position) = remainder;  % Store current digit
        temp = floor(temp / 4);
        position = position - 1;
    end
end
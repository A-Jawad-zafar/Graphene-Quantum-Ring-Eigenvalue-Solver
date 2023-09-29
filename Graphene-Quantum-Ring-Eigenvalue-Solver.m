% Author: Ahmal Jawad Zafar
% Created: 06/01/2021, 
% Email: azafar2@gsu.edu

% ====================================================================
% Graphene Quantum Ring Eigenvalue Solver
%
% This script numerically solves the eigenvalue equation of a graphene
% quantum ring using the Dirac effective model. The outer radius is fixed,
% and solutions are obtained for different inner radius values.
%
% The bisection method is employed to find solutions where the function
% changes its sign. Solutions exhibiting asymptotic behavior are filtered 
% out to retain only physical solutions of the equation.
% ====================================================================

% ------ Clear workspace and close figures ------
clear all;
close all;

% ========================= Parameters =========================
% Define the parameters needed for the calculations.

h = 0.65821;       % Reduced Planck constant (in eV·fs)
gamma = 3.03;      % The nearest neighbor hopping integral (in eV)
Rin = 5;           % Rin  - Inner radius of the quantum ring (in nm)
Rout = 15;         % Rout - Outer radius of the quantum ring (in nm)
m = -41/2:1:41/2;  % Magnetic quantum number
T = 1e3;           % Total number of input points for eigen function
S = 15;            % Number of solutions this code will go through

% ====================== Calculations ======================
% Perform the calculations for different inner and outer radius values.

% Initialize a parallel processing loop over different inner and outer radius values

Energy_k = zeros(2*S,length(m));
for i = 1:length(Rin)
    Ri = Rin(i);  % Current inner radius value
    % Calculate solutions and store them in E_k
    E_k = Solution(S, Ri, Rout, gamma, m, T); 
    Energy_k(:,:,i) = E_k;  % Store the current solutions in the Energy_k matrix
end

% Save the calculated energy values for the K valley
ENERGY_K_valley = Energy_k;
save('ENERGY_K_valley');


% ====================== Function Definitions ======================

function Energy_k = Solution(S, rin, rout, gamma, m, T)
    % Solution Function
    %
    % Computes the energy levels near the K-Valley of a Graphene Quantum Ring, 
    % considering the asymptotic solutions and zero crossings.
    % 
    % Inputs:
    %   S      - Number of solutions this code will go through.
    %   rin    - Inner radius (in nm)
    %   rout   - Outer radius (in nm)
    %   gamma  - The nearest neighbor hopping integral
    %   m      - Magnetic quantum number
    %   T      - Total Number of input points for Eigen function
    %
    % Output:
    %   Energy_k - Calculated energy levels (in eV)

    % Initialize Position and Negative arrays
    Pos = zeros(S, length(m));
    Neg = zeros(S, length(m));

    % Loop over each magnetic quantum number to determine the Asymptotic solutions.
    for i = 1:length(m)
        MM = m(i);
         [pos,~]=asymptotes(rin./rout,MM,T);
            if  length(Pos(:,i)) == length(pos) 
            Pos(:,i) = flip(pos);                  %epsilon for conduction band
        elseif length(Pos(:,i)) -length(pos) >0 
            a1 = length(Pos(:,i)) - length(pos);
            b1 = zeros(a1,1); 
            pos1 = [pos;b1];
            Pos(:,i) = flip(pos1);                 %epsilon for conduction band
            else
            if length(Pos(:,i)) -length(pos) <0 
            a3 = length(pos)- length(Pos(:,i));  
            pos2 = pos(1:end-a3);   
            Pos(:,i) = flip(pos2) ;                %epsilon for conduction band
            end
            end
         [~,neg]=asymptotes(rin./rout,MM,T);
            if  length(Neg(:,i)) == length(neg)
            Neg(:,i) = (neg);                      %epsilon for valence band
        elseif  length(Neg(:,i)) - length(neg) >0
            a2 = length(Neg(:,i)) - length(neg);
            b2 = zeros(a2,1);
            neg1 = [neg;b2];
            Neg(:,i) = (neg1);                     %epsilon for valence band
        else
            if length(Neg(:,i)) - length(neg) < 0
            a4 =  length(neg) - length(Neg(:,i));
            neg2 = neg(1:end-a4);
            Neg(:,i) = (neg2);                     %epsilon for valence band
            end
            end
    end

    % Merge positive and negative asymptotes
    Asymptotes = [Pos; Neg];
    e_C = zeros(S, length(m));
    e_V = zeros(S, length(m));

    % Adjust array lengths for all solutions to be stored in one matrix.
    for i = 1:length(m)
        MM = m(i);
         [E_c1,~]=zero_crossings(rin./rout,MM,T);
            if  length(e_C(:,i)) == length(E_c1) 
            e_C(:,i) = flip(E_c1);               %epsilon for conduction band
        elseif length(e_C(:,i)) -length(E_c1) >0 
            a5 = length(e_C(:,i)) - length(E_c1);
            b5 = zeros(a5,1); 
            E_c = [E_c1;b5];
            e_C(:,i) = flip(E_c);                 %epsilon for conduction band
            else
            if length(e_C(:,i)) -length(E_c1) <0 
            a6 = length(E_c1)- length(e_C(:,i));  
            E_c = E_c1(1:end-a6);   
            e_C(:,i) = flip(E_c) ;                %epsilon for conduction band
            end
            end
         [~,E_v1]=zero_crossings(rin./rout,MM,T);
            if  length(e_V(:,i)) == length(E_v1)
            e_V(:,i) = (E_v1);                    %epsilon for valence band
        elseif  length(e_V(:,i)) - length(E_v1) >0
            a7 = length(e_V(:,i)) - length(E_v1);
            b6 = zeros(a7,1);
            E_v = [E_v1;b6];
            e_V(:,i) = (E_v);                      %epsilon for valence band
        else
            if length(e_V(:,i)) - length(E_v1) < 0
            a4 =  length(E_v1) - length(e_V(:,i));
            E_v = E_v1(1:end-a4);
            e_V(:,i) = (E_v);                       %epsilon for valence band
            end
            end
    end

    epsilon1 = [e_C; e_V];

    % Filter out asymptotic solutions from the physical solutions.
    [Rows, Cols] = size(epsilon1);
    Epsilon1 = zeros(Rows, Cols);
    for j = 1:Cols
        a = Asymptotes(:, j);
        b = epsilon1(:, j);
        b(ismember(b, a)) = 0;  % Replace asymptotic solutions with zero
        Epsilon1(:, j) = b;     % Epsilon is a unitless variable for energy.
    end

    % Convert unitless variable to energy in eV

    Energy_k = (Epsilon1 .* gamma ./ rout);
end

% ====================== Function Definitions ======================

function [pos, neg] = asymptotes(beta, M, T)
    % Asymptotes Function
    %
    % This function finds the asymptotes of the energy levels in a Graphene 
    % Quantum Ring by looking for zero-crossings in the energy eigenvalue equation.
    %
    % Inputs:
    %   beta - Ratio of inner to outer radius
    %   M    - Magnetic quantum number
    %   T    - Total number of input points for Eigen function
    %
    % Outputs:
    %   pos - Positive energy eigenvalues corresponding to K-Valley
    %   neg - Negative energy eigenvalues corresponding to K-Valley

    E2 = linspace(0, -50, T);
    E1 = linspace(0, 50, T);
    y1 = real(energies(beta, E1, M));
    y2 = real(energies(beta, E2, M));

    % Initialize zero-crossing counts for positive and negative energies in K-valley
    kp1 = 1;
    kn1 = 1;

    for n = 1:T
        if n == T
            % Avoid exceeding the length of the array in (i+1)-th terms
            break
        else
            % K-Valley - Positive energy
            a1p1 = y1(n);
            a2p1 = y1(n + 1);
            Pap1 = a2p1 * a1p1;  % Product of consecutive terms - positive energy - K-Valley

            % K-Valley - Negative energy
            a1n1 = y2(n);
            a2n1 = y2(n + 1);
            Pan1 = a1n1 * a2n1;  % Product of consecutive terms - negative energy - K-Valley

            % Selection of zeros of the energy eigenvalue equation by finding the zero crossings

            % K-Valley Positive eigenenergies
            if Pap1 < 0 && imag(a1p1) == 0 && imag(a2p1) == 0
                Energy_valuesp1(kp1, 1) = (E1(n) + E1(n + 1)) / 2;  % Find energy eigenvalues
                kp1 = kp1 + 1;
            elseif Pap1 > 0
                Energy_valuesp1(kp1, 1) = 0;
                kp1 = kp1 + 1;
            end

            % K-Valley Negative eigenenergies
            if Pan1 < 0 && imag(a1n1) == 0 && imag(a2n1) == 0
                Energy_valuesn1(kn1, 1) = (E2(n) + E2(n + 1)) / 2;  % Find energy eigenvalues
                kn1 = kn1 + 1;
            elseif Pan1 > 0
                Energy_valuesn1(kn1, 1) = 0;
                kn1 = kn1 + 1;
            end
        end
    end

    % Find band gap between corresponding energy levels of conduction and valance band
    pos = nonzeros(Energy_valuesp1);
    neg = nonzeros(Energy_valuesn1);

    % Helper Function to Calculate the values of Denominator of EigenValue
    % Equation

function y = energies(beta, E, M)
    y = zeros(1, length(E));
    for j = 1:length(E)
        er = E(j);
        [Y] = deno(beta, er, M);
        y(j) = Y;
    end
end

% Helper Function to Calculate the Denominator of the Energy Eigenvalue Equation
function [Y] = deno(b, E, m)
    % Bessel function = besselj
    % Neumann function = bessely
    Y = (besselj(m + 1/2, b .* (E)) + besselj(m - 1/2, b .* (E)));

end


end

% ====================== Function Definitions ======================

function [E_c, E_v] = zero_crossings(beta, M, T)
    % zero_crossings Function
    %
    % Computes all values of zero crossing of Eigenvalue Equation at given input parameters.
    %
    % Inputs:
    %   beta - Ratio of inner and outer Radius
    %   M    - Magnetic Quantum number
    %   T    - Total Number of input points for Eigen Equation
    %
    % Outputs:
    %   E_c - Solution in positive energy region
    %   E_v - Solution in negative energy region

    % Initialize Energy Ranges and Values
    E2 = linspace(0, -50, T);
    E1 = linspace(0, 50, T);
    y1 = real(values(beta, E1, M));
    y2 = real(values(beta, E2, M));

    % Initialize zero-crossing counts for K-valley energies
    kp1 = 1;
    kn1 = 1;

    for n = 1:T
        if n == T  % Avoid exceeding the length of the array
            break
        else
            % K-Valley - Positive Energy
            a1p1 = y1(n);
            a2p1 = y1(n + 1);
            Pap1 = a2p1 * a1p1;  % Product of consecutive terms

            % K-Valley - Negative Energy
            a1n1 = y2(n);
            a2n1 = y2(n + 1);
            Pan1 = a1n1 * a2n1;  % Product of consecutive terms

            % Select zeros of the energy eigenvalue equation by finding zero crossings

            % K-Valley Positive Eigenenergies
            if Pap1 < 0 && imag(a1p1) == 0 && imag(a2p1) == 0
                Energy_valuesp1(kp1, 1) = (E1(n) + E1(n + 1)) / 2;  % Find energy eigenvalues
                kp1 = kp1 + 1;
            end

            % K-Valley Negative Eigenenergies
            if Pan1 < 0 && imag(a1n1) == 0 && imag(a2n1) == 0
                Energy_valuesn1(kn1, 1) = (E2(n) + E2(n + 1)) / 2;  % Find energy eigenvalues
                kn1 = kn1 + 1;
            end
        end
    end

    % Store the Corresponding Energy Levels of Conduction and Valance Band
    E_c = Energy_valuesp1;
    E_v = Energy_valuesn1;

    % Nested Function to Compute the Values of the Eigenvalue Equation
    function y = values(beta, E, M)
        y = zeros(1, length(E));
        for j = 1:length(E)
            er = E(j);
            [Y] = eigenequation(beta, er, M);
            y(j) = Y;
        end
    end

    % Nested Eigen Equation Function
    function [Y] = eigenequation(b, E, m)
        % b - beta (Ratio of inner and outer Radius)
        % m - Magnetic Quantum number
        % E - Energy Level
        
        a = (-bessely(m + 1/2, b .* (E)) - bessely(m - 1/2, b .* (E))) ./ ...
            (besselj(m + 1/2, b .* (E)) + besselj(m - 1/2, b .* (E)));

        x1 = a .* besselj(m - 1/2, E);
        x2 = a .* besselj(m + 1/2, E);
        x3 = bessely(m + 1/2, E);
        x4 = bessely(m - 1/2, E);

        Y = x1 + x2 + x3 + x4;
    end
end




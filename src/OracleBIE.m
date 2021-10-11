% Oracle for our toy boundary integral equation example
classdef OracleBIE < ProblemOracle
    properties
        % A function which takes in (x, y) and returns f(x, y)
        forceFunc
    end
    
    methods
        function obj = OracleBIE(nPointsCoarse, nPointsFine, forceFunc)
            obj@ProblemOracle(nPointsCoarse, nPointsFine);
            obj.forceFunc = forceFunc;
        end
        
         % Get the vector size of coarse solutions 
        function [n] = getCoarseSolutionSize(this)
            n = this.m_n_coarse;
        end
        
        % Get the vector size of fine solutions 
        function [n] = getFineSolutionSize(this)
            n = this.m_n_fine; 
        end
        
        % We need to add I / 2 to the operator M which we interpolate
        function [offset_op] = getAffineOffsetOperator(~, reduced_basis)
            n = size(reduced_basis, 2);
            offset_op = eye(n) / 2;
        end
        
        % Solve the coarse system for sample omega
        % Returns solution vector
        function [solution] = solveCoarse(this, omega)
            rho_c = omega';
            N = this.getCoarseSolutionSize();
            [p1, p2, ~, ~] = computeCurve(rho_c, N);
            f = this.forceFunc(p1, p2)';
            solution = solve(rho_c, f, N); 
        end
        
        % Solve the fine system for sample omega
        % Returns vector solution, left hand operator, and auxiliary information
        function [solution, loperator, aux] = solveFine(this, omega)
            rho_c = omega';
            N = this.getFineSolutionSize();
            [p1, p2, ~, ~] = computeCurve(rho_c, N);
            f = this.forceFunc(p1, p2)';
            [solution, loperator] = solve(rho_c, f, N);
            aux = []; % No auxiliary data needed
        end
        
        % Project an operator into the reduced basis
        function [proj_operator] = projectFineOperator(~, reduced_basis, operator)
            proj_operator = reduced_basis' * operator * reduced_basis;
        end
        
        % Sample important sections of all operators in the sample space
        function [samples] = sampleFineOperators(this, sample_space, n_samples)
            N = this.getFineSolutionSize();
            pck = randperm(N, n_samples);
            nRand = size(sample_space, 2);
            
           	samples = zeros(n_samples * N + N, nRand); 
            for g=1:nRand
                
                fprintf("Sampling Fine Operators... (%i / %i)\r", g, nRand);
                
                omega = sample_space(:, g);
                rho_c = omega';
                
                % Old sampling, construct entire integral operator
                %{
                op = -assembleIntegralKernel(rho_c, N);
                
                samples(1:N, g) = diag(op);
                vec_samp = op(:,pck);
                
                samples((N + 1):end, g) = vec_samp(:);
                %}
                
                [op_sub, op_diags] = assembleIntegralKernelColumnSubsampled(rho_c, ...
                    N, pck);
                
                op_sub = -op_sub;
                op_diags = -op_diags;
            
                samples(1:N, g) = op_diags;
                samples((N + 1):end, g) = op_sub(:);
            end

            fprintf("\n")
        end
        
        % No auxiliary data needed
        function computeAuxData(~, ~)
            
        end
        
        % Assemble the right hand of a reduced basis solve
        function [proj_f] = assembleRightHand(this, rb_object, omega, ~)
            rho_c = omega';
            N = this.getFineSolutionSize();
            [p1, p2] = computeCurve(rho_c, N);
            f = this.forceFunc(p1, p2)';
            proj_f = rb_object.m_reduced_basis' * f;
        end
        
        % View the reduced basis
        function viewBasis(this, reduced_basis, figure_id)
            figure(figure_id);
            N = this.getFineSolutionSize();
            theta = (-N/2:N/2)/N * pi * 2;
            mid_point = N/2+1;
            
            nSk = size(reduced_basis, 2);
            nSkrt = ceil(sqrt(nSk));
            ma = max(max(reduced_basis));
            mi = min(min(reduced_basis));
            for ind = 1:nSk
                subplot(nSkrt,nSkrt,ind);
                basis_vec = reduced_basis(:,ind);
                basis_vec_centered = ...
                    [ basis_vec(mid_point:end); basis_vec(1:mid_point)];
                
                plot(theta, basis_vec_centered);
                axis([0, 2 * pi, mi, ma]);
                xticks([-pi -pi/2 0 pi/2 pi]);
                xticklabels({'-\pi','-\pi/2' '0','\pi/2', '\pi'});
                xlim([-pi pi]);
            end
        end
        
        % View the reduced basis
        function viewBasisUncentered(this, reduced_basis, figure_id)
            figure(figure_id);
            N = this.getFineSolutionSize();
            theta = (0:N-1)/N * 2*pi;
            nSk = size(reduced_basis, 2);
            nSkrt = ceil(sqrt(nSk));
            ma = max(max(reduced_basis));
            mi = min(min(reduced_basis));
            for ind = 1:nSk
                subplot(nSkrt,nSkrt,ind);
                basis_vec = reduced_basis(:,ind);
                plot(theta, basis_vec);
                axis([0, 2 * pi, mi, ma]);
                xticks([0 pi/2 pi 3*pi/2 2*pi]);
                xticklabels({'0','\pi/2' '\pi','3\pi/2', '2\pi'});
            end
        end
       
        % Display a solution
        function viewResult(this, result_rb, result_fine, figure_id)
            N = this.getFineSolutionSize();
                       
            min_val = min(min(result_rb), min(result_fine));
            max_val = max(max(result_rb), max(result_fine));
            
            theta = (-N/2:N/2)/N * pi * 2;
            mid_point = N/2+1;
            result_rb_centered = ...
                [result_rb(mid_point:end); result_rb(1:mid_point)];
            result_fine_centered = ...
                [result_fine(mid_point:end); result_fine(1:mid_point)];
            
            figure(figure_id);
            
            subplot(1,2,1); plot(theta, result_rb_centered);
            xticks([-pi -pi/2 0 pi/2 pi]);
            xticklabels({'-\pi','-\pi/2' '0','\pi/2', '\pi'});
            xlim([-pi pi]);
            ylim([min_val, max_val]);
            title('Reduced Basis');
            
            subplot(1,2,2); plot(theta, result_fine_centered);
            xticks([-pi -pi/2 0 pi/2 pi]);
            xticklabels({'-\pi','-\pi/2' '0','\pi/2', '\pi'});
            xlim([-pi pi]);
            ylim([min_val, max_val]);
            title('Ground Truth');
        end
        
        function viewResultUncentered(this, result_rb, result_fine, figure_id)
            N = this.getFineSolutionSize();
            theta = (0:N-1)/N * 2*pi;
            
            min_val = min(min(result_rb), min(result_fine));
            max_val = max(max(result_rb), max(result_fine));
            
            figure(figure_id);
            subplot(1,2,1); plot(theta, result_rb); 
            xticks([0 pi/2 pi 3*pi/2 2*pi]);
            xticklabels({'0','\pi/2' '\pi','3\pi/2', '2\pi'});
            xlim([0 2*pi]);
            ylim([min_val, max_val]);
            title('Reduced Basis');
            subplot(1,2,2); plot(theta, result_fine);
            xticks([0 pi/2 pi 3*pi/2 2*pi]);
            xticklabels({'0','\pi/2' '\pi','3\pi/2', '2\pi'});
            xlim([0 2*pi]);
            ylim([min_val, max_val]);
            title('Ground Truth');
        end
        
        function viewF(this, omega)
            N = this.getFineSolutionSize();
            rho_c = omega';
            [p1, p2, ~, theta] = computeCurve(rho_c, N);
            f = this.forceFunc(p1, p2);
            plot(theta, f);
        end
        
        function viewFineSolution(this, omega)
            N = this.getFineSolutionSize();
            phi = this.solveFine(omega);
            viewSolution(omega, phi, N);
        end
        
        function viewCoarseSolution(this, omega)
            N = this.getCoarseSolutionSize();
            phi = this.solveCoarse(omega);
            viewSolution(omega, phi, N); 
        end
    end
end

function viewSolution(omega, phi, N)
    nK = length(omega);

    rho_c = omega';
    theta_c = (0:nK-1)/nK * 2*pi;

    [p1, p2, ~, theta] = computeCurve(rho_c, N);

    figure; %plot
    subplot(1,2,1); plot(p1,p2); hold on; plot(rho_c.*cos(theta_c), ...
        rho_c.*sin(theta_c),'r+');
    subplot(1,2,2); plot(theta,phi,'-+');
    xticks([0 pi/2 pi 3*pi/2 2*pi]);
    xticklabels({'0','\pi/2' '\pi','3\pi/2', '2\pi'});
    xlim([0 2*pi]);
end

function [solution, lop] = solve(rho_c, f, nPoints)

    M = assembleIntegralKernel(rho_c, nPoints);
    L = eye(nPoints) / 2.0 - M;

    solution = L\f; % solve
    lop = -M; % return the operator H
end

% Compute integral kernel
function [M] = assembleIntegralKernel(rho_c, nPoints)
    N = nPoints;
    [p1, p2, ~, theta] = computeCurve(rho_c, N);

    ks = [0:N/2-1,0,-N/2+1:-1]; 
    f1 = ifft(fft(p1).*(i*ks)); %use FFT to compute 1st derivatives
    f2 = ifft(fft(p2).*(i*ks));

    s1 = ifft(fft(p1).*(i*ks).^2); %use FFT to compute 2nd derivatives
    s2 = ifft(fft(p2).*(i*ks).^2);

    w = ones(size(theta)) * 2*pi/N; %weights
    J = sqrt(f1.^2+f2.^2); %jacobian

    e1 = f2 ./ J; %outward normal derivatives
    e2 = -f1 ./ J;

    crvt = (s2.*f1-s1.*f2) ./ J.^3; %curvature

    x1 = p1'; x2 = p2';
    y1 = p1; y2 = p2;

    r1 = x1*ones(1,N)-ones(N,1)*y1;
    r2 = x2*ones(1,N)-ones(N,1)*y2;
    n1 = ones(N,1) * e1;
    n2 = ones(N,1) * e2;
    r = sqrt( r1.^2 + r2.^2 );
    rn = r1.*n1 + r2.*n2;

    M = 1/(2*pi) * rn ./ r.^2; %form the regular part of the matrix
    bad = find(eye(N)==1);  
    M(bad) = -1/(2*pi) * crvt/2; %singularity correction along diagonal
    M = M .* (ones(N,1)*(J.*w)); %multiple with the correct quadrature weights
end

% Compute integral kernel
function [Msub, diags] = assembleIntegralKernelSubsampled(rho_c, nPoints, pck)
    N = nPoints;
    [p1, p2, ~, theta] = computeCurve(rho_c, N);
    
    [indx1, indx2] = meshgrid(1:nPoints, 1:nPoints);
    indx1 = indx1(pck);
    indx2 = indx2(pck);

    ks = [0:N/2-1,0,-N/2+1:-1]; 
    f1 = ifft(fft(p1).*(1i*ks)); %use FFT to compute 1st derivatives
    f2 = ifft(fft(p2).*(1i*ks));

    s1 = ifft(fft(p1).*(1i*ks).^2); %use FFT to compute 2nd derivatives
    s2 = ifft(fft(p2).*(1i*ks).^2);

    w = ones(size(theta)) * 2*pi/N; %weights
    J = sqrt(f1.^2+f2.^2); %jacobian

    e1 = f2 ./ J; %outward normal derivatives
    e2 = -f1 ./ J;

    crvt = (s2.*f1-s1.*f2) ./ J.^3; %curvature
    
    x1 = p1(indx2);
    x2 = p2(indx2);
    y1 = p1(indx1);
    y2 = p2(indx1);
    
    r1 = x1 - y1;
    r2 = x2 - y2;
    r_E2 = r1 .* r1 + r2 .* r2;
    
    n1 = e1(indx1);
    n2 = e2(indx1);
    rn = r1 .* n1 + r2 .* n2;
    
    Msub = 1/(2*pi) * rn ./ r_E2;
    
    diags = -1/(2*pi) * crvt/2;
    diag_entries = (indx1 == indx2);
    
    Msub(diag_entries) = diags(indx2(diag_entries));
    quad_w = J .* w;
    Msub = Msub .* quad_w(indx1);
    diags = (diags .* quad_w)';
end

% Compute integral kernel
function [Msub, diags] = assembleIntegralKernelColumnSubsampled(rho_c, nPoints, col_pck)
    N = nPoints;
    [p1, p2, ~, theta] = computeCurve(rho_c, N);
    
    [indx1, indx2] = meshgrid(col_pck, 1:nPoints);

    ks = [0:N/2-1,0,-N/2+1:-1]; 
    f1 = ifft(fft(p1).*(1i*ks)); %use FFT to compute 1st derivatives
    f2 = ifft(fft(p2).*(1i*ks));

    s1 = ifft(fft(p1).*(1i*ks).^2); %use FFT to compute 2nd derivatives
    s2 = ifft(fft(p2).*(1i*ks).^2);

    w = ones(size(theta)) * 2*pi/N; %weights
    J = sqrt(f1.^2+f2.^2); %jacobian

    e1 = f2 ./ J; %outward normal derivatives
    e2 = -f1 ./ J;

    crvt = (s2.*f1-s1.*f2) ./ J.^3; %curvature
    
    x1 = p1(indx2);
    x2 = p2(indx2);
    y1 = p1(indx1);
    y2 = p2(indx1);
    
    r1 = x1 - y1;
    r2 = x2 - y2;
    r_E2 = r1 .* r1 + r2 .* r2;
    
    n1 = e1(indx1);
    n2 = e2(indx1);
    rn = r1 .* n1 + r2 .* n2;
    
    Msub = 1/(2*pi) * rn ./ r_E2;
    
    diags = -1/(2*pi) * crvt/2;
    diag_entries = (indx1 == indx2);
    
    Msub(diag_entries) = diags(indx2(diag_entries));
    quad_w = J .* w;
    Msub = Msub .* quad_w(indx1);
    diags = (diags .* quad_w)';
end


function [p1, p2, rho, theta] = computeCurve(rho_c, nPoints)
    N = nPoints;
    rho = interpft(rho_c, N); %interpolate to all grid points
    theta = (0:N-1)/N * 2 * pi;

    p1 = rho.*cos(theta); %compute locations
    p2 = rho.*sin(theta); 
end
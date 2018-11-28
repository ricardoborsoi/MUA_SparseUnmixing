function [U,res2,rmse] = sunsal_tv_lw_sp(M,Y,varargin)


%%  =========================== Outputs ==================================
%
% U  =  [nxN] estimated  X matrix
%
%
% ----------------------------------------------------------------------

%%
%--------------------------------------------------------------
% test for number of required parametres
%--------------------------------------------------------------
if (nargin-length(varargin)) ~= 2
    error('Wrong number of required parameters');
end
% mixing matrix size
[LM,n] = size(M);
% data set size
[L,N] = size(Y);
if (LM ~= L)
    error('mixing matrix M and data set y are inconsistent');
end

%%
%--------------------------------------------------------------
% Set the defaults for the optional parameters
%--------------------------------------------------------------
%

% 'LAMBDA_1'
%  l1 regularization
reg_l1 = 0; % absent

reg_TV = 0; % absent
im_size = []; % image size
tv_type = 'niso'; % non-isotropic TV

% 'AL:ITERS'
% maximum number of AL iteration
% AL_iters = 1000;

% 'MU'
% AL weight
mu = 0.001;

% 'VERBOSE'
% display only sunsal warnings
verbose = 'off';

% 'POSITIVITY'
% Positivity constraint
positivity = 'no';
reg_pos = 0; % absent

% 'ADDONE'
%  Sum-to-one constraint
addone = 'no';
reg_add = 0; % absent

% initialization
U0 = 0;

% true X
true_x = 0;
rmse = 0;

% Read the optional parameters
%--------------------------------------------------------------
if (rem(length(varargin),2)==1)
    error('Optional parameters should always go by pairs');
else
    for i=1:2:(length(varargin)-1)
        switch upper(varargin{i})
            case 'LAMBDA_1'
                lambda_l1 = varargin{i+1};
                if lambda_l1 < 0
                    error('lambda must be positive');
                elseif lambda_l1 > 0
                    reg_l1 = 1;
                end
            case 'LAMBDA_TV'
                lambda_TV = varargin{i+1};
                if lambda_TV < 0
                    error('lambda must be non-negative');
                elseif lambda_TV > 0
                    reg_TV = 1;
                end
            case 'TV_TYPE'
                tv_type = varargin{i+1};
                if ~(strcmp(tv_type,'iso') | strcmp(tv_type,'niso'))
                    error('wrong TV_TYPE');
                end
            case 'IM_SIZE'
                im_size = varargin{i+1};
            case 'AL_ITERS'
                AL_iters = round(varargin{i+1});
                if (AL_iters <= 0 )
                    error('AL_iters must a positive integer');
                end
            case 'POSITIVITY'
                positivity = varargin{i+1};
                if strcmp(positivity,'yes')
                    reg_pos = 1;
                end
            case 'ADDONE'
                addone = varargin{i+1};
                if strcmp(addone,'yes')
                    reg_add = 1;
                end
            case 'MU'
                mu = varargin{i+1};
                if mu <= 0
                    error('mu must be positive');
                end
            case 'VERBOSE'
                verbose = varargin{i+1};
            case 'X0'
                U0 = varargin{i+1};
            case 'TRUE_X'
                XT = varargin{i+1};
                true_x = 1;
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end;
    end;
end

% test for true data size correctness
if true_x
    [nr nc] = size(XT);
    if (nr ~= n) | (nc ~= N)
        error('wrong image size')
    end
end
    n_lin = im_size(1);
    n_col = im_size(2);




%%
%---------------------------------------------
% just least squares
%---------------------------------------------
if ~reg_TV && ~reg_l1 && ~reg_pos && ~reg_add
    U = pinv(M)*Y;
    res = norm(M*X-Y,'fro');
    return
end
%---------------------------------------------
% just ADDONE constrained (sum(x) = 1)
%---------------------------------------------
SMALL = 1e-12;
B = ones(1,n);
a = ones(1,N);

if  ~reg_TV && ~reg_l1 && ~reg_pos && reg_add
    F = M'*M;
    % test if F is invertible
    if rcond(F) > SMALL
        % compute the solution explicitly
        IF = inv(F);
        U = IF*M'*Y-IF*B'*inv(B*IF*B')*(B*IF*M'*Y-a);
        res = norm(M*U-Y,'fro');
        return
    end
    % if M'*M is singular, let sunsal_tv run
end


%%
%---------------------------------------------
%  Constants and initializations
%---------------------------------------------

% number of regularizers
n_reg =  reg_l1 + reg_pos + reg_add + reg_TV;

IF = inv(M'*M + n_reg*eye(n));

%%
%---------------------------------------------
%  Initializations
%---------------------------------------------

% no intial solution supplied
if U0 == 0
    U = IF*M'*Y;
end


index = 1

% initialize V variables
V = cell(1 + n_reg,1);

% initialize D variables (scaled Lagrange Multipliers)
D = cell(1 + n_reg,1);


%  data term (always present)
reg(1) = 1;             % regularizers
V{index} = M*U;         % V1
D{1} = zeros(size(Y));  % Lagrange multipliers

% next V
index = index + 1;
% POSITIVITY
if reg_pos == 1
    reg(index) = 2;
    V{index} = U;
    D{index} = zeros(size(U));
    index = index +1;
end
% ADDONE
if reg_add == 1
    reg(index) = 3;
    V{index} = U;
    D{index} = zeros(size(U));
    index = index +1;
end
%l_{1,1}
if reg_l1 == 1
    reg(index) = 4;
    V{index} = U;
    D{index} = zeros(size(U));
    index = index +1;
end

%%
%---------------------------------------------
%  AL iterations - main body
%---------------------------------------------
tol1 = sqrt(N)*1e-5;
i=1;
res = inf;
res2=inf;
AL_iters2=60;
nc=n_col;
nr=n_lin;
np=n;
k=1;
while (k <= AL_iters2)
    
NU = zeros(np,nc*nr);
X2=reshape((V{3}-D{3})',nc,nr,np);

for i_p = 1:np
    image_temp = zeros(nc+2,nr+2);
    NU_temp = image_temp;
    image_temp(2:end-1,2:end-1) = X2(:,:,i_p);
    for i1 = 2:nc+1
        for j1 = 2:nr+1
            NU_temp(i1,j1) = ((1/sqrt(2))*image_temp(i1-1,j1-1) + image_temp(i1-1,j1) + (1/sqrt(2))*image_temp(i1-1,j1+1)+...
                image_temp(i1,j1-1) + image_temp(i1,j1+1)+...
                (1/sqrt(2))*image_temp(i1+1,j1-1) + image_temp(i1+1,j1) + (1/sqrt(2))*image_temp(i1+1,j1+1))/(4*(1/sqrt(2))+4);
          
        end
    end   
    NU_p = NU_temp(2:end-1,2:end-1);
    NU(i_p,:) = NU_p(:);
end

w=1./(0.01+abs(NU));

NU2 = sqrt(sum((V{3}-D{3}).^2,2));
b=1./NU2;
a2=repmat(b,1,size(V{3},2));
w1=a2.*w;

while (i <= AL_iters) && (sum(abs(res)) > tol1)
    % solve the quadratic step (all terms depending on U)
    Xi = M'*(V{1}+D{1});
    for j = 2:(n_reg+1)
        Xi = Xi+ V{j} + D{j};
    end
    U = IF*Xi;
    
    % Compute the Mourau proximity operators
    for j=1:(n_reg+1)
        %  data term (V1)
        if  reg(j) == 1
            V{j} = (1/(1+mu)*(Y+mu*(M*U-D{j})));
        end
        %  positivity   (V2)
        if  reg(j) == 2
         V{j} = max(U-D{j},0);

        end
        % addone  (project on the affine space sum(x) = 1)  (V3)
        if  reg(j) == 3
            nu_aux = U - D{j};
            V{j} = nu_aux + repmat((1-sum(nu_aux))/n,n,1);
        end
        % l1 norm  (V4)
        if  reg(j) == 4
              V{j} = soft(U-D{j},lambda_l1/mu.*w1); 

        end

        
    end
    
    % update Lagrange multipliers
    
    for j=1:(n_reg+1)
        if  reg(j) == 1
            D{j} = D{j} - (M*U-V{j});
        else
            D{j} = D{j} - (U-V{j});
        end
    end
    
    % compute residuals
    if mod(i,10) == 1
        st = [];
        for j=1:(n_reg+1)
            if  reg(j) == 1
                res(j) = norm(M*U-V{j},'fro');
                st = strcat(st,sprintf(' res(%i) = %2.6f',reg(j),res(j) ));
            else
                res(j) = norm(U-V{j},'fro');
                st = strcat(st,sprintf('  res(%i) = %2.6f',reg(j),res(j) ));
            end
        end
        if  strcmp(verbose,'yes')
            fprintf(strcat(sprintf('iter = %i -',i),st,'\n'));
        end
    end
    
    
    
    % compute RMSE
    if true_x
        rmse(i)= norm(U-XT,'fro');
        if  strcmp(verbose,'yes')
            fprintf(strcat(sprintf('iter = %i - ||Xhat - X|| = %2.3f',i, rmse(i)),'\n'));
        end
        
    end
    
    i=i+1;
end    
i=1;
k=k+1;
res2(k)=res(1);
end








% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

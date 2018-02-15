function [t, q, qd, qdd, qddd] = trajectory_generation(nPoints, trajectoryType, varargin)

switch trajectoryType
    case 'Linear'
        
        %% Trajectory parameters
        
        % Start time
        ts = varargin{1};
        
        % Start position
        qs = varargin{2};
        
        % Finish time
        tf = varargin{3};
        
        % Finish position
        qf = varargin{4};
        
        %% Compute trajectory
        
        t = linspace(ts,tf,nPoints);
        slope = (qf-qs)/(tf-ts);
        q    = qs + slope*t;
        qd   = slope*ones(1,nPoints);
        qdd  = zeros(1,nPoints);
        qddd = zeros(1,nPoints);
        
    case 'Poly3'
        
        %% Trajectory parameters
        
        % Start time
        ts = varargin{1};
        
        % Start position
        qs = varargin{2};
        
        % Start velocity
        qsDot = varargin{3};
        
        % Finish time
        tf = varargin{4};
        
        % Finish position
        qf = varargin{5};
        
        % Finish velocity
        qfDot = varargin{6};
        
        % Via time
        tv = [];
        
        % Via position
        qv = [];
        
        % Number of vias
        nv = 0;
        
        %% Compute coefficients
        
        a = poly_3_coefficients(ts, qs, qsDot, tf, qf, qfDot, nv, tv, qv);
        
        %% Generate trajectory
        
        % Time vector
        t = linspace(ts,tf,nPoints);
        
        % Compute trajectory
        [q,qd,qdd,qddd] = poly_3(a(1),a(2),a(3),a(4),t);
        
    case 'Poly3WithVia'
        
        %% Trajectory parameters
        
        % Start time
        ts = varargin{1};
        
        % Start position
        qs = varargin{2};
        
        % Start velocity
        qsDot = varargin{3};
        
        % Finish time
        tf = varargin{4};
        
        % Finish position
        qf = varargin{5};
        
        % Finish velocity
        qfDot = varargin{6};
        
        % Via time
        tv = varargin{7};
        
        % Via position
        qv = varargin{8};
        
        % Number of vias
        nv = numel(qv);
        
        %% Compute coefficients
        
        a = poly_3_coefficients(ts, qs, qsDot, tf, qf, qfDot, nv, tv, qv);
        
        %% Generate trajectory
        
        % Number of points per segments
        nPointsPerSeg = round(nPoints/(nv+1));
        
        % Time vector
        t = zeros(nv+1,nPointsPerSeg);
        for iSeg = 1:(nv+1)
            if iSeg == 1
                t(iSeg,:) = linspace(ts,tv(iSeg),nPointsPerSeg);
            elseif iSeg == nv+1
                t(iSeg,:) = linspace(tv(iSeg-1),tf,nPointsPerSeg);
            else
                t(iSeg,:) = linspace(tv(iSeg-1),tv(iSeg),nPointsPerSeg);
            end
        end
        
        % Compute trajectory
        tTmp = [];
        q = [];
        qd = [];
        qdd = [];
        qddd = [];
        for iSeg = 1:(nv+1)
            [qTmp,qdTmp,qddTmp,qdddTmp] = poly_3(...
                a(1+4*(iSeg-1)), a(2+4*(iSeg-1)), ...
                a(3+4*(iSeg-1)), a(4+4*(iSeg-1)), t(iSeg,:));
            tTmp = [tTmp,t(iSeg,:)];
            q    = [q,qTmp];
            qd   = [qd,qdTmp];
            qdd  = [qdd,qddTmp];
            qddd = [qddd,qdddTmp];
        end
        t = tTmp;
        
        % Remove duplicate point
        [t,iUnique,~] = unique(t);
        q = q(iUnique);
        qd = qd(iUnique);
        qdd = qdd(iUnique);
        qddd = qddd(iUnique);
        for iPoint = 1:(nPoints-numel(t))
            t(end+1)    = t(end);
            t(end-1)    = (t(end-2) + t(end-1))/2;
            q(end+1)    = q(end);
            q(end-1)    = (q(end-2) + q(end-1))/2;
            qd(end+1)   = qd(end);
            qd(end-1)   = (qd(end-2) + qd(end-1))/2;
            qdd(end+1)  = qdd(end);
            qdd(end-1)  = (qdd(end-2) + qdd(end-1))/2;
            qddd(end+1) = qddd(end);
            qddd(end-1) = (qddd(end-2) + qddd(end-1))/2;
        end
        
    case 'Poly5'
        
        %%
        %
        
    case 'Poly5WithVia'
        
        %%
        %
        
    case 'LinearWithParabolicBlends'
        
        %%
        %
        
end

end

%% poly_3_coefficients()
function a = poly_3_coefficients(ts, qs, qsDot, tf, qf, qfDot, nv, tv, qv)

% Intialize A and b
A = zeros(4+4*nv,4+4*nv);
b = zeros(4+4*nv,1);

% Populate A and b
A(1,1:4) = [1,ts,ts^2,ts^3];
A(2,1:4) = [0,1,2*ts,3*ts^2];
b(1) = qs;
b(2) = qsDot;
for i = 1:nv
    A(3+4*(i-1),4*(i-1)+1:4*i)     = [1,tv(i),tv(i)^2,tv(i)^3];
    A(4+4*(i-1),4*(i-1)+1:4*(i+1)) = [0,1,2*tv(i),3*tv(i)^2,0,-1,-2*tv(i),-3*tv(i)^2];
    A(5+4*(i-1),4*(i-1)+1:4*(i+1)) = [0,0,2,6*tv(i),0,0,-2,-6*tv(i)];
    A(6+4*(i-1),4*i+1:4*(i+1))     = [1,tv(i),tv(i)^2,tv(i)^3];
    b(3+4*(i-1)) = qv(i);
    b(4+4*(i-1)) = 0;
    b(5+4*(i-1)) = 0;
    b(6+4*(i-1)) = qv(i);
end
A(7+4*(nv-1),4*nv+1:4*(nv+1)) = [1,tf,tf^2,tf^3];
A(8+4*(nv-1),4*nv+1:4*(nv+1)) = [0,1,2*tf,3*tf^2];
b(7+4*(nv-1)) = qf;
b(8+4*(nv-1)) = qfDot;

% Solve a = A\b
a = A\b;

end

%% poly_3()
function [q,qd,qdd,qddd] = poly_3(a0,a1,a2,a3,t)
q    = a0 + a1.*t + a2.*t.^2 + a3.*t.^3;
qd   =      a1 +  2*a2.*t +  3*a3.*t.^2;
qdd  =            2*a2 +     6*a3.*t;
qddd =                       6*a3*ones(1,numel(t));
end

%% poly_5()
function [q,qd,qdd,qddd] = poly_5(a0,a1,a2,a3,a4,a5,t)
q    = a0 + a1.*t + a2.*t.^2 + a3.*t.^3 +   a4.*t.^4 +    a5.*t.^5;
qd   =      a1 +  2*a2.*t +  3*a3.*t.^2 + 4*a4.*t.^3 +  5*a5.*t.^4;
qdd  =            2*a2 +     6*a3.*t +   12*a4.*t.^2 + 20*a5.*t.^3;
qddd =                       6*a3 +      24*a4.*t +    60*a5.*t.^2;
end
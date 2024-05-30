function varargout = dome(flag, opt)
% pop-up dome generator
% https://www.youtube.com/watch?v=2STS0POwB7g

    arguments
        flag char = ''
        opt.wedges (1,1) double = 8
        opt.segments (1,1) double = nan
    end

    % return cell array of local functions, for unit testing
    if strcmp(flag, 'localfunctions')
        varargout{1} = localfunctions;
        return
    end

    if isnan(opt.segments), opt.segments = 2 + 4*round((opt.wedges - 2)/4); end
    [wedges, segments] = deal(opt.wedges, opt.segments);

    validateattributes(wedges,'numeric',{'integer','even','positive','>=',4})
    validateattributes(segments,'numeric',{'integer','even','positive'})
    assert(mod(segments-2, 4) == 0, 'Segments must be 2 + 4^n for n > 0, e.g. 6, 10, 14, ...')

    [azimuth, elevation, V] = vertices3d(wedges, segments);

    F = arrayfun(@(n) wedge_faces(n, segments, wedges), 1:wedges, UniformOutput=false);
    F = cat(1,F{:});

    figure(1); clf(); hold on;
    if contains(flag, 'debug')
        plot3(V(:,1),V(:,2),V(:,3),'go');
        text(V(:,1),V(:,2),V(:,3),string(1:size(V,1)))
    end
    patch('Faces',F,'Vertices', V, FaceColor='r', FaceAlpha=0.2)
    view(30,60)
    axis equal

    figure(2); clf(); hold on; axis equal
    % unit circle, and top projection
    ref_circles(1, 3, 3, 'color', 'b', 'LineWidth', 0.1)
    patch('Faces',F, 'Vertices', V + [3,3,0], FaceColor='none', EdgeColor='b')
    
    for j = 1:wedges/2
        [v,f] = flat_wedge(V, segments);
        r = vecnorm(v((1:segments/2)+2,:)');
        if j ==1
            ref_circles(r, 0, 0, 'color', [1,1,1]*0.8, 'LineWidth', 0.1)
        end
        R = rotmat(2*j*360/wedges - 180/wedges)';
        v = v*R;
        ref_circles(r, v(2,1), v(2,2), 'color', [1,1,1]*0.8, 'LineWidth', 0.1)
        [v,f] = wedge_flaps(v,f);
        
        patch('Faces',f,'Vertices', v, FaceColor='none')
    end

    annotations('color', 'b', 'LineWidth', 0.1)
end


function [azimuth, elevation, V] = vertices3d(wedges, segments)
% returns 2 + (S-1)*W vertices, arranged to work with WEDGE_FACES

    azimuth = linspace(-180, 180-360/wedges, wedges);
    elevation = linspace(-90 + 180/segments, 90 - 180/segments, segments/2);

    [az, el] = meshgrid(azimuth, elevation);
    [x,y,z] = sph2cart(az(:)*pi/180, el(:)*pi/180, 1);

    V = [[0,0,z(1)];
        [0,0,z(end)];
        [x,y,z]];

    % scale XY such that edges are tangent to unit sphere
    V(:,1:2) = V(:,1:2)/cosd(180/wedges);
end

function F = wedge_faces(offset, segments, wedges)
% return an S x 4 matrix of face-vertex indices for wedge w

    %     1+e ... 2e            3+e ... 2e+2
    %  -1            0  -->   1              2  
    %      1  ... e              3  ... e+2

    e = segments/2;
    ii = 1:e-1;
    s0 = (offset-1)*e;
    F = [[s0 + e + 1, -1, s0+1, NaN];
        [ii; ii+1; ii+e+1; ii+e]' + s0;
        [s0 + e, 0, s0 + 2*e, NaN]];
    
    if nargin > 2
        M = segments/2*wedges;
        F = rem(F-1,M) + 3;
    else
        F = F + 2;
    end
end

function [v,f] = flat_wedge(V, S)
% return a 2D "unfolded" projection of one wedge (pointing along x axis)

    f = wedge_faces(1, S);

    ii = (1:S/2) + 2;
    jj = ii + S/2;
    v = V([1,2,ii,jj],:);

    % points along center axis
    a = [v(1,:); (v(ii,:) + v(jj,:))/2; v(2,:)];

    xa = cumsum(rssq(diff(a,1,1),2));
    d = rssq(v(ii,:) - v(jj,:),2);

    v(ii,1) = xa(1:end-1);
    v(2,1) = xa(end);
    v(ii,2) = -d/2;
    v(jj,1) = v(ii,1);
    v(jj,2) = d/2;
    v(:,3) = [];
end

function [v,f] = wedge_flaps(v,f)
% add alternate "flaps" (neighboring faces) to wedge

    validateattributes(v, 'numeric', {'size', [NaN,2]})
    nans = isnan(f);
    f(nans) = 1;
    validateattributes(f, 'numeric', {'integer', 'positive', 'size', [NaN,4], '<=', size(v,1)})
    f(nans) = NaN;

    n = size(f,1);
    for j = 1:n
        if any(isnan(f(j,:)))
        % don't add triangular faces
            continue
        end
        if mod(j,2) == 0
            p = f(j,1);
            q = f(j,2);
            r = f(j,[3,4]);
        else
            p = f(j,3);
            q = f(j,4);
            r = f(j,[1,2]);
        end
        vr = reflect(v(r,:), v(p,:), v(q,:));
        idx = size(v,1) + (1:2);
        v(idx,:) = vr;
        f(end+1,:) = [p, q, idx]; %#ok<AGROW>
    end
end

function R = rotmat(x)
    R = [cosd(x), -sind(x); sind(x), cosd(x)];
end

function Q = reflection_matrix(x, y)
    Q = [(x+y)*(x-y), 2*x*y; 2*x*y, (x+y)*(y-x)]/(x^2 + y^2);
end

function Vr = reflect(V, p, q)
% reflect vertices V along line p-q
    Q = reflection_matrix(q(1)-p(1), q(2)-p(2));
    Vr = (V-p)*Q + p;
end

function ref_circles(r, x0, y0, varargin)

    x = cosd(0:360).*r(:) + x0;
    y = sind(0:360).*r(:) + y0;
    plot(x', y', varargin{:})
end

function fixedangle = stdangle(angle)
% Limit any angle(s) to values between -180° and 180°
    fixedangle = rem(rem(angle,360) + 540,360)-180;
end

function annotations(varargin)

    % internal tensors
    t = 0.2;
    a = (pi^2-8)/(4*pi - 8);
    x = [-t, 0, a, pi/2, pi/2 + t] + 2;
    y = -(0:4)*t - 2;
    
    arrayfun(@(x) plot([x,x], y([1,end]), varargin{:}), x)
    arrayfun(@(y) plot(x([1,end]), [y,y], varargin{:}), y)
end
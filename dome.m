function varargout = dome(wedges, segments, polygons, flag, opt)
% DOME() - pop-up dome generator
% DOME(wedges, segments) - Specify number of WEDGES and SEGMENTS
% DOME(wedges, segments, polygons) - Add projected POLYGONS: cell array of 
%   [az, el] (degrees) matrices, or 'earth' for world map.
% DOME(..., Name, Value)
%   'segmentdist' - {}'parallels' (default), 'regular', 'maxvolume'}, 
%   how to distribute segments.
%   'radius' - 1 (default), dome radius.
%   'flaps' - true (default), add alternate tabs cut-out.
%   'inset' - 1e-2 (default), inset for flaps.
%   'figsave' - false (default), save figures as PDF.
% DOME(..., 'debug') - Show debug information
% H = DOME([], [], [], 'localfunctions') - Return local function handles
% 
% REF: https://www.youtube.com/watch?v=2STS0POwB7g

    arguments
        wedges double = 12
        segments double = 2 + 4*round((wedges - 2)/4)
        polygons = 'earth'
        flag char = ''
        opt.segmentdist = 'parallels';
        opt.radius = 1;
        opt.flaps (1,1) logical = true
        opt.inset (1,1) double = 1e-2
        opt.figname = sprintf('%dx%d_',wedges, segments);
        opt.figsave = false;
    end

    % return cell array of local functions, for unit testing
    if strcmp(flag, 'localfunctions')
        handles = localfunctions;
        names = cellfun(@func2str, handles, 'UniformOutput', false);
        varargout{1} = cell2struct(handles, names, 1);
        return
    end
    opt.debug = contains(flag, 'debug');
    if ~opt.debug && ~isempty(flag)
        error('Unrecognized flag')
    end

    validateattributes(wedges,'numeric',{'scalar','integer','even','positive','>=', 4})
    validateattributes(segments,'numeric',{'scalar','integer','even','positive'})
    assert(mod(segments-2, 4) == 0, 'Segments must be 2 + 4^n for n > 0, e.g. 6, 10, 14, ...')

    if ischar(polygons) && strcmpi(polygons, 'earth')
        polygons = load('earth.mat', 'landpolygons').landpolygons;
    end

    [azimuth, elevation, V] = vertices(wedges, segments, opt);
    F = polyhedron_faces(segments, wedges);
    
    figure([opt.figname, 'polyhedron']); clf(); hold on; axis off;
    plot_polyhedron(V, F, opt.debug);

    for hemisphere = ['N', 'S']
        figure([opt.figname, hemisphere]); clf(); hold on; axis equal; axis off;

        % unit circle, and top projection
        plot_views(F, V, opt);
    
        if ~isempty(polygons)
            plot_polygons(polygons, azimuth, elevation, hemisphere, opt);
        end
        
        plot_cutout(V, F, azimuth, elevation, hemisphere, opt);

        % internal_structure()
        if opt.figsave
            set(gcf, 'PaperOrientation', 'landscape', 'Renderer', 'painters');
            set(gca,'Position', [0 0 1 1]);
            set(gcf,'PaperUnits', 'normalized', 'PaperPosition', [0 0 1 1]);
            saveas(gcf,[gcf().Name '.pdf']);
        end
    end

end

function args = styles(what)
    switch what
        case 'polyhedron'
            args = {'FaceColor', 'r', 'FaceAlpha', 0.2};
        case 'polyshape'
            % args = {'EdgeColor', 'k', 'FaceColor', 'none', 'LineWidth', 0.1};
            args = {'EdgeColor', [1,1,1]*0.8, 'FaceColor', 'none', 'LineWidth', 0.1};
        case 'polygon'
            args = {'color', [0.1,0.6,0.8] , 'LineWidth', 0.1};
        case 'annotation'
            args = {'EdgeColor', 'b', 'FaceColor', 'none', 'LineWidth', 0.1};
        case 'refline'
            % args = {'color', [1,1,1]*0.8, 'LineWidth', 0.1};
            args = {'color', 'none', 'LineWidth', 0.1};
        otherwise
            error('Unknown style: %s', what)
    end
end

function elevation_breaks = max_volume(segments)
% Return elevation breakpoints that maximize the frustum volume

    n = (segments - 2)/4;
    A = eye(n-1, n) - circshift(eye(n-1, n), 1, 2);
    b = zeros(n-1, 1);
    x0 = (1:n)'*360/segments;
    x = fmincon(@vol, x0, A, b, [], [], [], [], [], optimoptions('fmincon','Display','none'));
    elevation_breaks = [-flipud(x); 0; x]';
    
    function v = vol(breaks)
        c = [1; cosd(breaks)];
        s = [0; sind(breaks)];

        R = c(1:end-1);
        r = c(2:end);
        h = diff(s);

        v = 2/3*pi - sum(pi/3*h.*(R.^2 + R.*r + r.^2));
    end
end

function [azimuth, elevation, V] = vertices(wedges, segments, opt)
% returns 2 + (S-1)*W vertices, arranged to work with WEDGE_FACES

    azimuth = linspace(-180, 180-360/wedges, wedges);
    switch opt.segmentdist
        case 'maxvolume'
            elevation = max_volume(segments);
        case 'parallels'
            elevation = linspace(-90 + 360/(segments+2), 90 - 360/(segments+2), segments/2);
        case 'regular'
            elevation = linspace(-90 + 180/segments, 90 - 180/segments, segments/2);
    end

    [az, el] = meshgrid(azimuth, elevation);
    [x,y,z] = sph2cart(az(:)*pi/180, el(:)*pi/180, 1);

    V = [[0,0,z(1)];
        [0,0,z(end)];
        [x,y,z]]*opt.radius;

    % radius XY such that edges are tangent to unit sphere
    V(:,1:2) = V(:,1:2)/cosd(180/wedges);
end

function plot_polygons(polygons, azimuth_breaks, elevation_breaks, hemisphere, opt)

    validateattributes(polygons, {'cell'}, {'vector'})
    assert(all(cellfun(@(x) isnumeric(x) && ismatrix(x) && size(x,2) == 2, polygons)))

    P = cell(size(polygons));
    for j = 1:numel(polygons)
        [az, el] = deal(polygons{j}(:,1), polygons{j}(:,2));
        az = stdangle(az);
        if ~(az(end) == az(1) && el(end) == el(1))
            az(end+1) = az(1); %#ok<AGROW>
            el(end+1) = el(1); %#ok<AGROW>
        end
        P{j} = project(azimuth_breaks, elevation_breaks, az, el, hemisphere);
    end
    P = cat(1,P{:})*opt.radius;

    style = styles('polygon');
    plot(P(:,1), P(:,2), style{:})
end

function [edges, centers] = wedge_centers(azimuth_breaks, hemisphere)

    d = unique(diff(azimuth_breaks));
    assert(isscalar(d) && d > 0);

    switch upper(hemisphere)
        case 'N'
            centers = azimuth_breaks + d/2;
        case 'S'
            centers = azimuth_breaks - d/2;
    end
    edges = centers(1:2:end);
    centers = centers(2:2:end);
end

function G = project(azimuth_breaks, elevation_breaks, az, el, hemisphere)
% Rotated Cassini-like "Prismatic-Equidistant" projection, for each wedge

    validateattributes(azimuth_breaks,'numeric',{'vector','<=', 180, '>=', -180})
    validateattributes(elevation_breaks,'numeric',{'vector','<=', 90, '>=', -90})
    validateattributes(az,'numeric',{'vector','<=', 180, '>=', -180})
    validateattributes(el,'numeric',{'vector','<=', 90, '>=', -90})

    azimuth_breaks = azimuth_breaks(:);

    [edges, centers] = wedge_centers(azimuth_breaks, hemisphere);
    switch upper(hemisphere)
        case 'N'
            x = polygon_arc(elevation_breaks, -el); 
        case 'S'
            x = polygon_arc(elevation_breaks, el); 
    end

    bin = discretize(az, edges);
    bin(az < edges(1) | az > edges(end)) = length(edges);
    y = tand(az - centers(bin)) .* cosd(el);

    c = cosd(centers);
    s = sind(centers);

    P = [x.*c(bin) - y.*s(bin), x.*s(bin) + y.*c(bin)];

    % insert nans in bin jumps
    idx = (1:size(P,1))' + cumsum(bin ~= circshift(bin,1));
    G = nan(max(idx)+1, 2);
    G(idx,:) = P;
end

function x = polygon_arc(vertex_angles, query_angles)
% Arc-length along a unit-circle-inscribed polygon with vertices at VERTEX_ANGLES
% for QUERY_ANGLES, in measured from -90 degrees

    validateattributes(vertex_angles,'numeric',{'size', [1,NaN],'<', 90, '>', -90})
    validateattributes(query_angles,'numeric',{'vector','<=', 90, '>=', -90})
    assert(any(vertex_angles < 0) && any(vertex_angles > 0));

    vertex_angles = [nan, unique(vertex_angles), nan];
    vertex_angles(1) = stdangle(180 - vertex_angles(2));
    vertex_angles(end) = stdangle(180 - vertex_angles(end-1));

    side_angles = diff(vertex_angles);
    full_sides = 2*sind(side_angles/2);
    full_sides = [0, cumsum(full_sides)] - full_sides(1)/2;

    if size(query_angles, 2) == 1
        vertex_angles = vertex_angles';
        side_angles = side_angles';
        full_sides = full_sides';
    end

    bin = discretize(query_angles, vertex_angles);
    delta = query_angles - vertex_angles(bin);
    x = sind(delta)./cosd(side_angles(bin)/2 - delta);
    x = full_sides(bin) + x;
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

function F = polyhedron_faces(segments, wedges)

    F = cell(wedges,1);
    for j = 1:wedges
        f = wedge_faces(j, segments, wedges);
        F{j+1} = f(2:end-1,:);
    end
    F = cat(1,F{:});
    F(:,end+1:wedges) = NaN;
    F(end+1:end+2,:) = (0:wedges-1)*segments/2 + [3; 2+segments/2];
end

function [v, f] = flat_wedge(V, segments, opt)
% return a 2D "unfolded" projection of one wedge (pointing along x axis)

    f = wedge_faces(1, segments);

    ii = (1:segments/2) + 2;
    jj = ii + segments/2;
    v = V(1:segments+2,:);

    % points along center axis
    a = [v(1,:); (v(ii,:) + v(jj,:))/2; v(2,:)];

    xa = cumsum(rssq(diff(a,1,1),2));
    d = rssq(v(ii,:) - v(jj,:),2);

    % "unfold" wedge (should be equivalent to projection)
    v(ii,1) = xa(1:end-1);
    v(2,1) = xa(end);
    v(ii,2) = -d/2;
    v(jj,1) = v(ii,1);
    v(jj,2) = d/2;
    v(:,3) = [];

    % outline = [1, ii, 2, fliplr(jj), 1];
    % creases = reshape([ii; jj; nan(1,segments/2)], 1, []);

    % add alternate "flaps" (neighboring faces) to wedge
    if ~opt.flaps, return; end

    n = size(f,1);
    for j = 1:n
        if any(isnan(f(j,:)))
        % don't add triangular faces
            continue
        end
        if mod(j,2) == 0
            p = f(j,1);
            q = f(j,2);
            r = f(j,3);
            s = f(j,4);
        else
            p = f(j,3);
            q = f(j,4);
            r = f(j,1);
            s = f(j,2);
        end
        if j > (segments-2)/4 + 1
            c = v(2,:);
        else
            c = v(1,:);
        end
        vr = flap(v(p,:), v(q,:), v(r,:), v(s,:), c, opt.inset);
        if opt.inset > 0
            idx = size(v,1) + (1:4);
            f(end+1,:) = idx; %#ok<AGROW>
        else
            idx = size(v,1) + (1:2);
            f(end+1,:) = [p, idx, q]; %#ok<AGROW>
        end
        v(idx,:) = vr;
    end
end

function v = flap(p, q, r, s, c, inset)
% Returns vertices [p'] s' r' [q'], where:
%   s' = Q*s is s reflected around p-q
%   r' is r reflected around c-q, and extended to the line s'Q*r,
%      this is to make sure the flap lies within the flat circle centered at c
%   If inset > 0, it is applied to p', q' along pq, and to s', r' along s'r'.
%
%      p  --- q                s' -- r'      Qr
%       \     |               /       \  
%   c    \    |              /         \
%        s -- r         p - p'          q' - q

    sm = reflect(s, p, q);
    rm = reflect(r, p, q);

    if norm(c-p) < norm(c-q)
        x = reflect(r, c, q);
        rm = intersect(rm, sm, x, q);
    else
        x = reflect(s, p, c);
        sm = intersect(rm, sm, x, p);
    end
    [r, s] = deal(rm, sm);

    pq = (q - p)/norm(q-p);
    sr = (r - s)/norm(r-s);
    assert(sr*pq' > 0);

    if inset > 0
        assert(inset < norm(q-p));
        v = [p + pq*inset;
             s + sr*inset;
             r - sr*inset;
             q - pq*inset];
    else
        v = [s; r];
    end
end

function R = rotmat(x)
    R = [cosd(x), -sind(x); sind(x), cosd(x)];
end

function Q = reflection_matrix(x, y)
    Q = [(x+y)*(x-y), 2*x*y; 2*x*y, (x+y)*(y-x)]/(x^2 + y^2);
end

function x = intersect(a,b,c,d)
    ab = (b-a)/norm(b-a);
    cd = (d-c)/norm(d-c);
    t = [ab',-cd']\(c-a)';
    x = (t'*[ab;cd] + a + c)/2;
end

function Vr = reflect(V, p, q)
% reflect vertices V along line p-q
    Q = reflection_matrix(q(1)-p(1), q(2)-p(2));
    Vr = (V-p)*Q + p;
end

function ref_circles(r, x0, y0)

    style = styles('refline');
    x = cosd(0:360).*r(:) + x0;
    y = sind(0:360).*r(:) + y0;
    plot(x', y', style{:})
end

function fixedangle = stdangle(angle)
% Limit any angle(s) to values between -180° and 180°
    fixedangle = rem(rem(angle,360) + 540,360)-180;
end

function varargout = figure(name)
    h = findobj('type','figure','name',name);
    if isempty(h) || ~ishandle(h)
        builtin('figure', 'name', name)
    else
        builtin('figure', h)
    end
    if nargout > 0, varargout{1} = gcf(); end
end

function plot_polyhedron(V, F, debug)

    if debug
        plot3(V(:,1),V(:,2),V(:,3),'go');
        text(V(:,1),V(:,2),V(:,3),string(1:size(V,1)))
    end

    args = styles('polyhedron');
    patch('Faces',F,'Vertices', V, args{:})

    view(30,60)
    axis equal
end

function plot_cutout(V, F, azimuth, elevation, hemisphere, opt)

    style = styles('polyshape');
    segments = 2*numel(elevation);
    wedges = numel(azimuth);

    [~, centers] = wedge_centers(azimuth, hemisphere);

    r = polygon_arc(elevation, elevation);
    ref_circles(r, 0, 0);

    patch('Faces',F(end,:),'Vertices', V, style{:})

    for j = 1:wedges/2
        [v,f] = flat_wedge(V, segments, opt);

        R = rotmat(centers(j))';
        v = v*R;    
        patch('Faces',f(2:end,:),'Vertices', v, style{:})
    end

    if opt.debug
        [az, el] = meshgrid(azimuth, elevation);
        P = project(azimuth, elevation, az(:), el(:), hemisphere);
        plot(P(:,1), P(:,2), 'go');
    end
end

function plot_views(F, V, opt)
    style = styles('annotation');
    ref_circles(1, 4*opt.radius, 2*opt.radius)
    patch('Faces',F, 'Vertices', V + [4,2,0]*opt.radius, style{:})

    R = circshift(eye(3),-1,2);
    ref_circles(1, 4*opt.radius, -2*opt.radius)
    patch('Faces',F, 'Vertices', V*R + [4,-2,0]*opt.radius, style{:})
end
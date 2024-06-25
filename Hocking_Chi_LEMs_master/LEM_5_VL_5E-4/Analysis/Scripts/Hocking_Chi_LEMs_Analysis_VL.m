%% Companion script to analyze Hocking_Chi_LEMs.
%  Author: Chris Sheehan

% Reset script
clear all;
close all;
clc;

%% Handle directories

% Print status
disp('Handling directories...');

% Get cd and parent
script_directory = cd;
analysis_directory = erase(script_directory, '\Scripts');
model_directory = erase(analysis_directory, '\Analysis');
input_directory = [analysis_directory '\Input'];

% Import input parameters
Params = readtable([input_directory '\Parameters.csv']);

% Create output directories
mkdir([analysis_directory '\Output']);
mkdir([analysis_directory '\Output\Videos']);
mkdir([analysis_directory '\Output\Timeseries']);

% Get list of DEMs
contents = dir([model_directory '\Output\Main_run\topographic__elevation__tif']);
names = {contents.name};
List = {};
i = 1;
for c = 1 : size(names, 2)
    if contains(names{c}, '.tif') == 1 
        if contains(names{c}, '.aux') == 0 
            List{i} = names{c};
            i = i + 1;
        end
    end
end
List = erase(List, '.tif');
List = str2double(List);
List = sort(List);
ints = List / Params.plot_interval;
ints = rem(ints, 1) == 0;
Print_List = List(ints);
List_nums = List;
List = string(List);
Print_List = string(Print_List);
List = append(List, '.tif');
Print_List = append(Print_List, '.tif');

% Inport model parameters
spin_params = readtable([model_directory '\Input\Spin_up.csv']);
main_params = readtable([model_directory '\Input\Main_run.csv']);

% Import model uplift rate
u = readtable([model_directory '\Input\Spin_up.csv']);
u = u.u;

% Create colormap
cm = flipud(ttscm('roma'));  % load roma color map and reverse it

% Toggle off figure visibility
set(0, 'DefaultFigureVisible', 'off');

% Print space
disp(' ');

%% Start DEM loop

for t = 1 : size(List, 2)

%% Create directories

% Print status
disp(['DEM = ' List{t}]);

% Create directories
mkdir([analysis_directory '\Output\' erase(List{t}, '.tif') '/Images']);
if contains(List{t}, Print_List) == 1
    mkdir([analysis_directory '\Output\' erase(List{t}, '.tif') '/Grids']);
    mkdir([analysis_directory '\Output\' erase(List{t}, '.tif') '/Shapefiles']);
end
out_directory = [analysis_directory '\Output\' erase(List{t}, '.tif')];

% Print space
disp(' ');

%% Import and condition DEM

% Print status
disp('Importing and conditioning DEM...');

% DEM / dzdt Path
DEM_path = [model_directory '\Output\Main_run\topographic__elevation__tif\' List{t}];
dzdt_path = [model_directory '\Output\Main_run\dzdt__tif\' List{t}];

% Load DEM / dzdt
DEM_original = GRIDobj(DEM_path);
dzdt_grid = GRIDobj(dzdt_path);

% Block boundaries by increasing boundary elevations to the maximum
% elevation in the domain
DEM_original.Z(:, end) = max(max(DEM_original.Z));      % Eastern boundary
DEM_original.Z(1, :) = max(max(DEM_original.Z));        % Northern boundary
DEM_original.Z(:, 1) = max(max(DEM_original.Z));        % Western boundary
DEM_original.Z(end, :) = max(max(DEM_original.Z));      % Southern boundary

% Find breach location and enforce its elevation as 0
breach_r = size(DEM_original.Z, 1);
breach_c = main_params.breach_x;
DEM_original.Z(breach_r, breach_c) = 0;

% Preprocess DEM
if matches(Params.preprocess_method , 'carve') == true
    Flow_Direction = FLOWobj(DEM_original, 'preprocess', 'carve');
    DEM = imposemin(Flow_Direction, DEM_original); 
elseif matches(Params.preprocess_method , 'fill') == true
    DEM = fillsinks(DEM_original);
    Flow_Direction = FLOWobj(DEM,'preprocess', 'none');
elseif matches(Params.preprocess_method , 'none') == true
    DEM = DEM_original;
end

% Create DEM Colormap
zlimits = [min(DEM.Z(:)) max(DEM.Z(:))];
cc = demcmap(zlimits);

% Record Preproccessing Energy
preprocessing_change = DEM.Z;
preprocessing_change = DEM.Z - DEM_original.Z;
preprocessing_change(find(isnan(preprocessing_change))) = 0;
preprocessing_energy = sum(preprocessing_change(:)) * (DEM.cellsize)^2;
pc_min = abs(min(min(preprocessing_change)));
pc_max = abs(max(max(preprocessing_change)));
preprocessing_scale = [-max(pc_min, pc_max) max(pc_min, pc_max)];
preprocessing_max_change = max(pc_min, pc_max);
disp(['Preprocessing_energy = ' num2str(preprocessing_energy) ' m^3']);
disp(['Preprocessing_Max_Change = ' num2str(preprocessing_max_change) ' m']);

% Print space
disp(' ');

%% Create streams and divides

% Print status
disp('Creating streams and divides...');

% Streams
Streams = STREAMobj(Flow_Direction, 'minarea', Params.minarea, 'unit', 'map');

% Divides
lastwarn('');
Divides = DIVIDEobj(Flow_Direction, Streams, 'type', Params.divide_type, 'verbose', false);
[warnMsg, warnId] = lastwarn;
if contains(warnMsg, 'divide segments') == 0
    Divides_Lines = divorder(Divides);
end

% Print space
disp(' ');

%% Calculate grid morphometrics (during specified plot interval)

% Execute at specified interval
if contains(List{t}, Print_List) == 1

    % Print status
    disp('Calculating morphometrics...')
    
    % Morphometrics
    Slope_Radians = gradient8(DEM);
    Slope_Degrees = gradient8(DEM,'degree');  
    Relief = localtopography(DEM, Params.relief_radius);
    Aspect = aspect(DEM);
    hmean = localtopography(DEM, Params.relief_radius, 'type', 'mean');
    hmin =  localtopography(DEM, Params.relief_radius, 'type', 'min');
    hmax =  localtopography(DEM, Params.relief_radius, 'type', 'max');
    Hypsometric_Index = (hmean - hmin) / (hmax - hmin);
    
    % Print space
    disp(' ');

end

%% Calculate flow metrics

% Print status
disp('Calculating flow metrics')

% Flow metrics
Drainage_Cells  = flowacc(Flow_Direction);
Drainage_Area  = flowacc(Flow_Direction).*(Flow_Direction.cellsize^2);

% Print space
disp(' ');

%% Export morphometric images (during specified plot interval)

% Execute at specified interval
if contains(List{t}, Print_List) == 1

    % Print status
    disp('Exporting morphometric images...')
    
    % DEM_Map
    fig = figure();
    hold on;
    imageschs(DEM, [], 'colormap', cc, 'colorbarylabel', 'Elevation (m)');
    title(['Year = ' erase(List{t}, '.tif')]);
    exportgraphics(fig, [out_directory '\Images\DEM_Map.' char(Params.file_type)], 'Resolution', Params.file_dpi);
    close fig 1
    
    % Slope (degrees)
    fig = figure();
    hold on;
    imageschs(DEM, Slope_Degrees, 'colorbarylabel', 'Slope (degrees)');
    title(['Year = ' erase(List{t}, '.tif')]);
    exportgraphics(fig, [out_directory '\Images\Slope_Degrees_Map.' char(Params.file_type)], 'Resolution', Params.file_dpi);
    close fig 1
    
    % Relief
    fig = figure();
    hold on;
    imageschs(DEM, Relief, 'colorbarylabel', 'Relief (meters)');
    title(['Year = ' erase(List{t}, '.tif')]);
    exportgraphics(fig, [out_directory '\Images\Relief_Map.' char(Params.file_type)], 'Resolution', Params.file_dpi);
    close fig 1
    
    % Drainage Area
    fig = figure();
    hold on;
    imageschs(DEM, log10(Drainage_Area), 'colorbarylabel', 'log10(Drainage Area)');
    title(['Year = ' erase(List{t}, '.tif')]);
    exportgraphics(fig, [out_directory '\Images\Drainage_Area_Map.' char(Params.file_type)], 'Resolution', Params.file_dpi);
    close fig 1
    
    % Stream Map
    fig = figure();
    hold on;
    imageschs(DEM, [], 'colormap', cc, 'colorbarylabel', 'Elevation (m)');
    plot(Streams, 'b-')
    title(['Year = ' erase(List{t}, '.tif')]);
    exportgraphics(fig, [out_directory '\Images\Stream_Map.' char(Params.file_type)], 'Resolution', Params.file_dpi);
    close fig 1
    
    % Divides_Lines_Map
    fig = figure();
    hold on;
    imageschs(DEM, [], 'colormap', cc, 'colorbarylabel', 'Elevation (m)');
    plot(Divides_Lines, 'color','k');
    title(['Year = ' erase(List{t}, '.tif')]);
    exportgraphics(fig, [out_directory '\Images\Divide_Lines_Map.' char(Params.file_type)], 'Resolution', Params.file_dpi);
    close fig 1
    
    % Preprocessing_Change_Map
    fig = figure();
    hold on;
    imageschs(DEM, preprocessing_change, 'colormap', colormap(flipud(redblue)), 'caxis', preprocessing_scale, 'colorbarylabel', 'Elevation (m)');
    title(['Year = ' erase(List{t}, '.tif')]);
    exportgraphics(fig, [out_directory '\Images\Preprocessing_Change_Map.' char(Params.file_type)], 'Resolution', Params.file_dpi);
    close fig 1
    
    % Hypsometric_Index_Map
    fig = figure();
    hold on;
    imageschs(DEM, Hypsometric_Index, 'colorbarylabel', 'Hypsometric Index');
    title(['Year = ' erase(List{t}, '.tif')]);
    exportgraphics(fig, [out_directory '\Images\Hypsometric_Index_Map.' char(Params.file_type)], 'Resolution', Params.file_dpi);
    close fig 1
    
    % Print space
    disp(' ');

end

%% Export shapefiles (during specified plot interval)

% Execute at specified interval
if contains(List{t}, Print_List) == 1

    % Print status
    disp('Exporting stream and divide shapefiles...')
    
    % Streams
    MS = STREAMobj2mapstruct(Streams);
    shapewrite(MS, [out_directory '\Shapefiles\Streams.shp']);
    
    % Divides
    MS = DIVIDEobj2mapstruct(Divides_Lines, DEM, 1000);
    shapewrite(MS, [out_directory '\Shapefiles\Divides.shp']);
    
    % Print space
    disp(' ');

end

%% Channel metrics

% Print status
disp('Calculating channel metrics...')

% Smooth Streams
Z_smooth = smooth(Streams, DEM, 'K', Params.smoothing_factor);

% Concavity
if char(Params.calculate_theta) == true
    [mn, results] = mnoptim(Streams, DEM, Drainage_Area, 'crossval', false);
    theta_calc = mn.(1);
    close fig 1;
    close fig 2;
end

% Ksn
if char(Params.calculate_theta) == true
    Ksn = ksn(Streams, DEM, Drainage_Area, theta_calc, Params.smoothing_factor);
else
    Ksn = ksn(Streams, DEM, Drainage_Area, Params.theta, Params.smoothing_factor);
end
KsnMapped = mapfromnal(Flow_Direction, Streams, Ksn);

% Chi
if char(Params.calculate_theta) == true
    chi = chitransform(Streams, Drainage_Area, 'mn', theta_calc, 'correctcellsize', false);
else
    chi = chitransform(Streams, Drainage_Area, 'mn', Params.theta, 'a0', Params.minarea, 'correctcellsize', false);
end
chiMapped = mapfromnal(Flow_Direction, Streams, chi);

% Print space
disp(' ');

%% Export channel images (during specified plot interval)

% Execute at specified interval
if contains(List{t}, Print_List) == 1

    % Print status
    disp('Exporting channel images...')
    
    % Export Longitudinal Profile
    fig = figure();
    hold on;
    plotdz(Streams, Z_smooth);
    grid();
    title(['Year = ' erase(List{t}, '.tif')]);
    exportgraphics(fig, [out_directory '\Images\Longitudinal_Profile.' char(Params.file_type)], 'Resolution', Params.file_dpi);
    close fig 1;
    
    % Ksn_Map
    cm_gray = colormap('gray');
    close fig 1;
    fig = figure();
    hold on
    imageschs(DEM, [], 'colormap', cm_gray(1:255,:), 'colorbar', false);
    plotc(Streams, Ksn);
    colorbar;
    title(['Year = ' erase(List{t}, '.tif')]);
    exportgraphics(fig, [out_directory '\Images\Ksn_Map.' char(Params.file_type)], 'Resolution', Params.file_dpi);
    close fig 1;
    
    % Chi_Map
    fig = figure();
    hold on;
    imageschs(DEM, [], 'colormap', cm_gray(1:255,:), 'colorbar', false);
    plotc(Streams, chi, 'LineWidth', 3);
    colormap(cm);
    colorbar;
    title(['Year = ' erase(List{t}, '.tif')]);
    exportgraphics(fig, [out_directory '\Images\Chi_Map.' char(Params.file_type)], 'Resolution', Params.file_dpi);
    close fig 1;
    
    % Ksn_Mapped
    fig = figure(); 
    hold on;
    cm = flipud(ttscm('roma'));  % load roma color map and reverse it
    imageschs(DEM, KsnMapped, 'colormap', 'jet', 'colorbar', true, 'colorbarylabel', 'Ksn');
    plot(Divides_Lines, 'limit', [2000 inf], 'color', 'k');
    title(['Year = ' erase(List{t}, '.tif')]);
    exportgraphics(fig, [out_directory '\Images\Ksn_Mapped.' char(Params.file_type)], 'Resolution', Params.file_dpi);
    close fig 1;
    
    % Chi_Mapped
    fig = figure(); 
    hold on;
    imageschs(DEM, chiMapped, 'colormap', cm, 'colorbar', true);
    plot(Divides_Lines, 'limit', [2000 inf], 'color', 'k');
    title(['Year = ' erase(List{t}, '.tif')]);
    exportgraphics(fig, [out_directory '\Images\Chi_Mapped.' char(Params.file_type)], 'Resolution', Params.file_dpi);
    close fig 1;
    
    % Chi_Plot
    if char(Params.calculate_theta) == true
        cp = chiplot(Streams, Z_smooth, Drainage_Area, 'mn', theta_calc, 'color', 'k', 'plot', false);
    else
        cp = chiplot(Streams, Z_smooth, Drainage_Area, 'mn', Params.theta, 'color', 'k', 'plot', false);
    end
    close fig 1;
    fig = figure();
    hold on;
    plot(cp.chi, cp.elev);
    grid();
    title(['Year = ' erase(List{t}, '.tif')]);
    xlabel('\chi');
    ylabel('Elevation (m)');
    exportgraphics(fig, [out_directory '\Images\Chi_Plot.' char(Params.file_type)], 'Resolution', Params.file_dpi);
    close fig 1
    
    % Print space
    disp(' ');

end

%% Export channel grids (during specified plot interval)

% Execute at specified interval
if contains(List{t}, Print_List) == 1

    % Print status
    disp('Exporting channel grids...')
    
    % Export chi and ksn grids
    GRIDobj2geotiff(KsnMapped, [out_directory '\Grids\Ksn_Mapped_Grid.tif']);
    GRIDobj2geotiff(chiMapped, [out_directory '\Grids\Chi_Mapped_Grid.tif']);
    
    % Print space
    disp(' ');

end

%% Segment Picker (manual step for selecting channel heads)

% Sc = SegmentPicker_CS(DEM_original, Flow_Direction, Drainage_Area, Streams, 1, 'direction', 'down', 'plot_style', 'keep', 'ref_concavity', Params.theta, 'method', 'new_picks');
% PS = load('PickedSegments_1.mat');

%% Find / update channel heads

% Print status
disp('Handling channel heads...')

% During first iteration, find all channel heads within 1 km from eastern
% boundary
if t == 1

    % Find heads and add numeric ID
    Heads = streampoi(Streams);
    Heads = Heads(Heads(:, 1) >= 19050, :);
    [Heads(:, 2), id] = sort(Heads(:, 2), 'descend');
    Heads(:, 1) = Heads(id, 1);
    Heads(:, 3) = 1 : size(Heads, 1);

    % Create big distance array
    DISTANCE(1 : size(Heads, 1), 1) = 0;

% During each subsequent iteration, track changes in head positions    
elseif t > 1

    % Find new channel heads
    New_Heads = streampoi(Streams);

    % Initiate array to track distances between old and new heads
    distance = New_Heads(1 : end, 1) * NaN;

    DISTANCE = [DISTANCE transpose(((1 : size(DISTANCE, 1))) * NaN)];

    % Loop through old heads
    for i = 1 : size(Heads, 1)

        % Loop through new heads
        for j = 1 : size(New_Heads, 1)
        
            % For current old head, calculate distance to each new head
            distance(j, 1) = ( ( ( abs( Heads(i, 1) - New_Heads(j, 1) ) ) ^2 ) + ( ( abs( Heads(i, 2) - New_Heads(j, 2) ) ) ^ 0.5 ) ) ^ 0.5;

        end

        % Find index of closet new head
        r = find(distance == min(distance));

        % Update head
        Heads(i, 1 : 2) = New_Heads(r, :);

        % Record change in head position
        DISTANCE(i, t) = min(distance);

    end

end

% Print space
disp(' ');

%% Extract reversed channel

% Print status
disp('Extracting reversed channel...')

% Find channel head closest to northwestern domain corner
[x, y] = getcoordinates(DEM);
x = x(1);
y = y(1);
[x, y] = snap2stream(Streams, x, y);
rev_head = [x y 1];

% Create reversed channel stream object
Sc_rev = SegmentPicker_CS(DEM_original, Flow_Direction, Drainage_Area, Streams, 1, 'direction', 'down', 'plot_style', 'keep', 'ref_concavity', Params.theta, 'method', 'prev_picks', 'picks', rev_head);
close fig 1

% Print space
disp(' ');

%% Loop through channel heads and find nearest "knickpoint" downstream

% Print status
disp('Looping through channel heads along the eastern domain boundary and finding nearest "knickpoints" downstream...')

% Initialize knickpoints coordinate arrays
Xk = NaN;
Yk = NaN;
IXk = NaN;
dzdtk = NaN;
chik = NaN;
Xc = NaN;
Yc = NaN;
chic = NaN;

% Loop through heads
for i = 1 : size(Heads, 1)

    % Create segment
    Sc = SegmentPicker_CS(DEM_original, Flow_Direction, Drainage_Area, Streams, 1, 'direction', 'down', 'plot_style', 'keep', 'ref_concavity', Params.theta, 'method', 'prev_picks', 'picks', Heads(i, :));
    close fig 1

    % Sort Sc.IXgrid based on distance from outlet
    [distance_sorted, sortIdx] = sort(Sc.distance, 'descend');
    IXgrid_sorted = Sc.IXgrid(sortIdx);
    Y_sorted = Sc.y(sortIdx);
    X_sorted = Sc.x(sortIdx);
    
    % Initiate tickers
    j = 1;
    q = 1;

    % Loop through channel and find "knickpoint"
    while q == 1

        % If dz / dt at node is > 0, skip node
        if dzdt_grid.Z(IXgrid_sorted(j)) >= 0
            j = j + 1;

            % If node was the outlet, end loop
            if j > size(IXgrid_sorted, 1)
                q = 0;
            end

        % If dz / dt at node is negative but greater than u / 100, skip
        % node. This is implemented to account for some minor transience
        % that occurs near the northernmost channel head and is unrelated
        % to the base level signal.
        elseif dzdt_grid.Z(IXgrid_sorted(j)) >= -u / 100
            j = j + 1;

            % If node was the outlet, end loop
            if j > size(IXgrid_sorted, 1)
                q = 0;
            end

        % If dz / dt is < 0, identify knickpoint
        else
            q = 0;
        end

    end
    
    % Record knickpoint coordinates
    if j <= size(IXgrid_sorted, 1)
        Yk(i, 1) = Y_sorted(j);
        Xk(i, 1) = X_sorted(j);
        IXk(i, 1) = IXgrid_sorted(j);
        dzdtk(i, 1) = dzdt_grid.Z(IXgrid_sorted(j));
        chik(i, 1) = chiMapped.Z(IXgrid_sorted(j));

    % If no knickpoint was detected, record NaN
    elseif j > size(IXgrid_sorted, 1)
        Yk(i, 1) = NaN;
        Xk(i, 1) = NaN;
        IXk(i, 1) = NaN;
        dzdtk(i, 1) = NaN;
        chik(i, 1) = NaN;
    end

    % If recorded knickpoint is the channel head, consider the knickpoint
    % to have left the channel network and record NaN
    if j == 1
        Yk(i, 1) = NaN;
        Xk(i, 1) = NaN;
        IXk(i, 1) = NaN;
        dzdtk(i, 1) = NaN;
        chik(i, 1) = NaN;
    end

    % If some amount of filling occured along the segment during
    % preprocessing, consider the knickpoint to be the downstream-most
    % extent of filling (scheme assumes that preprocessing method was 
    % "fill" and not "carve")
    if sum(preprocessing_change(IXgrid_sorted)) > 0
        j2 = size(preprocessing_change(IXgrid_sorted), 1);
        q2 = 1;
        while q2 == 1
            if preprocessing_change(IXgrid_sorted(j2)) == 0
                j2 = j2 - 1;
            else
                q2 = 0;
            end
        end
        Yk(i, 1) = Y_sorted(j2);
        Xk(i, 1) = X_sorted(j2);
        IXk(i, 1) = IXgrid_sorted(j2);
        dzdtk(i, 1) = dzdt_grid.Z(IXgrid_sorted(j2));
        chik(i, 1) = chiMapped.Z(IXgrid_sorted(j2));
    end

    % Get chi value where segment crosses the contact. Scheme assumes that
    % contact is located at x = 15000 m.
    contact_node = find(X_sorted == 14950);   
    contact_node = max(contact_node);
    Xc(i, 1) = X_sorted(contact_node);
    Yc(i, 1) = Y_sorted(contact_node);
    chic(i, 1) = chiMapped.Z(IXgrid_sorted(contact_node));

    % Reinitiate tickers
    j = 1;
    q = 1;

    % Loop through channel and find junction with reversed channel
    while q == 1

        % If node is not junction, move on to next node downstream
        if ismember(IXgrid_sorted(j), Sc_rev.IXgrid) == 0
            j = j + 1;

        % If node is junction, end loop
        elseif ismember(IXgrid_sorted(j), Sc_rev.IXgrid) == 1
            q = 0;

        end
    
    end

    % Record junction information
    Xo(i, 1) = X_sorted(j);
    Yo(i, 1) = Y_sorted(j);
    chio(i, 1) = chiMapped.Z(IXgrid_sorted(j));
    IXo(i, 1) = IXgrid_sorted(j);

end

% Track chi changes through time
if t == 1
    CHIK = chik;
    CHIC = chic;
    CHIO = chio;
elseif t > 1
    CHIK = [CHIK, chik];
    CHIC = [CHIC, chic];
    CHIO = [CHIO, chio];
end

% Print space
disp(' ');

%% Analysis Plots

% Print status
disp('Exporting knickpoint analysis plots...')

% DEM_Knick_Map
fig = figure();
ax1 = axes;
hold on;
imageschs(DEM, [], 'colormap', cc, 'colorbarylabel', 'Elevation (m)');
plot(Streams, 'LineWidth', 1, 'Color', 'b');
ax2 = axes;
hold on;
scatter(Xk, Yk, 50, chik, 'filled', 'MarkerEdgeColor', 'k');
scatter(Xo, Yo, 50, chio, 'filled', "v", 'MarkerEdgeColor', 'k');
scatter(Xc, Yc, 50, chic, 'filled', "s", 'MarkerEdgeColor', 'k');
ax1.XTick = [];
ax1.YTick = [];
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
colormap(ax1, cc);
colormap(ax2, 'parula');
ax2.XLim =ax1.XLim;
ax2.YLim =ax1.YLim;
ax2.OuterPosition =ax1.OuterPosition;
ax2.InnerPosition =ax1.InnerPosition;
ax2.Position =ax1.Position;
ax2.PositionConstraint =ax1.PositionConstraint;
ax2.DataAspectRatio =ax1.DataAspectRatio;
ax2.DataAspectRatioMode =ax1.DataAspectRatioMode;
set([ax1, ax2], 'Position', [.17 .11 .685 .815]);        % [.17 .11 .685 .815]
cb1 = colorbar(ax1, 'Position', [.16 .11 .03 .815]);   % [.05 .11 .0675 .815]
cb1.Label.String = 'Elevation (m)';
cb2 = colorbar(ax2, 'Position', [.84 .11 .03 .815]);  % [.88 .11 .0675 .815]
cb2.Label.String = '\chi';
title(ax1, ['Year = ' erase(List{t}, '.tif')]);
exportgraphics(fig, [out_directory '\Images\DEM_Knick_Chi_Map.' char(Params.file_type)], 'Resolution', Params.file_dpi);
close fig 1

% ChiMapped_Knick_Chi_Map
fig = figure();
ax1 = axes;
hold on;
imageschs(DEM, chiMapped, 'colormap', cm, 'colorbar', true, 'colorbarylabel', '\chi');
plot(Streams, 'LineWidth', 1, 'Color', 'b');
ax2 = axes;
hold on;
scatter(Xk, Yk, 50, chik, 'filled', 'MarkerEdgeColor', 'k');
scatter(Xo, Yo, 50, chio, 'filled', "v", 'MarkerEdgeColor', 'k');
scatter(Xc, Yc, 50, chic, 'filled', "s", 'MarkerEdgeColor', 'k');
ax1.XTick = [];
ax1.YTick = [];
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
colormap(ax1, cm);
colormap(ax2, 'parula');
ax2.XLim =ax1.XLim;
ax2.YLim =ax1.YLim;
ax2.OuterPosition =ax1.OuterPosition;
ax2.InnerPosition =ax1.InnerPosition;
ax2.Position =ax1.Position;
ax2.PositionConstraint =ax1.PositionConstraint;
ax2.DataAspectRatio =ax1.DataAspectRatio;
ax2.DataAspectRatioMode =ax1.DataAspectRatioMode;
set([ax1, ax2], 'Position', [.17 .11 .685 .815]);        % [.17 .11 .685 .815]
cb1 = colorbar(ax1, 'Position', [.16 .11 .03 .815]);   % [.05 .11 .0675 .815]
cb1.Label.String = '\chi';
cb2 = colorbar(ax2, 'Position', [.84 .11 .03 .815]);  % [.88 .11 .0675 .815]
cb2.Label.String = '\chi at knickpoint';
title(ax1, ['Year = ' erase(List{t}, '.tif')]);
exportgraphics(fig, [out_directory '\Images\ChiMapped_Knick_Chi_Map.' char(Params.file_type)], 'Resolution', Params.file_dpi);
close fig 1

% chio_v_chik
fig = figure();
hold on;
mx = round(max(chio));
mn = round(min(chio));
line = mn : mx;
scatter(chio, chik, 50, Yk, 'filled', 'MarkerEdgeColor', 'k');
plot(line, line, 'k--');
colormap('plasma')
cb = colorbar();
ylabel(cb, 'Northing (m)');
grid();
xlabel('\chi at junction');
ylabel('Knickpoint \chi');
title(['Year = ' erase(List{t}, '.tif')]);
exportgraphics(fig, [out_directory '\Images\chio_v_chik.' char(Params.file_type)], 'Resolution', Params.file_dpi);
close fig 1;

% chic_v_chik
fig = figure();
hold on;
mx = round(max(chic));
mn = round(min(chic));
line = mn : mx;
scatter(chic, chik, 50, Yk, 'filled', 'MarkerEdgeColor', 'k');
plot(line, line, 'k--');
colormap('plasma')
cb = colorbar();
ylabel(cb, 'Northing (m)');
grid();
xlabel('\chi at lithologic contact downstream from knickpoint');
ylabel('Knickpoint \chi');
exportgraphics(fig, [out_directory '\Images\chic_v_chik.' char(Params.file_type)], 'Resolution', Params.file_dpi);
close fig 1;

% Print space
disp(' ');

%% End DEM loop

end

%% Time v chi plots

% Print status
disp('Exporting timeseries plots and data...')

% chik_timeseries
fig = figure();
hold on;
plot(List_nums, CHIK(1, :), 'r');
plot(List_nums, CHIK(end - 1, :), 'b');
grid();
xlabel('Model year');
ylabel('Knickpoint \chi');
legend('Knickpoint (north)', 'Knickpoint (south)', 'location', 'northwest');
exportgraphics(fig, [analysis_directory '\Output\Timeseries\chik_timeseries.' char(Params.file_type)], 'Resolution', Params.file_dpi);
close fig 1;

% dchik_dt_timeseries
dchik_dt_N = List_nums * NaN;
dchik_dt_S = List_nums * NaN;
for i = 1 : size(CHIK, 2)
    if i == 1
        dchik_dt_N(i) = CHIK(1, 1) / (List_nums(1));
        dchik_dt_S(i) = CHIK(end - 1, 1) / (List_nums(1));
    else
        dchik_dt_N(i) = (CHIK(1, i) - CHIK(1, i - 1)) / (List_nums(i) - List_nums(i - 1));
        dchik_dt_S(i) = (CHIK(end - 1, i) - CHIK(end - 1, i - 1)) / (List_nums(i) - List_nums(i - 1));
    end
end
ksp_1 = spin_params.ksp_1;
ksp_2 = spin_params.ksp_2;
rate_1 = ksp_1 * (Params.minarea ^ Params.theta);
rate_2 = ksp_2 * (Params.minarea ^ Params.theta);
fig = figure();
hold on;
plot(List_nums, dchik_dt_N, 'r');
plot(List_nums, dchik_dt_S, 'b');
yline(rate_1, 'k--');
yline(rate_2, 'k:');
grid();
xlabel('Model year');
ylabel('d\chi / dt');
legend('Knickpoint (north)', 'Knickpoint (south)', ['Predicted rate (Ksp = ' num2str(ksp_1) ')'], ['Predicted rate (Ksp = ' num2str(ksp_2) ')'])
exportgraphics(fig, [analysis_directory '\Output\Timeseries\dchik_dt_timeseries.' char(Params.file_type)], 'Resolution', Params.file_dpi);
close fig 1;

% chikN_div_chikS_timeseries
fig = figure();
hold on;
plot(List_nums, CHIK(1, :) ./ CHIK(end - 1, :), 'k');
grid();
xlabel('Model year');
ylabel('\chi northern knickpoint / \chi southern knickpoint');
exportgraphics(fig, [analysis_directory '\Output\Timeseries\chikN_div_chikS_timeseries.' char(Params.file_type)], 'Resolution', Params.file_dpi);
close fig 1;

% dchik_dt_N_div_dchik_dt_S_timeseries
fig = figure();
hold on;
plot(List_nums, dchik_dt_N ./ dchik_dt_S, 'k');
grid();
xlabel('Model year');
ylabel('Northern knickpoint (d\chi / dt) / southern knickpoint (d\chi / dt)');
exportgraphics(fig, [analysis_directory '\Output\Timeseries\dchik_dt_N_div_dchik_dt_S_timeseries.' char(Params.file_type)], 'Resolution', Params.file_dpi);
close fig 1;

% Export timeseries
csvwrite([analysis_directory '\Output\Timeseries\List_nums.csv'], List_nums);
csvwrite([analysis_directory '\Output\Timeseries\CHIK.csv'], CHIK);
csvwrite([analysis_directory '\Output\Timeseries\CHIC.csv'], CHIC);
csvwrite([analysis_directory '\Output\Timeseries\CHIO.csv'], CHIO);
csvwrite([analysis_directory '\Output\Timeseries\DISTANCE.csv'], DISTANCE);

% Print space
disp(' ');

%% Videos (images exported during plot interval)

% Print status
disp('Creating videos...')

% Get list of images
contents = dir([analysis_directory '\Output\' erase(Print_List{1}, '.tif') '\Images']);
names = {contents.name};
im_List = names(3 : end);

% Loop through videos
for j = 1 : size(im_List, 2)

    % Create images cells
    images = cell(size(Print_List, 2), 1);

    % Loop through images
    for i = 1 : size(Print_List, 2)
        
        % Open image
        images{i} = imread([analysis_directory '\Output\' erase(Print_List{i}, '.tif') '\Images\' im_List{j}]);

        % Enforce image size
        if i == 1
            px_sz = size(images{1});
            px_sz = px_sz(1, 1 : 2);
        else
            images{i} = imresize(images{i}, px_sz);
        end

    end
    
    % Initialize VideoWriter
    writerObj = VideoWriter([analysis_directory '\Output\Videos\' erase(im_List{j}, ['.' char(Params.file_type)]) '.avi']);
    writerObj.FrameRate = 1;
    secsPerImage = ones(1, size(Print_List, 2));
    open(writerObj);
    
    % Write frames
    for u=1:length(images)
        frame = im2frame(images{u});
        for v=1:secsPerImage(u) 
            writeVideo(writerObj, frame);
        end
    end
    
    % Close VideoWriter
    close(writerObj);

end

% Print space
disp(' ');

%% Videos (images exported every timestep)

% Print status
disp('Creating more videos...')

% Get list of images
contents = dir([analysis_directory '\Output\' erase(List{1}, '.tif') '\Images']);
names = {contents.name};
im_List = names(3 : end);

% Loop through videos
for j = 1 : size(im_List, 2)

    % Create images cells
    images = cell(size(List, 2), 1);

    % Loop through images
    for i = 1 : size(List, 2)
        
        % Open image
        images{i} = imread([analysis_directory '\Output\' erase(List{i}, '.tif') '\Images\' im_List{j}]);

        % Enforce image size
        if i == 1
            px_sz = size(images{1});
            px_sz = px_sz(1, 1 : 2);
        else
            images{i} = imresize(images{i}, px_sz);
        end

    end
    
    % Initialize VideoWriter
    writerObj = VideoWriter([analysis_directory '\Output\Videos\' erase(im_List{j}, ['.' char(Params.file_type)]) '.avi']);
    writerObj.FrameRate = 1;
    secsPerImage = ones(1, size(List, 2));
    open(writerObj);
    
    % Write frames
    for u=1:length(images)
        frame = im2frame(images{u});
        for v=1:secsPerImage(u) 
            writeVideo(writerObj, frame);
        end
    end
    
    % Close VideoWriter
    close(writerObj);

end

% Print space
disp(' ');

 %% Complete

% Print status
disp('Program complete. Happy trails!')
function fig = SixDOFanimation(varargin)

    %% Create local variables

    % Required arguments
    p = varargin{1};                % position of body
    Cor = varargin{2};                % rotation matrix of body
    Cell = varargin{3};
    [numSamples dummy] = size(p);

    % Default values of optional arguments
    SamplePlotFreq = 1;
    Trail = 'Off';
    LimitRatio = 1;
    Position = [];
    FullScreen = false;
    View = [30 20];
    AxisLength = 1;
    ShowArrowHead = 'on';
    Xlabel = 'X';
    Ylabel = 'Y';
    Zlabel = 'Z';
    Title = 'UWB Positioning';
    ShowLegend = true;
    CreateAVI = false;
    AVIfileName = 'UWB Positioning';
    AVIfileNameEnum = true;
    AVIfps = 30;

    for i = 4:2:nargin
        if  strcmp(varargin{i}, 'SamplePlotFreq'), SamplePlotFreq = varargin{i+1};
        elseif  strcmp(varargin{i}, 'Trail')
            Trail = varargin{i+1};
            if(~strcmp(Trail, 'Off') && ~strcmp(Trail, 'DotsOnly') && ~strcmp(Trail, 'All'))
                error('Invalid argument.  Trail must be ''Off'', ''DotsOnly'' or ''All''.');
            end
        elseif  strcmp(varargin{i}, 'LimitRatio'), LimitRatio = varargin{i+1};
        elseif  strcmp(varargin{i}, 'Position'), Position = varargin{i+1};
        elseif  strcmp(varargin{i}, 'FullScreen'), FullScreen = varargin{i+1};
        elseif  strcmp(varargin{i}, 'View'), View = varargin{i+1};
        elseif  strcmp(varargin{i}, 'AxisLength'), AxisLength = varargin{i+1};
        elseif  strcmp(varargin{i}, 'ShowArrowHead'), ShowArrowHead = varargin{i+1};
        elseif  strcmp(varargin{i}, 'Xlabel'), Xlabel = varargin{i+1};
        elseif  strcmp(varargin{i}, 'Ylabel'), Ylabel = varargin{i+1};
        elseif  strcmp(varargin{i}, 'Zlabel'), Zlabel = varargin{i+1};
        elseif  strcmp(varargin{i}, 'Title'), Title = varargin{i+1};
        elseif  strcmp(varargin{i}, 'ShowLegend'), ShowLegend = varargin{i+1};
        elseif  strcmp(varargin{i}, 'CreateAVI'), CreateAVI = varargin{i+1};
        elseif  strcmp(varargin{i}, 'AVIfileName'), AVIfileName = varargin{i+1};
        elseif  strcmp(varargin{i}, 'AVIfileNameEnum'), AVIfileNameEnum = varargin{i+1};
        elseif  strcmp(varargin{i}, 'AVIfps'), AVIfps = varargin{i+1};
        else error('Invalid argument.');
        end
    end;

    %% Reduce data to samples to plot only

    p = p(1:SamplePlotFreq:numSamples, :);
%     R = R(:, :, 1:SamplePlotFreq:numSamples) * AxisLength;
%     if(numel(View) > 2)
%         View = View(1:SamplePlotFreq:numSamples, :);
%     end
    [numPlotSamples dummy] = size(p);

    %% Setup AVI file

%     aviobj = [];                                                            	% create null object
%     if(CreateAVI)
%         fileName = strcat(AVIfileName, '.avi');
%         if(exist(fileName, 'file'))
%             if(AVIfileNameEnum)                                              	% if file name exists and enum enabled
%                 i = 0;
%                 while(exist(fileName, 'file'))                                  % find un-used file name by appending enum
%                     fileName = strcat(AVIfileName, sprintf('%i', i), '.avi');
%                     i = i + 1;
%                 end
%             else                                                                % else file name exists and enum disabled
%                 fileName = [];                                                  % file will not be created
%             end
%         end
%         if(isempty(fileName))
%             sprintf('AVI file not created as file already exists.')
%         else
%             aviobj = avifile(fileName, 'fps', AVIfps, 'compression', 'Cinepak', 'quality', 100);
%         end
%     end

    %% Setup figure and plot

    % Create figure
    fig = figure('NumberTitle', 'off', 'Name', 'UWB Positioning');
    hold on
    plot(Cor(:,2),Cor(:,3),'r*');
    hold on
    text_01_handle = text(Cor(1,2),Cor(1,3),' An1');
    text_02_handle = text(Cor(2,2),Cor(2,3),' An2');
    text_03_handle = text(Cor(3,2),Cor(3,3),' An3');
    text_04_handle = text(Cor(4,2),Cor(4,3),' An4');
    text_05_handle = text(Cor(5,2),Cor(5,3),' An5');
    text_06_handle = text(Cor(6,2),Cor(6,3),' An6');
    text_07_handle = text(Cor(7,2),Cor(7,3),' An7');
    text_08_handle = text(Cor(8,2),Cor(8,3),' An8');
    text_09_handle = text(Cor(9,2),Cor(9,3),' An9');
    text_10_handle = text(Cor(10,2),Cor(10,3),' An10');
    hold on
    
    if(FullScreen)
        screenSize = get(0, 'ScreenSize');
        set(fig, 'Position', [0 0 screenSize(3) screenSize(4)]);
    elseif(~isempty(Position))
        set(fig, 'Position', Position);
    end
    set(gca, 'drawmode', 'fast');
    lighting phong;
    set(gcf, 'Renderer', 'zbuffer');
    hold on;
    axis equal;
    grid on;
%     view(View(1, 1), View(1, 2));
    title(i);
    xlabel(Xlabel);
    ylabel(Ylabel);
    zlabel(Zlabel);

    % Create plot data arrays
    if(strcmp(Trail, 'DotsOnly') || strcmp(Trail, 'All'))
        x = zeros(numPlotSamples, 1);
        y = zeros(numPlotSamples, 1);
        z = zeros(numPlotSamples, 1);
    end
%     if(strcmp(Trail, 'All'))
%         ox = zeros(numPlotSamples, 1);
%         oy = zeros(numPlotSamples, 1);
%         oz = zeros(numPlotSamples, 1);
%         ux = zeros(numPlotSamples, 1);
%         vx = zeros(numPlotSamples, 1);
%         wx = zeros(numPlotSamples, 1);
%         uy = zeros(numPlotSamples, 1);
%         vy = zeros(numPlotSamples, 1);
%         wy = zeros(numPlotSamples, 1);
%         uz = zeros(numPlotSamples, 1);
%         vz = zeros(numPlotSamples, 1);
%         wz = zeros(numPlotSamples, 1);
%     end
    x(1) = p(100,1);
    y(1) = p(100,2);
    z(1) = p(100,3);
%     ox(1) = x(1);
%     oy(1) = y(1);
%     oz(1) = z(1);
%     ux(1) = R(1,1,1:1);
%     vx(1) = R(2,1,1:1);
%     wx(1) = R(3,1,1:1);
%     uy(1) = R(1,2,1:1);
%     vy(1) = R(2,2,1:1);
%     wy(1) = R(3,2,1:1);
%     uz(1) = R(1,3,1:1);
%     vz(1) = R(2,3,1:1);
%     wz(1) = R(3,3,1:1);

    % Create graphics handles
    orgHandle = plot(x, y, 'k.');
%     if(ShowArrowHead)
%         ShowArrowHeadStr = 'on';
%     else
%         ShowArrowHeadStr = 'off';
%     end
%     quivXhandle = quiver3(ox, oy, oz, ux, vx, wx,  'r', 'ShowArrowHead', ShowArrowHeadStr, 'MaxHeadSize', 0.999999, 'AutoScale', 'off');
%     quivYhandle = quiver3(ox, oy, oz, uy, vy, wy,  'g', 'ShowArrowHead', ShowArrowHeadStr, 'MaxHeadSize', 0.999999, 'AutoScale', 'off');
%     quivZhandle = quiver3(ox, ox, oz, uz, vz, wz,  'b', 'ShowArrowHead', ShowArrowHeadStr, 'MaxHeadSize', 0.999999, 'AutoScale', 'off');

    % Create legend
%     if(ShowLegend)
%         legend('Origin', 'X', 'Y');
%     end
    
    % Set initial limits
%     Xlim = [x(1)-AxisLength x(1)+AxisLength] * LimitRatio;
%     Ylim = [y(1)-AxisLength y(1)+AxisLength] * LimitRatio;
%     Zlim = [z(1)-AxisLength z(1)+AxisLength] * LimitRatio;
    set(gca, 'Xlim', [-5 25], 'Ylim', [-8 12]);
    
    % Set initial view
%     view(View(1, :));

    %% Plot one sample at a time

    for i = 1:numPlotSamples

        % Update graph title
        if(strcmp(Title, ''))
            titleText = sprintf('Sample %i of %i', 1+((i-1)*SamplePlotFreq), numSamples);
        else
            titleText = strcat(Title, ' (', sprintf('Sample %i of %i', 1+((i-1)*SamplePlotFreq), numSamples), ')');
        end
        title(titleText);

        % Plot body x y z axes
        if(strcmp(Trail, 'DotsOnly') || strcmp(Trail, 'All'))
            x(1:i) = p(1:i,1);
            y(1:i) = p(1:i,2);
            z(1:i) = p(1:i,3);
        else
            x = p(i,1);
            y = p(i,2);
            z = p(i,3);
        end
%         if(strcmp(Trail, 'All'))
%             ox(1:i) = p(1:i,1);
%             oy(1:i) = p(1:i,2);
%             oz(1:i) = p(1:i,3);
%             ux(1:i) = R(1,1,1:i);
%             vx(1:i) = R(2,1,1:i);
%             wx(1:i) = R(3,1,1:i);
%             uy(1:i) = R(1,2,1:i);
%             vy(1:i) = R(2,2,1:i);
%             wy(1:i) = R(3,2,1:i);
%             uz(1:i) = R(1,3,1:i);
%             vz(1:i) = R(2,3,1:i);
%             wz(1:i) = R(3,3,1:i);
%         else
%             ox = p(i,1);
%             oy = p(i,2);
%             oz = p(i,3);
%             ux = R(1,1,i);
%             vx = R(2,1,i);
%             wx = R(3,1,i);
%             uy = R(1,2,i);
%             vy = R(2,2,i);
%             wy = R(3,2,i);
%             uz = R(1,3,i);
%             vz = R(2,3,i);
%             wz = R(3,3,i);
%         end
        set(orgHandle, 'xdata', x, 'ydata', y);
        
        if Cell(1,i)==1
            set(text_01_handle,'Color','green');
            set(text_02_handle,'Color','green');
            set(text_03_handle,'Color','green');
            set(text_04_handle,'Color','green');
            set(text_05_handle,'Color','black');
            set(text_06_handle,'Color','black');
            set(text_07_handle,'Color','black');
            set(text_08_handle,'Color','black');   
            set(text_09_handle,'Color','black');
            set(text_10_handle,'Color','black'); 
        else if Cell(1,i)==3
                set(text_01_handle,'Color','black');
                set(text_02_handle,'Color','black');
                set(text_03_handle,'Color','green');
                set(text_04_handle,'Color','green');
                set(text_05_handle,'Color','green');
                set(text_06_handle,'Color','green');
                set(text_07_handle,'Color','black');
                set(text_08_handle,'Color','black');   
                set(text_09_handle,'Color','black');
                set(text_10_handle,'Color','black');                 
            else if Cell(1,i)==7
                    set(text_01_handle,'Color','black');
                    set(text_02_handle,'Color','black');
                    set(text_03_handle,'Color','black');
                    set(text_04_handle,'Color','black');
                    set(text_05_handle,'Color','black');
                    set(text_06_handle,'Color','black');
                    set(text_07_handle,'Color','green');
                    set(text_08_handle,'Color','green');   
                    set(text_09_handle,'Color','green');
                    set(text_10_handle,'Color','green');                       
                end
            end
        end
%         set(quivXhandle, 'xdata', ox, 'ydata', oy, 'zdata', oz,'udata', ux, 'vdata', vx, 'wdata', wx);
%         set(quivYhandle, 'xdata', ox, 'ydata', oy, 'zdata', oz,'udata', uy, 'vdata', vy, 'wdata', wy);
%         set(quivZhandle, 'xdata', ox, 'ydata', oy, 'zdata', oz,'udata', uz, 'vdata', vz, 'wdata', wz);

        % Adjust axes for snug fit and draw
%         axisLimChanged = false;
%         if((p(i,1) - AxisLength) < Xlim(1)), Xlim(1) = p(i,1) - LimitRatio*AxisLength; axisLimChanged = true; end
%         if((p(i,2) - AxisLength) < Ylim(1)), Ylim(1) = p(i,2) - LimitRatio*AxisLength; axisLimChanged = true; end
%         if((p(i,3) - AxisLength) < Zlim(1)), Zlim(1) = p(i,3) - LimitRatio*AxisLength; axisLimChanged = true; end
%         if((p(i,1) + AxisLength) > Xlim(2)), Xlim(2) = p(i,1) + LimitRatio*AxisLength; axisLimChanged = true; end
%         if((p(i,2) + AxisLength) > Ylim(2)), Ylim(2) = p(i,2) + LimitRatio*AxisLength; axisLimChanged = true; end
%         if((p(i,3) + AxisLength) > Zlim(2)), Zlim(2) = p(i,3) + LimitRatio*AxisLength; axisLimChanged = true; end
%         if(axisLimChanged), set(gca, 'Xlim', Xlim, 'Ylim', Ylim); end
        drawnow;

        % Adjust view
%         if(numel(View) > 2)
%             view(View(i, :));
%         end

        % Add frame to AVI object
%         if(~isempty(aviobj))
%             frame = getframe(fig);
%             aviobj = addframe(aviobj, frame);
%         end

    end

    hold off;

    % Close AVI file
%     if(~isempty(aviobj))
%         aviobj = close(aviobj);
%     end

end
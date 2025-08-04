clear; clc; format short g; format compact;
% This was added by Jake
zap

%% Add Function Path and Lookup .TIF Files
Lookup = FileLookup('tif');

%% Read TIF Files
for i = 1:Lookup.FileCount
    Image.All{i} = double(imread(fullfile(Lookup.FolderAddress, Lookup.FolderInfo(i).name)));
end

%% Separate Green/Red Channels
for i = 1:Lookup.FileCount
    if contains(Lookup.FolderInfo(i).name, 'green')
        Image.Green.Data{i} = Image.All{i};
        Image.Green.Name{i} = Lookup.FolderInfo(i).name;
    elseif contains(Lookup.FolderInfo(i).name, 'red')
        Image.Red.Data{i} = Image.All{i};
        Image.Red.Name{i} = Lookup.FolderInfo(i).name;
    end
end

Image.Green.Data = Image.Green.Data(~cellfun(@isempty, Image.Green.Data));
Image.Green.Name = Image.Green.Name(~cellfun(@isempty, Image.Green.Name));
Image.Red.Data = Image.Red.Data(~cellfun(@isempty, Image.Red.Data));
Image.Red.Name = Image.Red.Name(~cellfun(@isempty, Image.Red.Name));

ImageSet(Image.Green, "Green Channel - Unaligned", "Data");

%% Registeration
% Register Green Channel
[Image.Green, GreenTransforms, RegistrationPoints] = ManualRegisterByPoint(Image.Green);

% Apply Registration to Red Channel
Image.Red = ApplyTransforms(Image.Red, GreenTransforms, size(Image.Green.Registered{1}));

% Display Registered Images
ImageSet(Image.Green, "Green Aligned", "Registered")
ImageSet(Image.Red, "Red Aligned", "Registered")

% Draw ROI on Registered Green Image
figure;
[ROI.Green, ~] = DefineSingleROIOnRegistered(Image.Green, 'Registered');

%save the ROI
GreenROIInfo.MaskSize = {size(ROI.Green.Mask)};
GreenROIInfo.NumPixels = sum(ROI.Green.Mask(:));
GreenROIInfo.Coordinates = {ROI.Green.UserDefined.Position};
GreenROIInfo.Channel = "Green";
GreenROIInfo.ROI_Name = "Green_ROI";
GreenROIInfo.Type = "Freehand";
GreenROIInfoTable = struct2table(GreenROIInfo, 'AsArray', true);

% Draw Multiple Named ROIs on Registered Red Image
figure;
ROI.Red = DefineMultipleROIsOnRegistered(Image.Red, 'Registered');

%save the ROIs
for r = 1:numel(ROI.Red)
    RedROIInfo(r).MaskSize = {size(ROI.Red(r).Mask)};
    RedROIInfo(r).NumPixels = sum(ROI.Red(r).Mask(:));
    RedROIInfo(r).Coordinates = {ROI.Red(r).Coordinates};
    RedROIInfo(r).Channel = "Red";
    RedROIInfo(r).ROI_Name = ROI.Red(r).Name;
    RedROIInfo(r).Type = "Freehand";
end

RedROIInfoTable = struct2table(RedROIInfo, 'AsArray', true);

% Overlay Green ROI on Registered Green Images
OverlayROIOnImages(Image.Green, ROI.Green.Mask, [.1 .8 .8], 'Green Channel with ROI');

% Overlay Red ROIs on Registered Red Images
OverlayMultipleROIs(Image.Red, ROI.Red, 'Red Channel with ROIs');

%% Extract ROI Intensity from Green Images
for i = 1:numel(Image.Green.Registered)
    region = Image.Green.Registered{i}(ROI.Green.Mask);
    ROI.GreenInfo(i).Channel = "Green";
    ROI.GreenInfo(i).Name = Image.Green.Name{i};
    ROI.GreenInfo(i).ROI_Name = "Green_ROI";
    ROI.GreenInfo(i).Mean = mean(region);
    ROI.GreenInfo(i).StdDev = std(region);
    ROI.GreenInfo(i).Pixels = numel(region);
end

%% Extract ROI Intensities from Each Red ROI
idx = 1;
for r = 1:numel(ROI.Red)
    mask = ROI.Red(r).Mask;
    for i = 1:numel(Image.Red.Registered)
        region = Image.Red.Registered{i}(mask);
        ROI.RedInfo(idx).Channel = "Red";
        ROI.RedInfo(idx).Name = Image.Red.Name{i};
        ROI.RedInfo(idx).ROI_Name = ROI.Red(r).Name;
        ROI.RedInfo(idx).Mean = mean(region);
        ROI.RedInfo(idx).StdDev = std(region);
        ROI.RedInfo(idx).Pixels = numel(region);
        idx = idx + 1;
    end
end

%% Combine and Save Results
AllResults = [ROI.GreenInfo, ROI.RedInfo];
AllResultsTable = struct2table(AllResults);
disp(AllResultsTable);

filename = uniqueFilename(Lookup.FolderAddress, 'ROI_Intensity_Results', '.xlsx');
writetable(AllResultsTable, filename);
fprintf('Saved results to %s\n', filename);

%save ROIs
AllROITable = [GreenROIInfoTable; RedROIInfoTable];
roiFilename = uniqueFilename(Lookup.FolderAddress, 'ROI_Metadata', '.xlsx');
writetable(AllROITable, roiFilename);
fprintf('Saved ROI metadata to %s\n', roiFilename);

%save registration
RegistrationTable = struct2table(RegistrationPoints);
regFilename = uniqueFilename(Lookup.FolderAddress, 'Registration_Points', '.xlsx');
writetable(RegistrationTable, regFilename);
fprintf('Saved registration points to %s\n', regFilename);

%% ---------- FUNCTIONS ----------

function Lookup = FileLookup(extension)
    Lookup.FolderAddress = uigetdir(pwd, sprintf("Select folder with .%s files", extension));
    Lookup.FolderInfo = dir(fullfile(Lookup.FolderAddress, ['*.', extension]));
    Lookup.FileCount = numel(Lookup.FolderInfo);
end

function ImageSet(ImageStruct, Channel, ImageType)
    figure('Name', Channel)
    tiledlayout('flow');
    for i = 1:numel(ImageStruct.(ImageType))
        ax = nexttile;
        imagesc(ax, ImageStruct.(ImageType){i});
        clim([0 5000]);
        colormap gray; axis image off;
        title(ax, string(ImageStruct.Name{i}), 'Interpreter', 'none');
    end
end

function [ImageStruct, Transforms, RegistrationPoints] = ManualRegisterByPoint(ImageStruct)
    referenceImage = ImageStruct.Data{1};
    RegistrationPoints = struct(); 

    figure('Name', 'Select reference point');
    imagesc(referenceImage); colormap gray; axis image;
    clim([0 5000]);
    title('Click on reference point in the FIRST image');
    h = drawpoint('Color', 'cyan');
    pos = h.Position;
    xRef = pos(1); yRef = pos(2);
    close;

    RegistrationPoints(1).ImageName = ImageStruct.Name{1};
    RegistrationPoints(1).X = xRef;
    RegistrationPoints(1).Y = yRef;

    Transforms = cell(1, numel(ImageStruct.Data));
    ImageStruct.Registered = cell(1, numel(ImageStruct.Data));
    ImageStruct.Registered{1} = referenceImage;
    Transforms{1} = [0, 0];

    for i = 2:numel(ImageStruct.Data)
        currImage = ImageStruct.Data{i};
        figure('Name', sprintf('Select point in image #%d', i));
        imagesc(currImage); colormap gray; axis image;
        clim([0 5000]);
        title(sprintf('Click on same structure in image #%d', i));
        h = drawpoint('Color', 'cyan');
        pos = h.Position;
        xCurr = pos(1); yCurr = pos(2);
        close;

        RegistrationPoints(i).ImageName = ImageStruct.Name{i};
        RegistrationPoints(i).X = xCurr;
        RegistrationPoints(i).Y = yCurr;

        dx = xRef - xCurr;
        dy = yRef - yCurr;
        Transforms{i} = [dx, dy];

        tform = affine2d([1 0 0; 0 1 0; dx dy 1]);
        Rfixed = imref2d(size(referenceImage));
        ImageStruct.Registered{i} = imwarp(currImage, tform, 'OutputView', Rfixed);
    end
end

function ImageStruct = ApplyTransforms(ImageStruct, Transforms, referenceSize)
    ImageStruct.Registered = cell(1, numel(ImageStruct.Data));
    for i = 1:numel(ImageStruct.Data)
        dx = Transforms{i}(1);
        dy = Transforms{i}(2);
        tform = affine2d([1 0 0; 0 1 0; dx dy 1]);
        Rfixed = imref2d(referenceSize);
        ImageStruct.Registered{i} = imwarp(ImageStruct.Data{i}, tform, 'OutputView', Rfixed);
    end
end

function [ROI, ReferenceImage] = DefineSingleROIOnRegistered(ImageStruct, ImageType)
    tiledlayout('flow');
    for i = 1:numel(ImageStruct.(ImageType))
        ax = nexttile;
        imagesc(ax, ImageStruct.(ImageType){i}); colormap gray; axis image off;
        clim([0 5000]);
        title(ax, ImageStruct.Name{i}, 'Interpreter', 'none');
    end
    disp("Click on the image you'd like to draw the ROI on...");
    [~,~] = ginput(1);
    ax = gca;
    ReferenceImage = 1;  % simplified
    disp("Draw your ROI...");
    ROI.UserDefined = drawfreehand('Parent', ax, 'FaceAlpha', 0.2, 'Color', [.1 .8 .8]);
    ROI.Mask = createMask(ROI.UserDefined);
end

function ROIArray = DefineMultipleROIsOnRegistered(ImageStruct, ImageType)
    tiledlayout('flow');
    for i = 1:numel(ImageStruct.(ImageType))
        ax = nexttile;
        imagesc(ax, ImageStruct.(ImageType){i}); colormap gray; axis image off;
        clim([0 20000]);
        title(ax, ImageStruct.Name{i}, 'Interpreter', 'none');
    end
    ROIArray = [];
    done = false;
    while ~done
        disp("Click the image to draw a new ROI or press Enter in figure to finish.");
        try
            [~,~] = ginput(1);
        catch
            break;
        end
        ax = gca;
        h = drawfreehand('Parent', ax, 'FaceAlpha', 0.2, 'Color', 'm');
        mask = createMask(h);
        name = inputdlg('Enter a name for this ROI:', 'ROI Name', 1, {sprintf('ROI_%d', numel(ROIArray)+1)});
        if isempty(name), break; end
        ROIArray(end+1).Mask = mask;
        ROIArray(end).Name = name{1};
        ROIArray(end).Coordinates = h.Position;
    end
end

function OverlayROIOnImages(ImageStruct, mask, color, titleStr)
    figure('Name', titleStr);
    tiledlayout('flow');
    boundary = bwboundaries(mask);
    for i = 1:numel(ImageStruct.Registered)
        ax = nexttile;
        imagesc(ax, ImageStruct.Registered{i}); colormap gray; axis image off;
        clim([0 5000]);
        title(ax, ImageStruct.Name{i}, 'Interpreter', 'none');
        hold on;
        for b = 1:length(boundary)
            plot(ax, boundary{b}(:,2), boundary{b}(:,1), 'Color', color, 'LineWidth', 2);
        end
    end
end

function OverlayMultipleROIs(ImageStruct, ROIArray, titleStr)
    figure('Name', titleStr);
    tiledlayout('flow');
    for i = 1:numel(ImageStruct.Registered)
        ax = nexttile;
        imagesc(ax, ImageStruct.Registered{i}); colormap gray; axis image off;
        clim([0 20000]);
        title(ax, ImageStruct.Name{i}, 'Interpreter', 'none');
        hold on;
        for r = 1:numel(ROIArray)
            boundary = bwboundaries(ROIArray(r).Mask);
            for b = 1:length(boundary)
                plot(ax, boundary{b}(:,2), boundary{b}(:,1), 'm', 'LineWidth', 2);
            end
        end
    end
end

function outFilename = uniqueFilename(folder, baseName, ext)
    % Generate a unique filename by appending _1, _2, ... if needed
    % folder: directory where file will be saved
    % baseName: desired filename without extension
    % ext: extension with dot, e.g., '.xlsx'
    
    outFilename = fullfile(folder, [baseName, ext]);
    counter = 1;
    while exist(outFilename, 'file')
        outFilename = fullfile(folder, sprintf('%s_%d%s', baseName, counter, ext));
        counter = counter + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%
%% Animation Exporter %%
%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate GIF frame delay
framedelay = animationtime/length(im);

% Process all captured frames
for i = 1:length(im)

    % GIF does not support 3-D frame data. rgb2ind converts to indexed
    % image A and a colormap map
    [A, map] = rgb2ind(im{i}, 256);

    % Overwrite file if first step, otherwise append with delay
    if i == 1
        imwrite(A, map, animationfile, WriteMode="overwrite", LoopCount=Inf, DelayTime=framedelay)
    else
        imwrite(A, map, animationfile, WriteMode="append", DelayTime=framedelay)
    end
end
%% Wei's 90° Stereo Vision for Cerana/Hornet Interaction
% size:  720  1280  3  240fps      4,17 - 4,18 ms / frame
% 80 min auf 32 GB card
% side cam: wifi name ends on 26  name side26
% top cam: wifi ...88, name top88

%% Config Variables
load('stereoParamWEI.mat')      % camera config from stereocal (with checkerboard paper)
 videoFileLeft = 'Top88.mp4';   % Video from top
vvideoFileRight = 'Side26.mp4';  % Video from the side
      duration = 10000;         % frames you want to analyse (~ 2000 wasted by calibration)
                                % calculate back by multiplying seconds * 240 and round down
  LeftCam_left = 290;           % Amount of Pixels to cut from the left side of the top cam
 LeftCam_right = 940;           % Amount of Pixels to cut from the right side of the top cam
 RightCam_left = 365;           % Amount of Pixels to cut from the left side of the side cam
RightCam_right = 960;           % Amount of Pixels to cut from the right side of the side cam
  triggerLeft1 = 40;            % 4 corners were the sync-Light is positioned, can be rather big
  triggerLeft2 = 70;
  triggerLeft3 = 300;
  triggerLeft4 = 340;
 triggerRight1 = 570;           % 4 corners for the side cam, the other, upper is the top cam
 triggerRight2 = 600;
 triggerRight3 = 360;
 triggerRight4 = 420;
       MinBlob = 6;             % smalest group of Pixels still a bee
 shortestTrack = 40;            % least amound of coordinates that is a continuus track

%% Initiation Stuff
 tic
 profile on
  stereoParams = stereoParamsWEI;
showExtrinsics(stereoParams);
    readerLeft = vision.VideoFileReader(videoFileLeft, 'VideoOutputDataType', 'uint8');
   readerRight = vision.VideoFileReader(videoFileRight, 'VideoOutputDataType', 'uint8');
        player = vision.DeployableVideoPlayer('Location', [20, 400]);

%% were my trigger at?
 triggerLeft = zeros(1000,1);
triggerRight = zeros(1000,1);
for i=1:1000
          frameLeft = readerLeft.step();
         frameRight = readerRight.step();
     triggerLeft(i) = mean(mean(frameLeft(triggerLeft1:triggerLeft2,triggerLeft3:triggerLeft4,1)));
    triggerRight(i) = mean(mean(frameRight(triggerRight1:triggerRight2,triggerRight3:triggerRight4,1)));
    
end
 normed_triggerLeft = triggerLeft(1:end-1)-triggerLeft(2:end);
normed_triggerRight = triggerRight(1:end-1)-triggerRight(2:end);
       leftStartLED = find(normed_triggerLeft>10,1);   % this frame the LED turns on top cam,
      RightStartLED = find(normed_triggerRight>10,1);  % this frame the LED turns on side cam,

%% Sync the Frames
for i=1:200+leftStartLED-RightStartLED
    readerLeft.step();
end
for i=1:200
   readerRight.step();
end

%% Get average Background
                j = 1;
        frameLeft = readerLeft.step();
       frameRight = readerRight.step();
   frameLeftGray  = rgb2gray(frameLeft);
   frameRightGray = rgb2gray(frameRight);
 frameLeftGrayCut = frameLeftGray(:,LeftCam_left:LeftCam_right);
frameRightGrayCut = frameRightGray(:,RightCam_left:RightCam_right);
   bufferPicRight = double(frameRightGrayCut);
    bufferPicLeft = double(frameLeftGrayCut);
for i=10:10:1000  % actually 1000 till 2000
            frameLeft = readerLeft.step();
           frameRight = readerRight.step();
       frameLeftGray  = rgb2gray(frameLeft);
       frameRightGray = rgb2gray(frameRight);
     frameLeftGrayCut = frameLeftGray(:,LeftCam_left:LeftCam_right);
    frameRightGrayCut = frameRightGray(:,RightCam_left:RightCam_right);
        bufferPicLeft = bufferPicLeft+double(frameLeftGrayCut);
       bufferPicRight = bufferPicRight+double(frameRightGrayCut);
                    j = j+1;
end
 bufferPicLeft = uint8(round(bufferPicLeft/100));    % background 4 2 be substracted
bufferPicRight = uint8(round(bufferPicRight/100));
            kk = 1;


%% Substract BG and Track
               j = 1;
   centroidsLeft = nan(duration,2,20);
  centroidsRight = nan(duration,2,20);
        AreaLeft = nan(duration,1,20);
       AreaRight = nan(duration,1,20);
 OrientationLeft = nan(duration,1,20);
OrientationRight = nan(duration,1,20);
      LeftBWList = nan(duration,1);
     RightBWList = nan(duration,1);
for i=1:duration
              frameLeft = readerLeft.step();
             frameRight = readerRight.step();
       frameLeftGrayCut = frameLeft(:,LeftCam_left:LeftCam_right,1);
      frameRightGrayCut = frameRight(:,RightCam_left:RightCam_right,1);
     frameLeftGrayCutBG = bufferPicLeft-frameLeftGrayCut;
    frameRightGrayCutBG = bufferPicRight-frameRightGrayCut;
                 LeftBW = imbinarize(frameLeftGrayCutBG);
                RightBW = imbinarize(frameRightGrayCutBG);

    if mean(mean(LeftBW)) < 0.001
        tracksLeft = regionprops(LeftBW,'centroid','Area','Orientation');
        centroidsLeft(j,:,1:size(cat(1,tracksLeft.Centroid),1)) = cat(1, tracksLeft.Centroid)';
                 AreaLeft(j,1,1:size(cat(1,tracksLeft.Area),1)) = cat(1, tracksLeft.Area)';
   OrientationLeft(j,1,1:size(cat(1,tracksLeft.Orientation),1)) = cat(1, tracksLeft.Orientation)';
    end
    if mean(mean(RightBW)) < 0.001
        tracksRight = regionprops(RightBW,'centroid','Area','Orientation');
        centroidsRight(j,:,1:size(cat(1,tracksRight.Centroid),1)) = cat(1, tracksRight.Centroid)';
                 AreaRight(j,1,1:size(cat(1,tracksRight.Area),1)) = cat(1, tracksRight.Area)';
   OrientationRight(j,1,1:size(cat(1,tracksRight.Orientation),1)) = cat(1, tracksRight.Orientation)';
    end
     LeftBWList(j) = mean(mean(LeftBW));
    RightBWList(j) = mean(mean(RightBW));
    j = j+1;
end

figure              % both plots must be under 0.001 meaning the B/W pic is mostly dark 
subplot(2,1,1)      %- all white is supposted to be bees
plot(LeftBWList)
subplot(2,1,2)
plot(RightBWList)

%% Dataforming (boring)
                       centroidsRightCopy = centroidsRight;
                          AreaRightCutter = [AreaRight AreaRight];
centroidsRight(AreaRightCutter < MinBlob) = nan;
           AreaRight(AreaRight < MinBlob) = nan;
    OrientationRight(AreaRight < MinBlob) = nan;
                           MaxTracksRight = size(AreaRight,3);
for i = 1:MaxTracksRight
if sum(isnan(centroidsRight(:,1,i))) == duration
    centroidsRight = centroidsRight(:,:,1:i-1);
    MaxTracksRight = i-1;
  break;
end
end
                      centroidsLeftCopy = centroidsLeft;
                         AreaLeftCutter = [AreaLeft AreaLeft];
centroidsLeft(AreaLeftCutter < MinBlob) = nan;
           AreaLeft(AreaLeft < MinBlob) = nan;
    OrientationLeft(AreaLeft < MinBlob) = nan;
                          MaxTracksLeft = size(AreaLeft,3);
for i = 1:MaxTracksLeft
if sum(isnan(centroidsLeft(:,1,i))) == duration
    centroidsLeft = centroidsLeft(:,:,1:i-1);
    MaxTracksLeft = i-1;
    break;
end
end


%% Sort Tracks to Trajectories
for cameraSide = 1:2
    if cameraSide == 1
         testData = centroidsRight;
        MaxTracks = MaxTracksRight+200;
         testData = cat(3,testData, nan(duration,2,200));
    else
         testData = centroidsLeft;
        MaxTracks = MaxTracksLeft+200;
         testData = cat(3,testData, nan(duration,2,200));
 sortedTrackRight = testData2;
    end
    
           testData2 = nan(duration,2,MaxTracks);
    testData2(1,:,:) = testData(1,:,:);
            testData = testData(2:end,:,:);
    
    for i = 1:duration-1
            [n,d] = knnsearch(squeeze(testData(i,:,:))',squeeze(testData2(i,:,:))','k',1,'distance','minkowski','p',5);
      n(isnan(d)) = nan;
        testData2(testData2==0)=nan;
        for j = 1:MaxTracks % sum(isfinite(n))%    MaxTracks-sum(isnan(n))   do it as often as there are coordinates
            if length(n)>j
            if isfinite(n(j))     % if this one is                          %
                testData2(i+1,:,j) = testData(i,:,n(j)); % copy nearest match into sorted array
                testData(i,:,n(j)) = nan;                % delete this data in the old array
            end
            end
        end 
        
        if     sum(sum(isfinite(testData(i,:,:)))) > 0
            if sum(sum(isfinite(testData(i,:,:)))) == 2
                             buff = (testData(i,:,:));
                buff(isnan(buff)) = [];
                if sum(~isnan(squeeze(testData2(i+1,1,:)))) == 0
                    testData2(i+1,:,1) = buff;
                else
                    testData2(i+1,:,find(isfinite(squeeze(testData2(i+1,1,:))),1,'last')+1) = buff;
                end
                testData(i,:,:) = nan;
            else
                buffX = squeeze((testData(i,1,:)));
                buffY = squeeze((testData(i,2,:)));
                buffX(isnan(buffX)) = [];
                buffY(isnan(buffY)) = [];
                buff = [buffX buffY]';
                if sum(~isnan(squeeze(testData2(i+1,1,:)))) == 0
                    testData2(i+1,:,1:length(buffY)) = buff;
                else
                    testData2(i+1,:,find(isfinite(squeeze(testData2(i+1,1,:))),1,'last')+1:find(isfinite(squeeze(testData2(i+1,1,:))),1,'last')+length(buffY)) = buff;
                end
                testData(i,:,:) = nan;
            end
        end
    end
end
sortedTrackLeft = testData2;


%% cut out the empty nan-rows at the end
for i = 1:size(sortedTrackLeft,3)
if sum(isnan(sortedTrackLeft(:,1,size(sortedTrackLeft,3)-i))) ~= duration
    sortedTrackLeft = sortedTrackLeft(:,:,1:i-1);
  break;
end
end

for i = 1:size(sortedTrackRight,3)
if sum(isnan(sortedTrackRight(:,1,size(sortedTrackRight,3)-i))) ~= duration
    sortedTrackRight = sortedTrackRight(:,:,1:i-1);
  break;
end
end


%% camera cleaning from stereoParams
sortedTrackRightCopy = sortedTrackRight;
 sortedTrackLeftCopy = sortedTrackLeft;
   undistortedPoints = nan(size(sortedTrackRight));
for cameraSide = 1:2 % first Right, then Left
    if cameraSide == 2
        undistortedPointsRight = undistortedPoints;
              sortedTrackRight = sortedTrackLeft;
    end
    for i = 1:duration
        sortedTrackRight1x = squeeze(sortedTrackRight(i,1,:));
        sortedTrackRight1y = squeeze(sortedTrackRight(i,2,:));
        sortedTrackRight1x(isnan(sortedTrackRight1x(:,1))) = [];
        sortedTrackRight1y(isnan(sortedTrackRight1y(:,1))) = [];
        sortedTrackRight1 = cat(2,sortedTrackRight1x,sortedTrackRight1y);
        if ~isempty(sortedTrackRight1x)
            if size(sortedTrackRight1x,1) == 1
                undistortedPoints(i,:,1) = undistortPoints(sortedTrackRight1,stereoParams.CameraParameters1);
            else
                undistortedPoints(i,:,1:length(sortedTrackRight1x)) = undistortPoints(sortedTrackRight1,stereoParams.CameraParameters1)';
            end
        else
            undistortedPoints(i,:,:) = nan;
        end
        
    end
end
     undistortedPointsLeft = undistortedPoints;
          sortedTrackRight = sortedTrackRightCopy;
undistortedPointsRightCopy = undistortedPointsRight;
 undistortedPointsLeftCopy = undistortedPointsLeft;

                undistortedPointsRight = undistortedPointsRightCopy;
undistortedPointsRight(end-10:end,:,:) = 0;
                 undistortedPointsLeft = undistortedPointsLeftCopy;
 undistortedPointsLeft(end-10:end,:,:) = 0;

%%###############
%%von hier auf xyz erst dann tracks indexen

undistortedPointsRight = undistortedPointsRightCopy;
 undistortedPointsLeft = undistortedPointsLeftCopy;

xyz = nan(duration,3,size(undistortedPointsLeft,3));

for i = 1:size(undistortedPointsLeft,1)
    for j = 1:size(undistortedPointsLeft,3)
        [c, index] = min(abs(undistortedPointsRight(i,1,:)-undistortedPointsLeft(i,1,j)));
        if abs(c) < 60
                        xyz(i,1:2,j) = undistortedPointsLeft(i,1:2,j);
                          xyz(i,3,j) = undistortedPointsRight(i,2,index);
        undistortedPointsLeft(i,:,j) = nan;
   undistortedPointsRight(i,:,index) = nan;
        end
    end
end

figure   % 3D plot of trajectories
hold on
for i = 1:size(xyz,3)
plot3(xyz(:,1,i),xyz(:,2,i),-xyz(:,3,i),'*')
end
%%###############


for cameraSide = 1:2 % first Right, then Left
    if cameraSide == 2
        undistortedPointsLeft = undistortedPointsRight;
        SuperSortLeftCopy = SuperSortLeft;
    end
%
tracksLeft = nan(10000,3,200);
         n = 1;
undistortedPointsLeft(undistortedPointsLeft == 0) = nan;
for k = 1:size(undistortedPointsLeft,3)
    indexStart = 1;
    indexStop = 1;
    while undistortedPointsLeft(indexStart,1,1) ~= 0 && undistortedPointsLeft(indexStop,1,1) ~= 0
        indexStart = find(isfinite(undistortedPointsLeft(indexStop+1:end,1,k)),1)+indexStop;
        if isempty(indexStart)
            break
        end
        indexStop = find(isnan(undistortedPointsLeft(indexStart:end,1,k)),1)+indexStart-2;
        if isempty(indexStop)
            break
        end
        if -indexStart+indexStop > shortestTrack
            tracksLeft(1:-indexStart+indexStop+1,1:2,n) = undistortedPointsLeft(indexStart:indexStop,:,k);
            tracksLeft(1:-indexStart+indexStop+1,3,n) = [indexStart:1:indexStop];
            undistortedPointsLeft(indexStart:indexStop,:,k) = nan;
            n = n+1;
        end
    end
end

tracksLeft(tracksLeft == 0) = nan;

for i = 1:size(tracksLeft,1)-1
    if sum(isnan(tracksLeft(i,1,:))) == size(tracksLeft,3)
        tracksLeft = tracksLeft(1:i-1,:,:);
        break;
    end
end
for i = 1:size(tracksLeft,3)-1
    if sum(isnan(tracksLeft(:,1,size(tracksLeft,3)-i))) == size(tracksLeft,1)
        tracksLeft = tracksLeft(:,:,1:size(tracksLeft,3)-(i-1));
        break;
    end
end

tracksLeftCopy = tracksLeft;
tracksLeft = tracksLeftCopy;

SuperSortLeft = nan(size(tracksLeft));
            n = 1;
     minDelta = 30;
   minIndexes = 10;
   for k = 1:size(tracksLeft,3)
       deltaTrack = tracksLeft(1:end-1,1,k)-tracksLeft(2:end,1,k);
       deltaTrack(1) = 100;
       deltaTrack(isnan(deltaTrack)) = 0;
       deltaTrack(find(deltaTrack ~= 0, 1, 'last')) = 100;
       deltaTrack(abs(deltaTrack)<minDelta) = 0;
       deltaTrack(deltaTrack~=0) = 1;
       indexes = find(deltaTrack == 1);
       if length(indexes) > 1
           for i = 1:length(indexes)-1
               if indexes(i+1)-indexes(i) > minIndexes
                   SuperSortLeft(1:indexes(i+1)-indexes(i)+1,:,n) = tracksLeft(indexes(i):indexes(i+1),:,k);
                   tracksLeft(indexes(i):indexes(i+1),:,k) = nan;
                   n = n+1;
               end
           end
       end
   end

SuperSortLeft(SuperSortLeft == 0) = nan;

for i = 1:size(SuperSortLeft,1)-1
    if sum(isnan(SuperSortLeft(i,1,:))) == size(SuperSortLeft,3)
        SuperSortLeft = SuperSortLeft(1:i-1,:,:);
        break;
    end
end
for i = 1:size(SuperSortLeft,3)-1
    if sum(isnan(SuperSortLeft(:,1,size(SuperSortLeft,3)-i))) == size(SuperSortLeft,1)
        SuperSortLeft = SuperSortLeft(:,:,1:size(SuperSortLeft,3)-(i-1));
        break;
    end
end

end
SuperSortRight = SuperSortLeft;
SuperSortLeft = SuperSortLeftCopy;

%#################################################

twistedLeft = SuperSortLeft;   % y-achsis upsidedown, twisted, corrected
twistedLeft(:,2,:) = (SuperSortLeft(:,2,:)-800)*(-1);
twistedRight = SuperSortRight;
twistedRight(:,2,:) = (SuperSortRight(:,2,:)-800)*(-1);

xyz = nan(duration,3,size(twistedLeft,3));

for i = 1:size(twistedLeft,1)
    for j = 1:size(twistedLeft,3)
        [c, index] = min(abs(twistedRight(i,1,:)-twistedLeft(i,1,j)));
        if abs(c) < 60
                       xyz(i,1:2,j) = twistedLeft(i,1:2,j);
                         xyz(i,3,j) = twistedRight(i,2,index);
                 twistedLeft(i,:,j) = nan;
            twistedRight(i,:,index) = nan;
        end
    end
end

toc
profile off


%% Ploterino
figure  % Left unsorted (Left, Top88, cam2) 
%imshow(LeftBW)
hold on
for i=1:size(centroidsLeftCopy,3)
    plot(centroidsLeftCopy(:,1,i),centroidsLeftCopy(:,2,i),'*')
end
hold off

figure  
hold on
for i=1:size(sortedTrackLeft,3)
    plot(sortedTrackLeft(:,1,i),sortedTrackLeft(:,2,i),'*')
end
hold off

figure  % left sorted
imshow(RightBW)
hold on
for i=1:size(centroidsRightCopy,3)
    plot(centroidsRightCopy(:,1,i),centroidsRightCopy(:,2,i),'*')
end
hold off

figure  % right sorted
%imshow(RightBW)
hold on
for i=1:size(sortedTrackRightCopy,3)
    plot(sortedTrackRightCopy(:,1,i),-sortedTrackRightCopy(:,2,i),'*')
end
hold off

figure  % Right cam distrotion fixed in blue
hold on
imshow(RightBW)
for i=1:size(sortedTrackRight,3)
%plot(undistortedPointsRight(:,1,i),undistortedPointsRight(:,2,i),'b*')
plot(sortedTrackRight(:,1,i),sortedTrackRight(:,2,i),'ro')
end
hold off

figure    % Left cam distrotion fixed in blue
hold on
imshow(LeftBW)
hold on
for i=1:size(sortedTrackLeft,3)
%plot(undistortedPointsLeft(:,1,i),undistortedPointsLeft(:,2,i),'*')
plot(sortedTrackLeft(:,1,i),sortedTrackLeft(:,2,i),'*')
end
hold off


figure  
hold on
for i=1:size(SuperSortLeftCopy,3)
  %  if sum(isfinite(SuperSortLeft(:,1,i)))>1000
    plot(SuperSortLeftCopy(:,1,i),SuperSortLeftCopy(:,2,i),'*')
   % end
end
hold off


figure  
hold on
for i=1:size(tracksLeftCopy,3)
    plot(tracksLeftCopy(:,1,i),tracksLeftCopy(:,2,i),'*')
end
hold off

      
      
figure  
hold on
for i=1:size(undistortedPointsLeft,3)
    plot(undistortedPointsLeft(:,1,i),undistortedPointsLeft(:,2,i),'*')
end
hold off  
      


figure   % 3D plot of trajectories
hold on
for i=1:size(xyz,3)
plot3(xyz(:,1,i),xyz(:,2,i),-xyz(:,3,i),'*')
end
%groundPlane
[X,Y] = meshgrid(-250:100:650,-250:100:650);
Z=zeros(size(X))-650;
surface(X,Y,Z,'FaceColor','none')
view(3)
% backplate
[X,Y] = meshgrid(-250:100:650,-650:100:250);
Z=zeros(size(X))+650;
surface(X,Z,Y,'FaceColor','none')
view(3)




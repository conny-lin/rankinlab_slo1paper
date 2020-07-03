

%% define paths
% program path
pProgram = '/Users/connylin/OneDrive/MATLAB/Functions_Developer';
pJava = '/Users/connylin/OneDrive/MATLAB/Java_Programs';
% add program path
addpath(pProgram);
% data path
pData = '/Volumes/ParahippocampalGyrus/MWT_Data';
% output
pO = '/Users/connylin/OneDrive/Dance_Output';


%% choose a set of data to play
ExpName = '20141116B_JS_300s0x0s0s_24hrsPE_nahMT14480';
% ExpName = '20141106C_SM_100s30x10s10s';
pExp = [pData,'/',ExpName];


%% USER INPUTS
display ' ';
display 'Name your analysis output folder';
AnalysisName = input(': ','s');


%% GENERATE CHOR JAVA SYNTAX (chorescript)

% GENERATE JAVA ARGUMENTS
% path to java programs
%javapath = [strrep(userpath,pathsep,''),'/MATLAB MWT/SubFun_Java'];
% javapath = [pProgram,'/Java'];
javapath = pJava;
choreversion = 'Chore_1.3.0.r1035.jar';

b = blanks(1); % blank
% call java 
javacall = 'java -jar'; javaRAM = '-Xmx7G'; javaRAM7G = '-Xmx7G';
beethoven = ['''',javapath,'/Beethoven_v2.jar','''']; % call beethoven 
chor = ['''',javapath,'/',choreversion,'''']; % call chor 
% chor calls 
map = '--map';
% settings 
pixelsize = '-p 0.027'; speed = '-s 0.1'; 
mintime = '-t 20'; minmove = '-M 2'; shape = '--shadowless -S';
% plugins 
preoutline = '--plugin Reoutline::exp';  
prespine = '--plugin Respine';
% plugins (reversals) 
revbeethoven_trv = '--plugin MeasureReversal::tap::dt=1::collect=0.5::postfix=trv';
revignortap_sprevs = '--plugin MeasureReversal::postfix=sprevs';
rev_ssr = '--plugin MeasureReversal::tap::collect=0.5::postfix=ssr';

% dat output 
odrunkposture = '-O drunkposture -o nNslwakb';
odrunkposture2 = '-O drunkposture2 -o nNslwakbcemM';

oconny = '-o 1nee#e*ss#s*SS#S*ll#l*LL#L*ww#w*aa#a*mm#m*MM#M*kk#k*bb#b*pp#p*dd#d'; % Conny's 
obeethoven = '-o nNss*b12M'; % standard for Beethoven
oshanespark = '-O shanespark -o nNss*b12M'; % standard for Beethoven
oevan = '-O evan -o nNss*b12'; % Evan's dat output
oevanall = '-O evanall -N all -o nNss*b12';
oswanlakeall = '-O swanlakeall -N all -o tnNemMawlkcspbd1';
oswanlake = '-O swanlake -o tnNemMawlkcspbd1e#m#M#a#w#l#k#c#s#p#b#d#e-m-M-a-w-l-k-c-s-p-b-d-e*m*M*a*w*l*kvc*s*p*b*d*';
ostarfish = '-O starfish -N all -o nNss*b12xyMmeSakcr';

chorscript{1} = [javacall,b,javaRAM,b,chor,b,pixelsize,b,speed,b,...
            mintime,b,minmove,b,shape,b,ostarfish,b,...
            preoutline,b,prespine,b]; 
fval = {'*starfish*'};


%% Choreography
% get all MWT folder paths
[~,pMWT] = dircontentmwt(pExp);

% validate if need to chor
needChorVal = true(size(pMWT));
for x = 1:numel(pMWT);
    % this only checks the first fval
    if isempty(dircontent(pMWT{x},fval{1})) == 0 
        needChorVal(x) = false;
    end
end
pMWTc = pMWT(needChorVal);

% chor
if isempty(pMWTc) == 0
display(sprintf('Need to chor %d MWT files',numel(pMWTc)));
str = 'Chor-ing MWTfile [%s]...';
for x = 1:numel(pMWTc); 
    [~,fn] = fileparts(pMWTc{x}); file = strcat('''',pMWTc{x},''''); 
    display(sprintf(str,fn));
    for x = 1:numel(chorscript) 
        system([chorscript{x} file], '-echo'); 
    end  
end
end
display 'Chor Completed';  


%% create output folder
display ' ';
display 'Name your analysis output folder';
name = ['StarFish_',generatetimestamp,'_',AnalysisName];
% cd(pO);
pOFolder = [pO,'/',name];
if isdir(pOFolder) == 0; mkdir(pOFolder); end


%% STARFISH
% For Matlab set current folder to the time-date stamped experiment folder,
% then use this script (start and finish variables define the time window 
% of interest):

% DEFINE START-FINISH TIME SERIES
startTime = [40 [100:10:290]];
finishTime = [60 [110:10:300]];

% LOAD DATA TO MATLAB (EFFICIENCY REASON)
display ' ';
display 'converting raw data...';
for MWTfnX = 1:numel(pMWT)
    cd(pMWT{MWTfnX}); % cd to MWT folder
    [~,a] = fileparts(pMWT{MWTfnX});
    display(sprintf('converting [%s]',a));
    dirData = dir('*starfish*');
    for i = 1:numel(dirData);
        DataImport{MWTfnX,i} = dlmread(dirData(i).name);
    end
end

% save Data
% cd(pOFolder);
% save('matlab.mat','DataImport');


%% Making graphs
display ' ';
display 'making starfish graphs';

for t = 1:numel(startTime)
    
    for MWTfnX = 1:numel(pMWT)
        [p,MWTfn] = fileparts(pMWT{MWTfnX}); % get MWT file name
        
        % create output sub folder
        [~,GroupName] = fileparts(p);
        % make group folder if don't have one
        pOGroup = [pOFolder,'/',GroupName];
        if isdir(pOGroup) == 0; mkdir(pOGroup); end

        %% STAR FISH CODE -------------------------------------------------
%         displaceTog = [0];
        start = startTime(t);
        finish = finishTime(t);
        figure1 = figure('Visible','off');
        hold on
        k=1;
        for i = 1:sum((cellfun(@numel,DataImport(MWTfnX,:)) >0));
            storedData = DataImport{MWTfnX,i};
            startime = storedData(1,1);
            endtime = storedData(end,1);
            time = storedData(:,1);
            xpos = storedData(:,9);
            xpos1 = xpos(1);
            ypos = storedData(:,10);
            ypos1 = ypos(1);
            
            if startime < start && endtime > finish;
                k=1+k;
                irTimes = time<start|time>finish;
                storedData(irTimes,:)=[];
                conx = isnan(xpos); 
                storedData(conx,:)=[];
                
                cony = isnan(ypos); 
                storedData(cony,:)=[];

                x = (xpos(1)-xpos);
                y = (ypos(1)-ypos);
                
                randCol = [rand rand rand]; % random color generator;
                plot(x,y,'color',randCol,'LineWidth',1.5);

%                 displace = sqrt(x(end)*x(end) + y(end)*y(end));
%                 displaceTog = [displaceTog;displace];

            end
            clear conx cony irTimes storedData displace
        end
        %% STAR FISH CODE ENDS --------------------------------------------


    % save figure
    figureName = [MWTfn,'[',num2str(startTime(t)),'_', num2str(finishTime(t)),']'];
    savefigepsOnly(figureName,pOGroup);
    end

end

display 'Starfish done';


%%  Chrysanthenum
% plot starfish with line has speed as color
% tap as black crosses or black circle
% X = each row is time(i), col1 = position at t(i-1), col2 = position at t(i)
%% speed color table
% create from yellow to red
maxSpeed = 0.6;
a = (1:0.05:0)'; % 20 colors
colorScale = [ones(size(a)), a zeros(size(a))]; % yellow to red
step = 0.6/numel(a);
speedScale = (0:step:maxSpeed)';


%%

X = [xpos(1:end-1) xpos(2:end)];
Y = [ypos(1:end-1) ypos(2:end)];
S = speed(2:end);
tap = tap(2:end);
% translate speed
SC = nan(numel(S),3);
SC(S > max(SpeedScale)) = colorScale(end,:);
SC(S == 0) = colorScale(1,:);
for i = 1:numel(speedScale)-1
    SC(S>speedScale(i-1) & S<=speedScale(i)) = colorScale(i,:);
end

%% plot line
fig1 = figure('Visible','on'); 
hold on;
for i =  1:size(X,1)
    line(X(i,:),Y(i,:),SC(i,:));
end

%% plot tap
xtap = X(tap==1,2);
ytap = Y(tap==1,2);
plot(xtap,ytap,'k','LineStyle','none');















    









































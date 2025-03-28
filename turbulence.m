function [output] = turbulence(SUB, TS, TR, CoG_file, outdir)

%% Function that computes the Turbulence-measures for a given subject
% Input data:
% SUB:  String of Subject-ID
% TS:   Timeseries as matrix with [nParcels, nTimepoints]
% TR:   TR of resting state scan

%% Variables setting
%Utils from TTB toolbox. Demean function is used during Hilbert transform
%cd('C:/Users/bened/OneDrive - Fondazione Istituto Italiano Tecnologia/Desktop/University/PhD IIT/Turbulence_Autisms');

disp('Running turbulence analysis');

% COG: Center of gravity of all 1000 schaefer parcels, Matrix [1000,3] with
% parcels as rows and coordinates for cog of each parcel as columns
% CoG = load(fullfile(root, 'postproc/Files/schaefercog.mat'));
CoG = load(CoG_file);
CoG = CoG.Schaefer1054_CoG_updated;

% RSN: Information about which parcel belongs to which resting state
% Define label groups and their sizes (info from .txt file of Schaefer1000
% parcellation of 7networks_2mm_MNI
group_sizes = [81, 91, 61, 55, 29, 57, 126];
hemisphere_labels = repelem(1:7, group_sizes)';
labels = [repmat(hemisphere_labels, 2, 1); repelem(8, 54)'];
RSNs = ["Visual", "Somatomotor", "Dorsal_Attention", "Saliency_Ventral_Attention", "Limbic", "Control", "Default_Mode", "Subcortical"];


%Subject specific parameters
subID = SUB; 
ts=TS; %Ts-Matrix for this subject [nParcels,nTimepoints]
tr = TR; %TR of resting state
NPARCELLS = length(ts(:,1)); %Number of nodes of the parcellation
Tmax = length(ts(1,:)); %Timepoints

% Turbulence specific parameters
NR = 400; %Hardcoded variabled from Turbu paper. Steps of spatial range.
%The following two variables are the upper and lower bound of where to look
%for the power law (in which range)
NRini = 20; %Lower bound of power law
NRfin = 80; %Upper bound of power law
iLambda = 0.01; %Lower bound of lambda
fLambda = 0.3;  %Upper bound of lambda
stepsLambda = 0.03; %Steps of Lambda
LAMBDA = iLambda:stepsLambda:fLambda; %Range of lambdas
LAMBDA = flip(LAMBDA); %The Lambda goes from high to low, so small to large distances
NLAMBDA= length(LAMBDA); %How many lambdas we wanna examine


% Bandpass filter settings
fnq=1/(2*tr);                 % Nyquist frequency
flp = 0.008;                  % lowpass frequency of filter (Hz)
fhi = 0.08;                   % highpass
Wn=[flp/fnq fhi/fnq];         % butterworth bandpass non-dimensional frequency
k=2;                          % 2nd order butterworth filter
[bfilt,afilt]=butter(k,Wn);   % construct the filter

%% Computations

% Calculating the eucledian distance between every node. Resulting matrix
% rr is a matrix which encapsulates the eucledian distance between each
% pair of nodes
disp('Computing euclidean distance');
for i=1:NPARCELLS
    
    for j=1:NPARCELLS
        
        rr(i,j)=norm(CoG(i,:)-CoG(j,:));

    end % for j=1:NPARCELLS

end % for i=1:NPARCELLS


% Here, we fix a scale on which we want to
% assess the spatial information flow. This is what delta does. We
% extract the maximum range in the brain, split it up (with the
% random number 400), and then for every step, see how well
% "correlated" the brain is (see information transfer).

range=max(max(rr)); %Maximum range/distance in the brain
delta=range/NR;     %Steps of analysis (in mm)
xrange = (delta/2):delta:range; %Vector of all the distances, starting from delta/2 to max_range in steps of delta

% Shrinkage factor between two parcels (as a function of lambda).
% The eucledian distance for each pair of parcels is inserted in
% the Kernel function. The "G" part of the equation
disp('Figuring out shrinkage factors');
C1=zeros(NLAMBDA,NPARCELLS,NPARCELLS);
ilam=1;

for lambda=LAMBDA

    for i=1:NPARCELLS

        for j=1:NPARCELLS

            C1(ilam,i,j)=exp(-lambda*rr(i,j));

        end %for j=1:NPARCELLS

    end % for i=1:NPARCELLS

    ilam=ilam+1;

end % for lambda=LAMBDA

%Preallocation
Turbulence=zeros(NLAMBDA,1);
Turbulence_Node=zeros(NPARCELLS, NLAMBDA);
Turbulence_RSN = zeros(length(1:max(labels)), NLAMBDA);
Information_Cascadeflow=zeros(NLAMBDA,1);
Information_Transfer = zeros(NLAMBDA,1);
fclam=zeros(NLAMBDA,NPARCELLS,NPARCELLS);
enstrophy=zeros(NLAMBDA,NPARCELLS,Tmax);
signal_filt=zeros(NPARCELLS,Tmax);
Phases=zeros(NPARCELLS,Tmax);

%Filtering, Hilbert-Transform and Kuramoto LOP for each parcel
disp('Filtering...');
for seed=1:NPARCELLS

    %Detrending is used to remove any long-term fluctuations
    ts(seed,:)=detrend(ts(seed,:)-nanmean(ts(seed,:)));

    %Bandpass filtering (except for the parcels, where there are the
    %nans). FiltFilt is a forward and backwards filter, removing phase
    %distortions, with the arguments from the butterworth before.
    if sum(isnan(ts(seed,:)))<1
        signal_filt(seed,:) =filtfilt(bfilt,afilt,ts(seed,:));
    else
        signal_filt(seed,:) = ts(seed,:);
    end % if sum(isnan(ts(seed,:)))<1

    %Calculating the analytic signal and extracting the phases
    Xanalytic = hilbert(demean(signal_filt(seed,:)));
    Phases(seed,:) = angle(Xanalytic);

end % for seed=1:NPARCELLS

%Calculating the Kuramoto local order parameter. The resulting
%matrix enstrophy1 is the final matrix on which all other analysis steps will
%be done [10,1000,180], with lambdas ar rows, parcels as columns
%and timepoints as the third dimension
disp('Computing Kuramoto...');
for i=1:NPARCELLS

    for ilam=1:NLAMBDA

        C1lam=squeeze(C1(ilam,:,:));
        sumphases=nansum(repmat(C1lam(i,:)',1,Tmax).*complex(cos(Phases),sin(Phases)))/nansum(C1lam(i,:));
        enstrophy1(ilam,i,:)=abs(sumphases);

    end % for ilam=1:NLAMBDA

end % for i=1:NPARCELLS

%% Here the analysis of the different measure begin.

%Calculating Turbulence for each RSN. For each lambda, this is done by
%choosing the parcels that correspond to the network and taking all of
%the values inside to calculate the sd of the KLOP of each parcel.
disp('Computing turbulence for each RSN...');
for ii=1:max(labels)
    
    clear enstrophy;

    for ilam=1:NLAMBDA

        enstrophy(ilam,:,:) = enstrophy1(ilam,labels(:)==ii,:);
        Turbulence_RSN(ii, ilam)=nanstd(squeeze(enstrophy(ilam,:)));

    end % for ilam=1:NLAMBDA

end % for ii=1:max(labels)

%Turbulence and Node-Turbulence
disp('Computing node turbulence over lambdas...');
for ilam=1:NLAMBDA

    Turbulence(ilam)=nanstd(squeeze(enstrophy1(ilam,:)));

    for j=1:NPARCELLS

        Turbulence_Node(j, ilam) = nanstd(squeeze(enstrophy1(ilam,j,:)));
    
    end % for j=1:NPARCELLS

end % for ilam=1:NLAMBDA

%Metastability
disp('Computing metastability...');
gKoP = nanmean(abs(nansum(complex(cos(Phases(:,:)),sin(Phases(:,:))),1))/NPARCELLS);
Meta = nanstd(abs(nansum(complex(cos(Phases(:,:)),sin(Phases(:,:))),1))/NPARCELLS);

%Information Cascade Flow
disp('Computing information cascade flow...');
Information_Cascadeflow(1)=NaN; %Because we have one lambda-pair less

%Basically the equation for cascade flow
for ilam=1:NLAMBDA-1

    [cc, pp]=corr(squeeze(enstrophy1(ilam+1,:,2:end))',squeeze(enstrophy1(ilam,:,1:end-1))');
    Information_Cascadeflow(ilam+1)=nanmean(abs(cc(pp(:)<0.05)));

end % for ilam=1:NLAMBDA-1

% Information-Transfer across space
% fclam is a 3D-matrix that for each lambda has a corrcoef for each
% pair of nodes inserted (symmetrical to diagonal)
disp('Computing information transfer...');
for ilam=1:NLAMBDA

    fclam(ilam,:,:)=corrcoef(squeeze(enstrophy1(ilam,:,:))');

end % for ilam=1:NLAMBDA

% So for each lambda we want to
% assess the spatial information flow. This is what delta does. We
% basically extract the maximum range in the brain, split it up (with the
% random number 400), and then for each of this step, see how well
% "correlated" the brain is (see information transfer).
for lam=1:NLAMBDA %For each lambda

    numind=zeros(1,NR);
    fcra=zeros(1,NR);

    for i=1:NPARCELLS

        for j=1:NPARCELLS

            r=rr(i,j); %Extracting the Eucledian distance for this pair
            index=floor(r/delta)+1; %This is still the distance, but scaled (so instead of looking at 1cm, we look at delta cm)
            
            if index==NR+1
                index=NR;
            end % if index==NR+1
            
            mcc=fclam(lam,i,j); %Extracting the corrcoeff
            
            if ~isnan(mcc)
                fcra(index)=fcra(index)+mcc; %For this specific distance, add the corrcoef
                numind(index)=numind(index)+1; %And count the number of distances that have been added
            end % if ~isnan(mcc)

        end % for j=1:NPARCELLS

    end % for i=1:NPARCELLS


    % Fitting the power law.
    grandcorrfcn=fcra./numind; %Calculating the mean corrcoef for each distance
    clear xcoor;
    clear ycoor;
    nn=1;

    %We are only interested in the range between 20*delta and 80*delta
    for k=NRini:NRfin

        if grandcorrfcn(k)>0

            xcoor(nn)=log(xrange(k));
            ycoor(nn)=log(grandcorrfcn(k)/grandcorrfcn(NRini));
            nn=nn+1;

        end % if grandcorrfcn(k)>0

    end % for k=NRini:NRfin

    linfunc = @(A, x)(A(1)*x+A(2));
    options=optimset('MaxFunEvals',10000,'MaxIter',1000,'Display','off');
    A0=[-1 1];
    [Afit, ~]= lsqcurvefit(linfunc,A0,xcoor,ycoor,[-4 -10],[4 10],options);
    Information_Transfer(lam)=abs(Afit(1));  %Extract the slope of that resulting graph

end % for lam=1:NLAMBDA

%Information cascade as the mean over all scales of cascade flow
Information_Cascade = nanmean(Information_Cascadeflow);

%% Saving the output files
% First as a Matlab structure
disp('Saving data...');
output.Turbulence=Turbulence;
output.Turbulence_Node=Turbulence_Node;
output.Turbulence_RSN=Turbulence_RSN;
output.Information_Transfer=Information_Transfer;
output.Information_Cascade=Information_Cascade;
output.Information_Cascadeflow=Information_Cascadeflow;
output.GlobalKuramoto=gKoP;
output.Metastability=Meta;
output.Lambda = LAMBDA; 
output.LocalKuramoto = enstrophy1;
output.Phases = Phases;
save(fullfile(outdir, sprintf('%s_turbulence_measures.mat',subID)), 'output');

% Then each variable as a csv file 
% Adding the header 
Turbulence_Header = [LAMBDA',Turbulence]; 
Turbulence_Header = [["Lambda", "Value"] ;Turbulence_Header];
Information_Transfer_Header = [LAMBDA',Information_Transfer]; 
Information_Transfer_Header = [["Lambda", "Value"] ;Information_Transfer_Header];
Information_Cascadeflow_Header = [LAMBDA',Information_Cascadeflow]; 
Information_Cascadeflow_Header = [["Lambda", "Value"] ;Information_Cascadeflow_Header];
Turbulence_Node_Header = [(1:NPARCELLS)',Turbulence_Node]; 
Turbulence_Node_Header = [["Lambda",LAMBDA];Turbulence_Node_Header];
Turbulence_RSN_Header = [RSNs',Turbulence_RSN];
Turbulence_RSN_Header = [["Lambda",LAMBDA];Turbulence_RSN_Header];

%Creating cell with all the measures with headers 
measures = {Turbulence_Header, Turbulence_Node_Header, Turbulence_RSN_Header, ...
    Information_Transfer_Header, Information_Cascade, ...
    Information_Cascadeflow_Header, gKoP, Meta};
fext = '.csv'; % file extension to use (e.g., .csv)
% filenames for the output files
fnames = {'turbulence', 'turbulence_node', ...
    'turbulence_rsn', 'turbulence_information_transfer', ...
    'turbulence_information_cascade', 'turbulence_information_cascadeflow', ...
    'turbulence_global_kuramoto', 'turbulence_metastability'};

% loop over measures and save output files
for ifile = 1:length(fnames)

    data2write = measures(ifile);
    outfile2use = fullfile(outdir, sprintf('%s_%s%s', subID, fnames{ifile}, fext));
    writematrix(data2write{1}, outfile2use);

end % for ifile

end % function turbulence


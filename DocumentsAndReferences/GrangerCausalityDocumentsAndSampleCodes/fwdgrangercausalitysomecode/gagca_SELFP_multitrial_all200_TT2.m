% N Killian, compute and save Granger causality values
% uses methods by Zu-Xiang Liu, which uses methods based on the  
% Autoregressive modeling toolbox of Stijn deWaele which is available on matlab exchange

% clar
% layersets = [1 2;1 3;2 3;];
% for layi = 1:size(layersets,1)
% nexdir     = 'S:\NEX files\Peepers\';%where to load data
% nexdir     = 'D:\Timmy\EC_vpc2\';%where to load data
nexdir     = 'D:\raw\';%where to load data
% nexdir = 'D:\My Dropbox\ECG2010\Data\';
nexfiletype   = ['TT*vpc1*SUA.nex'];%
% nexfiletype   = ['TT*vpc1*LFP.nex'];
matdir     = 'D:\AD\ECVPC1\';%where to load analyzed data from
% matdir = 'D:\My Dropbox\';
newmatdir = 'D:\My Dropbox\ECG2010\Data\';

% tag1 = '_CSDprepfixations200_enc';%FT processed/prepared data saved to .mat
% tag2 = '_CSDprepfixations200_rec';
% tag1 = '_LFP_TL_100_enc';
% tag1 = '_SUA_TL_SOn_nofilt_enc';
tag1 = '_SUA_TL_SOn_vtl_dft_CSD_enc';
tag2 = '_SUA_TL_SOn_vtl_dft_CSD_rec';

% tag1 = '_SUA_TL_SOn_nodft_CSD_enc';
% tag2 = '_SUA_TL_SOn_nodft_CSD_rec';
dataformat = 'timelock';
hionly = 0;
%for the doCSD_ script
dofilteredversions = 1;
goodonly = 1;contgoodonly = 0;
fixbad = 0;notip = 1;
%for the plotCSD_ script
savefigs = 1;
plotallELFP = 1;
plotedgeELFP = 1;
microns = 0;
CSDofAvg = 0;

structure      = 'EC';
% get data lists and parameters-------------------------------------------
[Num Txt Raw] = xlsread(['D:\My Dropbox\DB\ArrayInfo.xls']);
valid = Raw(find(strcmp(structure,Raw(2:end,8))&[Raw{2:end,11}]')+1,1);

leavein = 1;    keep = {'100216'};
pullout = 0;    remove = {'090320' '090730'};
dnex = dir(strcat(nexdir,nexfiletype));

count = 1;
for k=1:length(dnex)
    if pullout
        if any(strcmp(dnex(k).name(3:8),remove))
            continue
        end
    end
    if leavein
        if any(strcmp(dnex(k).name(3:8),keep))
            nexfiles{count,1}=dnex(k).name;
            count = count + 1;
        end
    else
        if any(strcmp(dnex(k).name(3:8),valid))
            nexfiles{count,1}=dnex(k).name;
            count = count + 1;
        end
    end
end
%------------------------------------------------------------------------

lastnexchar = 13;%15,17 for loading the .mat file

fs = 1000;
%filter cutoffs and order........
wcotheta    = [2 100]/(fs/2);
wcogamma    = [30 150]/(fs/2);
wcohigamma  = [100 250]/(fs/2);
nth = 4;ng = 4;nhig = 4;
%................................
% N Killian
% compute and save Granger Causality values,
% uses methods coded by  
el_pos = [0.5:.15:2.15]*1e-3;%mm->m
cond = 0.3;%S/m
cond_top = 0.3;
diam = 1*1e-3;%mm->m
gauss_sigma = 0.1*1e-3;%mm->m
filter_range = 5*gauss_sigma;

tset = [-.275 0.425];

% sheet = 'LayersMix2';
sheet = 'manual101103';
% load('D:\My Dropbox\SDVPC1data_LFP_all');
% [Num Txt RawLayers] = xlsread('D:\Dropbox\DB\Layers.xls','PCA_ITC_3grp');
[Num Txt RawLayers] = xlsread('D:\Dropbox\DB\Layers.xls',sheet);

fset=1:length(nexfiles);
% fset = 1;
nexfiles(fset)
for filenum=fset
    clear f1 f2 f3
    fid = nexfiles{filenum}(1:lastnexchar);
    if isempty(find(strcmp(fid(3:8),Raw(:,1)))),disp(['update arrayinfo.xls for ' fid]),continue;end;
    disp(fid)
    
    getlayers_all;%returns L cell mat, needs Raw from xlsread
    
    %         lay = layers(2,find(layers(1,:) == chan));
    %         otherchansinlayer = layers(1,find(layers(2,:)==lay));
    %         chansnotinlayer = layers(1,find(layers(2,:)~=lay));
    
    
    %get bad channels from arrayinfo spreadsheet--------------
    %     arraynumber = Raw{find(strcmp(fid(3:8),Raw(:,2))),3};
    bad = str2num(Raw{find(strcmp(fid(3:8),Raw(:,1))),4});
    %             good = setdiff(1:12,[bad 13]);
    super = [L{1} L{2} L{3}];
    super = setdiff(super,[bad 13]);
    %     deep = L{2};deep = setdiff(deep,[bad 13]);
    %     if isempty(super)|isempty(deep), continue;end
    if isempty(super), continue;end
    good = [super];
    old = [good];
    %---------------------------------------------------------
    
    
    
    f1   =   load(strcat(matdir,fid,tag1,'.mat'));
    f2   =   load(strcat(matdir,fid,tag2,'.mat'));
    
    for trltype = 1
        
        %         data  = eval(sprintf('f%g.data1',trltype));%load the data that's chunked into trials
        %         data  = eval(sprintf('f%g.timelock',trltype));
        data1 = f1.timelock;
        data2 = f2.timelock;
        %         time  = eval(sprintf('f%g.timelock.time',trltype));
        time = f1.timelock.time;
        labels = data1.label;
        
        %fix bad trlmats
        oldtrlmat = f1.cfg1.trl;
        f1.cfg1.trl = f1.cfg1.trlold;
        [dum oldind] = sort(f1.cfg1.trl(:,4));
        f1.cfg1.trl(:,1:4) = f1.cfg1.trl(oldind,1:4);
        %put in chronological order
        [dum Index] = sort(f1.cfg1.trl(:,2),'ascend');
        f1.cfg1.trl = f1.cfg1.trl(Index,:);
        f1.cfg1 = cliptrials(f1.cfg1);
        numcorrect = length(find(oldtrlmat(:,4) == f1.cfg1.trl(:,4)));
        nt = size(oldtrlmat,1);
        if nt~=numcorrect, error('checktrials');end
        
        oldtrlmat = f2.cfg1.trl;
        f2.cfg1.trl = f2.cfg1.trlold;
        [dum oldind] = sort(f2.cfg1.trl(:,4));
        f2.cfg1.trl(:,1:4) = f2.cfg1.trl(oldind,1:4);
        %put in chronological order
        [dum Index] = sort(f2.cfg1.trl(:,2),'ascend');
        f2.cfg1.trl = f2.cfg1.trl(Index,:);
        f2.cfg1 = cliptrials(f2.cfg1);
        numcorrect = length(find(oldtrlmat(:,4) == f2.cfg1.trl(:,4)));
        nt = size(oldtrlmat,1);
        if nt~=numcorrect, error('checktrials');end
        
        
        trlmat = f1.cfg1.trl;
        trlsel = find(trlmat(:,5)>=tset(2)*fs);
        data1.trial = data1.trial(trlsel,:,:);
        trlmat1 = trlmat(trlsel,:);
        trls{1} = 1:size(trlmat1,1);
        
        trlmat = f2.cfg1.trl;
        trlsel = find(trlmat(:,6)>=tset(2)*fs);
        data2.trial = data2.trial(trlsel,:,:);
        trlmat2 = trlmat(trlsel,:);
        trls{2} = 1:size(trlmat2,1);
        Nr1 = size(trlmat1,1);Nr2 = size(trlmat2,1);
        Nrm = min([Nr1 Nr2]);
        Nr1
        Nr2
        %         [dum sortedtrls] = sort(trlmat(:,7),'descend');%sort by NP
        
        %         hitrls = sortedtrls(1:30);
        %         lotrls = sortedtrls(end-29:end);
        %         Nr = 30;
        
        %     if hionly
        %         trls = hitrls;
        %     else
        %         trls = 1:size(trlmat,1);
        %     end
        %     Nr = length(trls);
        
        % to accomodate differing dx's, should be dxleft*dxright
        %important for sorting channels:
        locstmp        = [500 650 800 950 1100 1250 1400 1550 1700 1850 2000 2150 0]/1000;
        
        timei = [nearest(time,tset(1)):nearest(time,tset(2))];
        
        
        
        %         good = [1:3 10:12];
        %         super = [1:3];
        %         deep = [10:12];
        
        %         good = [3 12];
        
        %         data.avg = zeros(size(data.trial{1},1),size(data.trial{1},2));
        %         for ch = 1:size(data.trial{1},1)
        %             alltrls = zeros(length(data.trial),size(data.trial{1},2));
        %             for kd = 1:length(data.trial)
        %                 alltrls(kd,:) = data.trial{kd}(ch,:);
        %             end
        %             data.avg(ch,:) = nanmean(alltrls,1);
        %         end
        
        rows = [];chi = 1;superinds = [];
        %         deepinds = [];
        for Cvs = good
            Aname = ['AD' leadz(Cvs,2)];
            ACv = Cvs+100;
            row = findpowrow(labels,Aname);
            rows(chi) = row;
            if ismember(Cvs,super), superinds = [superinds chi];end
            %             if ismember(Cvs,deep), deepinds = [deepinds chi];end
            chi = chi+1;
        end
        
        el_pos = locstmp(good);
        [el_pos si] = sort(el_pos);
        rows = rows(si);
        
        switch dataformat
            %             case 'rawdata'
            %                 nvar = length(rows);
            %                 %                 Nr = length(data.trial);
            %                 %                 Nl = size(data.trial{1},2);
            %                 Nl = length(timei);
            %                 data.cat = zeros(nvar,Nr*Nl);
            %                 %concatenate trials
            %                 for chi = 1:nvar
            %                     for ri = 1:Nr
            %                         trli = trls(ri);
            %                         data.cat(chi,Nl*(ri-1)+1:Nl*ri) = data.trial{trli}(rows(chi),timei);
            %                     end
            %                 end
            case 'timelock'
                nvar = length(rows);
                %                 Nr = size(data.trial,1);
                %                 Nl = size(data.trial,3);
                Nl = length(timei);
                %                 data.cat = zeros(nvar,Nr*Nl);
                
                data_all = [];
                %concatenate trials
                for cond = 1:2
                    for chi = 1:nvar
                        for ri = 1:length(trls{cond})
                            trli = trls{cond}(ri);
                            %                             data.cat(chi,Nl*(ri-1)+1:Nl*ri) = data.trial(trli,rows(chi),timei);
                            if cond == 1
                                data_all{cond}(:,chi,ri) = squeeze(data1.trial(trli,rows(chi),timei));
                            elseif cond == 2
                                data_all{cond}(:,chi,ri) = squeeze(data2.trial(trli,rows(chi),timei));
                            end
                        end
                    end
                end
        end
        
        %         X = data.avg(rows(si),1:10000);
        
        %         X = data.cat;
        %         good
%         superinds;
        disp('doing causality analysis')
        if length(super)>1
            Cs = nchoosek(superinds,2);
            for superset = 1:size(Cs,1)
                ch1 = Cs(superset,1);
                ch2 = Cs(superset,2);
                if ch1 == ch2, continue;end
                
                fname=[newmatdir fid  '_CSDencrec_' num2str(good(ch1)) '_' num2str(good(ch2)) '.mat'];
                disp(fname)
                
                
                winsize = 150;
                windowstep = 10;
                maxsteps = floor((diff(tset)*1000-winsize)/windowstep);
                
                stage=maxsteps; % total window steps
                for ia=1:stage %do entire analysis separately for each window and condition
                    for ja=1:2 %condition
                        if ja == 1, Nr = Nr1;elseif ja == 2, Nr = Nr2;end
                        
                        LFPsignalA = squeeze(data_all{ja}(:,ch1,:));%samples x trials
                        LFPsignalB = squeeze(data_all{ja}(:,ch2,:));
                        %150 datapts, Data is samples x channel x trial
                        
                        Data=zeros(winsize,2,Nr);
                        for ka=1:Nr
                            Data(:,1,ka)=LFPsignalA(1+(ia-1)*windowstep:1+(ia-1)*windowstep+149,ka);%sig 1, up/top
                            Data(:,2,ka)=LFPsignalB(1+(ia-1)*windowstep:1+(ia-1)*windowstep+149,ka);%sig 2, down/bottom
                        end
                        
                        gagca;
                    end
                end
                
                
                
                clear OrderT gbottomupsel gtopdownsel gseltmp hsel pcsel R0hat nspec pcsellAll R0hatAll v logres cic pchat psel fsic ARloop Data
                clear data;
                clear LFPsignalA;
                clear LFPsignalB;
                clear TypeAStartin TypeAStartout TypeAStart Startin Startout;
                clear TypeAEnd;
                
                save(fname, 'TypeACoh', 'TypeAGrangerDown', 'TypeAGrangerUp', 'TypeAOrder','tset','super','good')
                
                clear TypeACoh TypeAGrangerDown TypeAGrangerUp TypeAOrder tmp LFP1signal LFP2signal AD15 AD16 LFP1filename LFP2filename LFP1fname LFP2fname
                
                
                disp('finished analyzing pair')
                
            end
            
            %             save([newmatdir fid tag1(1:end-4) '_iCSD_trltyp' num2str(trltype)],'CSD_cs','zs','time','el_pos','good')
        end
        
        %         if length(deep)>1
        %             Cd = nchoosek(deepinds,2);
        %             for deepset = 1:size(Cd,1)
        %                 ch1 = Cd(deepset,1);
        %                 ch2 = Cd(deepset,2);
        %
        %                 if ch1 == ch2, continue;end
        %
        %                 fname=[newmatdir fid  '_' num2str(good(ch1)) '_' num2str(good(ch2)) '.mat'];
        %                 disp(fname)
        %
        %
        %                 winsize = 150;
        %                 windowstep = 10;
        %                 maxsteps = floor((diff(tset)*1000-winsize)/windowstep);
        %
        %                 stage=maxsteps; % total window steps
        %                 for ia=1:stage %do entire analysis separately for each window and attention condition
        %                     for ja=1:2 %condition
        %
        %                         LFPsignalA = squeeze(data_all{ja}(:,ch1,:));%samples x trials
        %                         LFPsignalB = squeeze(data_all{ja}(:,ch2,:));
        %                         %150 datapts, Data is samples x channel x trial
        %
        %                         Data=zeros(winsize,2,Nr);
        %                         for ka=1:Nr
        %                             Data(:,1,ka)=LFPsignalA(1+(ia-1)*windowstep:1+(ia-1)*windowstep+149,ka);%sig 1, up/top
        %                             Data(:,2,ka)=LFPsignalB(1+(ia-1)*windowstep:1+(ia-1)*windowstep+149,ka);%sig 2, down/bottom
        %                         end
        %
        %                         gagca;
        %                     end
        %                 end
        %
        %
        %
        %                 clear OrderT gbottomupsel gtopdownsel gseltmp hsel pcsel R0hat nspec pcsellAll R0hatAll v logres cic pchat psel fsic ARloop Data
        %                 clear data;
        %                 clear LFPsignalA;
        %                 clear LFPsignalB;
        %                 clear TypeAStartin TypeAStartout TypeAStart Startin Startout;
        %                 clear TypeAEnd;
        %
        %                 save(fname, 'TypeACoh', 'TypeAGrangerDown', 'TypeAGrangerUp', 'TypeAOrder','tset','super','deep','good')
        %
        %                 clear TypeACoh TypeAGrangerDown TypeAGrangerUp TypeAOrder tmp LFP1signal LFP2signal AD15 AD16 LFP1filename LFP2filename LFP1fname LFP2fname
        %
        %
        %                 disp('finished analyzing pair')
        %
        %             end
        %
        %             %             save([newmatdir fid tag1(1:end-4) '_iCSD_trltyp' num2str(trltype)],'CSD_cs','zs','time','el_pos','good')
        %
        %         end
        
    end%end trltype
end
% end
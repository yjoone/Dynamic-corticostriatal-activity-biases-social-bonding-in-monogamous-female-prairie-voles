
% clar
for remedges = [0]
    %tN = 630;filetag = '';tset = [-0.5 0.75];bls = [-275 -75];
    tN = 350;filetag = '_ER';tset = [-.275 .425];bls = [-200 -75];
    
    figdir = 'D:\My Dropbox\ECG2010\SD_Coh_GC\';
    stat = 1;p=0.05;
    alpha = 0.05;numcont = 5;
    maxfreq = 150;
    consbin=10;
    slidingwindowsize=150;
    step=10;
    winsize = 150;
    windowstep = 10;
    maxsteps = floor((diff(tset)*1000-winsize)/windowstep);
    startpoint = -(tset(1)*1000 - step);
    
    time=1:maxsteps;time=time*step-startpoint+slidingwindowsize/2;
    %         tind = nearest(time,-300):nearest(time,tN);t = time(tind);time2 = t;
    tind = nearest(time,-200):nearest(time,tN);t = time(tind);time2 = t;baselineset = nearest(time2,-200):nearest(time2,-75);
    bli = nearest(time,bls(1)):nearest(time,bls(2));%time bins
    baselineset = nearest(time2,bls(1)):nearest(time2,bls(2));
    figure(201);clf(201);figure(201);
    %% --------------------------------------------------------
    sheet = 'manual101103';
    %     laysets = [1 2;1 3;2 3];
    %     laysets = [1 3];
    %     laysets = [1 1;2 2;3 3];
    laysets = [1 2];
    for layi = 1:size(laysets,1)
        
        lay1 = laysets(layi,1);lay2 = laysets(layi,2);
        FT = ['D:\' '_lay' num2str(lay1) num2str(lay2) '_' sheet '_remedges' num2str(remedges) filetag];
        FT2 = ['_p' num2str(alpha) '_' num2str(numcont) 'cont' '_' num2str(tN) 'ms_lay' num2str(lay1) num2str(lay2) '_' sheet '_remedges' num2str(remedges) filetag];
        load(FT);
        clear a1 a2 a3 a4 a
        
        GrangerUpInAll = real(GrangerUpInAll);
        GrangerUpOutAll = real(GrangerUpOutAll);
        GrangerDownInAll = real(GrangerDownInAll);
        GrangerDownOutAll = real(GrangerDownOutAll);
        
        figure(101);clf(101);figure(101);
        
        CohInAll=abs(CohInAll);
        CohOutAll=abs(CohOutAll);
        [steps  freq n]=size(CohInAll);
        CohInMean=mean(CohInAll,3);%avg across pairs
        CohOutMean=mean(CohOutAll,3);%avg across pairs
        
        CohAtt=CohInMean-CohOutMean;
        
        
        freq=1:maxfreq;
        %     subplot(3,1,1)
        tmp=CohAtt(:,freq)';
        [freqs steps]=size(tmp);
        
        %%
        
        figure(201)
        
        
        fsets = [30 90;100 140;];
        norm2bl = 1;samebl = 1;
        for fi = 1:size(fsets,1)
            
            subplot(3,2,sub2ind([2 3],fi,layi))
            line([0 0],[-15 15],'color','k');hold on;
            line([-200 600],[1 1],'color','k')
            fset2 = fsets(fi,1):fsets(fi,2);
            CohInAll=abs(CohInAll);
            CohOutAll=abs(CohOutAll);
            if norm2bl
                if samebl
                    bl2 = repmat(nanmean(cat(1,CohInAll(bli,fset2,:),CohOutAll(bli,fset2,:)),1),[length(tind) 1 1]);
                    blin = bl2;blout = bl2;
                else
                    blin = repmat(nanmean(CohInAll(bli,fset2,:),1),[length(tind) 1 1]);
                    blout = repmat(nanmean(CohOutAll(bli,fset2,:),1),[length(tind) 1 1]);
                end
                CohInMean2=squeeze(mean(CohInAll(tind,fset2,:)./blin,2));%avg across freqs
                CohOutMean2=squeeze(mean(CohOutAll(tind,fset2,:)./blout,2));%avg across freqs
            else
                CohInMean2=squeeze(mean(CohInAll(tind,fset2,:),2));%avg across freqs
                CohOutMean2=squeeze(mean(CohOutAll(tind,fset2,:),2));%avg across freqs
            end
            
            plot(t,mean(CohInMean2,2),'r');hold on
            dofill(t,CohInMean2','r',1,1)
            plot(t,mean(CohOutMean2,2),'b');hold on
            dofill(t,CohOutMean2','b',1,1)
            
            s1 = CohInMean2;s2 = CohOutMean2;
            m1 = max((mean(s1,2)));m2 = max((mean(s2,2)));
            mval = max([m1 m2]);
            minval = min([mean(s1,2);mean(s2,2)]);
            axis([t(1) t(end) minval-.1 mval+.1])
            bw = 1;%number of bins to average for ttest computation
            bini = 1;
            hs = [];ps = [];ind = [];
            for ti = 1:bw:length(t)-bw+1
                [hs(bini) ps(bini)] = ttest(mean(s1(ti:ti+bw-1,:),1),mean(s2(ti:ti+bw-1,:),1));
                
                ind(bini) = ti;
                bini = bini + 1;
            end
            below = [0 hs];
            % bins = ind(find(hs));
            possibleonset = find(diff(below)==1);
            possibleonset(possibleonset>(length(below)-numcont+1)) = [];
            possibleonset = [possibleonset (length(hs)+1)];
            sums = [];
            for po = 1:length(possibleonset)-1
                %             sums(po) = sum(hs(possibleonset(po):possibleonset(po+1)-1));%numcontinuous after this onset
                width = find(diff(hs(possibleonset(po):possibleonset(po+1)-1)),1,'first');%numcontinuous after this onset
                if isempty(width), width = length(hs(possibleonset(po):possibleonset(po+1)-1));end
                sums(po) = width;
            end
            pos = possibleonset(sums>=numcont);
            sums = sums(sums>=numcont);
            % plot(t,smooth(nanmean(s1,2),5),'r'),hold on
            % x = t;y = s1;color = 'r';dim = 2;smval = 5;
            % dofill(x,y',color,1,smval);
            patchi = 1;
            for bin = pos
                hi = ind(bin);
                %    xl = t(hi)-0.005;xr = t(hi+bw-1)+0.004;
                dt = t(2)-t(1);
                xl = t(hi)-dt/2;
                xr = t(hi+bw*sums(patchi)-1)+dt/2;
                yb = -15;yt = 15;
                xvert = [xl xl xr xr xl];yvert = [yb yt yt yb yb];
                fill(xvert,yvert,'k', 'FaceAlpha',.3, 'EdgeColor','k','EdgeAlpha',0);
                patchi = patchi+1;
            end
            line([0 0],[-15 15],'color','k')
            line([-200 600],[1 1],'color','k')
            xlabel('time from stimulus onset (ms)')
            title(['avg LFP coh.,cross-layers, ' num2str(fset2(1)) '-' num2str(fset2(end)) ' Hz, n = ' num2str(size(s1,2))])
        end
        
        
        %% 2D PLOTS
        figure(105);clf(105);figure(105);
        figure(106);clf(106);figure(106);
        figure(501);clf(501);figure(501);
        figure(502);clf(502);figure(502);
        figure(601);clf(601);figure(601);
        fsets = [30 90;100 140];
        norm2bl = 1;samebl = 1;
        for fi = 1:size(fsets,1)
            freqset = fsets(fi,1):fsets(fi,2);
            
            %     bli = nearest(time,-275):nearest(time,-75);%time bins
            %----------------------------------
            %mean across all freqs for all pairs
            GrangerUpOutgamma=mean(GrangerUpOutAll(tind,freqset,:),2);
            GrangerUpIngamma=mean(GrangerUpInAll(tind,freqset,:),2);
            % sterrgrangerupoutgamma=std(squeeze(GrangerUpOutgamma(:,:)),0,2)/sqrt(k);
            % sterrgrangerupingamma=std(squeeze(GrangerUpIngamma(:,:)),0,2)/sqrt(k);
            
            for i=1:size(GrangerUpOutgamma,3)
                GrangerUpOutgammanormal(:,i)=squeeze(GrangerUpOutgamma(:,1,i))/mean([mean(squeeze(GrangerUpOutgamma(baselineset,1,i)),1) mean(squeeze(GrangerUpIngamma(baselineset,1,i)),1)]);
                GrangerUpIngammanormal(:,i)=squeeze(GrangerUpIngamma(:,1,i))/mean([mean(squeeze(GrangerUpOutgamma(baselineset,1,i)),1) mean(squeeze(GrangerUpIngamma(baselineset,1,i)),1)]);
            end
            % sterrgrangerupoutgammanormal=std(GrangerUpOutgammanormal,0,2)/sqrt(k);
            % sterrgrangerupingammanormal=std(GrangerUpIngammanormal,0,2)/sqrt(k);
            GrangerUp=(mean(squeeze(GrangerUpIngamma),2)-mean(squeeze(GrangerUpOutgamma),2));
            GrangerUpnormal=(mean(GrangerUpIngammanormal,2)-mean(GrangerUpOutgammanormal,2));
            
            %     figure(104);clf(104);figure(104)%---------------------------------------------------------
            %     % subplot(2,1,2)
            %     plot(time,mean(squeeze(GrangerUpIngamma),2),'r')
            %     hold on
            %     dofill(time,squeeze(GrangerUpIngamma)','r',1,1)
            %     plot(time, mean(squeeze(GrangerUpOutgamma),2),'b')
            %     dofill(time,squeeze(GrangerUpOutgamma)','b',1,1)
            %
            %     % fill([time fliplr(time)] ,[(mean(squeeze(GrangerUpOutgamma),2)+sterrgrangerupoutgamma)' fliplr((mean(squeeze(GrangerUpOutgamma),2)-sterrgrangerupoutgamma)')],'b', 'FaceAlpha',.15, 'EdgeColor','b','EdgeAlpha',.1);
            %     xlabel('time (ms)'),ylabel('normalized GC')
            %     title(['b/r ' num2str(lay2) '->' num2str(lay1) ', m/g ' num2str(lay1) '->' num2str(lay2) ''])
            
            GrangerDownOutgamma=mean(GrangerDownOutAll(tind,freqset,:),2);
            GrangerDownIngamma=mean(GrangerDownInAll(tind,freqset,:),2);
            % sterrgrangerdowningamma=std(squeeze(GrangerDownIngamma(:,:)),0,2)/sqrt(k);
            % sterrgrangerdownoutgamma=std(squeeze(GrangerDownOutgamma(:,:)),0,2)/sqrt(k);
            
            for i=1:size(GrangerDownOutgamma,3)
                GrangerDownOutgammanormal(:,i)=squeeze(GrangerDownOutgamma(:,1,i))/mean([mean(squeeze(GrangerDownOutgamma(baselineset,1,i)),1) mean(squeeze(GrangerDownIngamma(baselineset,1,i)),1)]);
                GrangerDownIngammanormal(:,i)=squeeze(GrangerDownIngamma(:,1,i))/mean([mean(squeeze(GrangerDownOutgamma(baselineset,1,i)),1) mean(squeeze(GrangerDownIngamma(baselineset,1,i)),1)]);
            end
            % sterrgrangerdownoutgammanormal=std(GrangerDownOutgammanormal,0,2)/sqrt(k);
            % sterrgrangerdowningammanormal=std(GrangerDownIngammanormal,0,2)/sqrt(k);
            
            
            GrangerDown=(mean(squeeze(GrangerDownIngamma),2)-mean(squeeze(GrangerDownOutgamma),2));
            
            GrangerDownnormal=(mean(GrangerDownIngammanormal,2)-mean(GrangerDownOutgammanormal,2));
            
            %     % subplot(2,1,1)
            %     plot(time,mean(squeeze(GrangerDownIngamma),2),'m')
            %     hold on
            %     dofill(time,squeeze(GrangerDownIngamma)','m',1,1)
            %     plot(time, mean(squeeze(GrangerDownOutgamma),2),'g')
            %     dofill(time,squeeze(GrangerDownOutgamma)','g',1,1)
            %     title(['GC ' num2str(lay1) '->' num2str(lay2) ', hi:red, lo:blue, unnormed'])
            %     xlabel('time (ms)'),ylabel('normalized GC')
            
            % axis([-300 500 0.015 .055])
            figure(105)
            %             pn = sub2ind([2 2],fi,2);
            subplot(2,2,sub2ind([2 2],fi,2))
            
            plot(time2,mean(abs(GrangerDownIngammanormal),2),'r')
            hold on
            dofill(time2,GrangerDownIngammanormal','r',1,1)
            plot(time2, mean(GrangerDownOutgammanormal,2),'b')
            dofill(time2,GrangerDownOutgammanormal','b',1,1)
            dopptt(abs(GrangerDownOutgammanormal),abs(GrangerDownIngammanormal),time2,1,5,0)
            xlabel('time (ms)'),ylabel('normalized GC')
            title(['layer a -> layer b,' num2str(fsets(fi,1)) '-' num2str(fsets(fi,2)) ' Hz, n = ' num2str(size(GrangerDownIngammanormal,2))])
            %             legend('strong encoding', 'weak encoding','location','best')
            line([-200 650],[1 1],'color','k');hold on;            line([0 0],[-15 15],'color','k')
            %             axis([-100 600 0.8 2])
            %             subplot(3,2,sub2ind([2 3],fi,layi))
            m1 = max(mean(abs(GrangerDownIngammanormal),2));m2 = max(mean(abs(GrangerDownOutgammanormal),2));
            mval = max([m1 m2]);
            minval = min([mean(abs(GrangerDownIngammanormal),2)' mean(abs(GrangerDownOutgammanormal),2)']);
            axis([-100 tN minval-.1 mval+.1])
            
            
            subplot(2,2,sub2ind([2 2],fi,1))
            plot(time2,mean(GrangerUpIngammanormal,2),'r')
            hold on
            dofill(time2,GrangerUpIngammanormal','r',1,1)
            plot(time2, mean(GrangerUpOutgammanormal,2),'b')
            dofill(time2,GrangerUpOutgammanormal','b',1,1)
            dopptt(abs(GrangerUpOutgammanormal),abs(GrangerUpIngammanormal),time2,1,5,0)
            title(['layer b -> layer a,' num2str(fsets(fi,1)) '-' num2str(fsets(fi,2)) ' Hz, n = ' num2str(size(GrangerUpIngammanormal,2))])
            %             title(['GC ' num2str(lay2) '->' num2str(lay1) ', hi:red, lo:blue, theta, Normalized to baseline']);
            xlabel('time (ms)'),ylabel('normalized GC')
            %             legend('strong encoding', 'weak encoding','location','best')
            line([-200 650],[1 1],'color','k');hold on;            line([0 0],[-15 15],'color','k')
            axis([-100 tN 0.8 2])
            
            m1 = max(mean(abs(GrangerUpIngammanormal),2));m2 = max(mean(abs(GrangerUpOutgammanormal),2));
            mval = max([m1 m2]);
            minval = min([mean(abs(GrangerUpIngammanormal),2)' mean(abs(GrangerUpOutgammanormal),2)']);
            axis([-100 tN minval-.1 mval+.1])
            
            
            
            %put layers on same plot
            figure(501);
            %             pn = sub2ind([2 2],fi,1);
            subplot(2,2,sub2ind([2 2],fi,1))
            
            plot(time2,mean(abs(GrangerDownIngammanormal),2),'r')
            hold on
            dofill(time2,GrangerDownIngammanormal','r',1,1)
            
            
            
            plot(time2, mean(GrangerUpIngammanormal,2),'b')
            dofill(time2,GrangerUpIngammanormal','b',1,1)
            
            dopptt(abs(GrangerDownIngammanormal),abs(GrangerUpIngammanormal),time2,1,5,0)
            
            
            xlabel('time (ms)'),ylabel('normalized GC')
            title(['strong encoding' num2str(fsets(fi,1)) '-' num2str(fsets(fi,2)) ' Hz, n = ' num2str(size(GrangerDownIngammanormal,2))])
            %             legend('a to b', 'b to a','location','best')
            line([-200 650],[1 1],'color','k');hold on;            line([0 0],[-15 15],'color','k')
            %             axis([-100 600 0.8 2])
            %             subplot(3,2,sub2ind([2 3],fi,layi))
            m1 = max(mean(abs(GrangerDownIngammanormal),2));m2 = max(mean(abs(GrangerUpIngammanormal),2));
            mval = max([m1 m2]);
            minval = min([mean(abs(GrangerDownIngammanormal),2)' mean(abs(GrangerUpIngammanormal),2)']);
            axis([-100 tN minval-.1 mval+.1])
            
            
            subplot(2,2,sub2ind([2 2],fi,2))
            plot(time2,mean(GrangerDownOutgammanormal,2),'r')
            hold on
            dofill(time2,GrangerDownOutgammanormal','r',1,1)
            
            plot(time2, mean(GrangerUpOutgammanormal,2),'b')
            dofill(time2,GrangerUpOutgammanormal','b',1,1)
            
            dopptt(abs(GrangerUpOutgammanormal),abs(GrangerDownOutgammanormal),time2,1,5,0)
            
            
            title(['weak encoding,' num2str(fsets(fi,1)) '-' num2str(fsets(fi,2)) ' Hz, n = ' num2str(size(GrangerUpIngammanormal,2))])
            %             title(['GC ' num2str(lay2) '->' num2str(lay1) ', hi:red, lo:blue, theta, Normalized to baseline']);
            xlabel('time (ms)'),ylabel('normalized GC')
            %             legend('a to b', 'b to a','location','best')
            line([-200 650],[1 1],'color','k');hold on;            line([0 0],[-15 15],'color','k')
            axis([-100 tN 0.8 2])
            
            
            m1 = max(mean(abs(GrangerDownOutgammanormal),2));m2 = max(mean(abs(GrangerUpOutgammanormal),2));
            mval = max([m1 m2]);
            minval = min([mean(abs(GrangerDownOutgammanormal),2)' mean(abs(GrangerUpOutgammanormal),2)']);
            axis([-100 tN minval-.1 mval+.1])
            
            %             figure(502)
            
            %             m1 = max(mean(abs(GrangerDownIngammanormal,2)));m2 = max(mean(abs(GrangerDownOutgammanormal,2)));
            %             mval = max([m1 m2]);
            %             minval = min([mean(abs(GrangerDownIngammanormal,2)) mean(abs(GrangerDownOutgammanormal,2))]);
            %             axis([t(1) t(end) minval-.1 mval+.1])
            % title(['GC ' num2str(lay1) '->' num2str(lay2) ', hi:red, lo:blue, Normalized to baseline');
            %             title(['b/r ' num2str(lay2) '->' num2str(lay1) ', m/g ' num2str(lay1) '->' num2str(lay2) ''])
            figure(601)
            %              GrangerUpOutgamma=mean(GrangerUpOutAll(tind,freqset,:),2);
            %             GrangerUpIngamma=mean(GrangerUpInAll(tind,freqset,:),2);
            subplot(1,2,sub2ind([2 1],fi,1))
            GrangerUp = mean(cat(2,GrangerUpOutAll(tind,freqset,:),GrangerUpInAll(tind,freqset,:)),2);
            
            for i=1:size(GrangerUpOutgamma,3)
                %                 GrangerUpOutgammanormal(:,i)=squeeze(GrangerUpOutgamma(:,1,i))/mean([mean(squeeze(GrangerUpOutgamma(baselineset,1,i)),1) mean(squeeze(GrangerUpIngamma(baselineset,1,i)),1)]);
                %                 GrangerUpIngammanormal(:,i)=squeeze(GrangerUpIngamma(:,1,i))/mean([mean(squeeze(GrangerUpOutgamma(baselineset,1,i)),1) mean(squeeze(GrangerUpIngamma(baselineset,1,i)),1)]);
                GrangerUpnormal(:,i)=squeeze(GrangerUp(:,1,i))/mean([mean(squeeze(GrangerUp(baselineset,1,i)),1)]);
            end
            
            %                          GrangerDownOutgamma=mean(GrangerDownOutAll(tind,freqset,:),2);
            %             GrangerDownIngamma=mean(GrangerDownInAll(tind,freqset,:),2);
            GrangerDown = mean(cat(2,GrangerDownOutAll(tind,freqset,:),GrangerDownInAll(tind,freqset,:)),2);
            
            for i=1:size(GrangerDownOutgamma,3)
                %                 GrangerDownOutgammanormal(:,i)=squeeze(GrangerDownOutgamma(:,1,i))/mean([mean(squeeze(GrangerDownOutgamma(baselineset,1,i)),1) mean(squeeze(GrangerDownIngamma(baselineset,1,i)),1)]);
                %                 GrangerDownIngammanormal(:,i)=squeeze(GrangerDownIngamma(:,1,i))/mean([mean(squeeze(GrangerDownOutgamma(baselineset,1,i)),1) mean(squeeze(GrangerDownIngamma(baselineset,1,i)),1)]);
                GrangerDownnormal(:,i)=squeeze(GrangerDown(:,1,i))/mean([mean(squeeze(GrangerDown(baselineset,1,i)),1)]);
            end
            
            plot(time2,mean(GrangerDownnormal,2),'r')
            hold on
            dofill(time2,GrangerDownnormal','r',1,1)
            plot(time2, mean(GrangerUpnormal,2),'b')
            dofill(time2,GrangerUpnormal','b',1,1)
            dopptt(abs(GrangerUpnormal),abs(GrangerDownnormal),time2,1,5,0)
            
            xlabel('time (ms)'),ylabel('normalized GC')
            title(['all trials' num2str(fsets(fi,1)) '-' num2str(fsets(fi,2)) ' Hz, n = ' num2str(size(GrangerDownnormal,2))])
            %             legend('a to b', 'b to a','location','best')
            line([-200 650],[1 1],'color','k');hold on;            line([0 0],[-15 15],'color','k')
            %             axis([-100 600 0.8 2])
            %             subplot(3,2,sub2ind([2 3],fi,layi))
            m1 = max(mean(abs(GrangerDownnormal),2));m2 = max(mean(abs(GrangerUpnormal),2));
            mval = max([m1 m2]);
            minval = min([mean(abs(GrangerDownnormal),2)' mean(abs(GrangerUpnormal),2)']);
            axis([-100 tN minval-.1 mval+.1])
            
        end
        
        
        fsets = [30 140];
        
        for fi = 1:size(fsets,1)
            freqset = fsets(fi,1):fsets(fi,2);
            
            %     bli = nearest(time,-275):nearest(time,-75);%time bins
            %----------------------------------
            %mean across all freqs for all pairs
            GrangerUpOutgamma=mean(GrangerUpOutAll(tind,freqset,:),2);
            GrangerUpIngamma=mean(GrangerUpInAll(tind,freqset,:),2);
            % sterrgrangerupoutgamma=std(squeeze(GrangerUpOutgamma(:,:)),0,2)/sqrt(k);
            % sterrgrangerupingamma=std(squeeze(GrangerUpIngamma(:,:)),0,2)/sqrt(k);
            
            for i=1:size(GrangerUpOutgamma,3)
                GrangerUpOutgammanormal(:,i)=squeeze(GrangerUpOutgamma(:,1,i))/mean([mean(squeeze(GrangerUpOutgamma(baselineset,1,i)),1) mean(squeeze(GrangerUpIngamma(baselineset,1,i)),1)]);
                GrangerUpIngammanormal(:,i)=squeeze(GrangerUpIngamma(:,1,i))/mean([mean(squeeze(GrangerUpOutgamma(baselineset,1,i)),1) mean(squeeze(GrangerUpIngamma(baselineset,1,i)),1)]);
            end
            % sterrgrangerupoutgammanormal=std(GrangerUpOutgammanormal,0,2)/sqrt(k);
            % sterrgrangerupingammanormal=std(GrangerUpIngammanormal,0,2)/sqrt(k);
            GrangerUp=(mean(squeeze(GrangerUpIngamma),2)-mean(squeeze(GrangerUpOutgamma),2));
            GrangerUpnormal=(mean(GrangerUpIngammanormal,2)-mean(GrangerUpOutgammanormal,2));
            
            %     figure(104);clf(104);figure(104)%---------------------------------------------------------
            %     % subplot(2,1,2)
            %     plot(time,mean(squeeze(GrangerUpIngamma),2),'r')
            %     hold on
            %     dofill(time,squeeze(GrangerUpIngamma)','r',1,1)
            %     plot(time, mean(squeeze(GrangerUpOutgamma),2),'b')
            %     dofill(time,squeeze(GrangerUpOutgamma)','b',1,1)
            %
            %     % fill([time fliplr(time)] ,[(mean(squeeze(GrangerUpOutgamma),2)+sterrgrangerupoutgamma)' fliplr((mean(squeeze(GrangerUpOutgamma),2)-sterrgrangerupoutgamma)')],'b', 'FaceAlpha',.15, 'EdgeColor','b','EdgeAlpha',.1);
            %     xlabel('time (ms)'),ylabel('normalized GC')
            %     title(['b/r ' num2str(lay2) '->' num2str(lay1) ', m/g ' num2str(lay1) '->' num2str(lay2) ''])
            
            GrangerDownOutgamma=mean(GrangerDownOutAll(tind,freqset,:),2);
            GrangerDownIngamma=mean(GrangerDownInAll(tind,freqset,:),2);
            % sterrgrangerdowningamma=std(squeeze(GrangerDownIngamma(:,:)),0,2)/sqrt(k);
            % sterrgrangerdownoutgamma=std(squeeze(GrangerDownOutgamma(:,:)),0,2)/sqrt(k);
            
            for i=1:size(GrangerDownOutgamma,3)
                GrangerDownOutgammanormal(:,i)=squeeze(GrangerDownOutgamma(:,1,i))/mean([mean(squeeze(GrangerDownOutgamma(baselineset,1,i)),1) mean(squeeze(GrangerDownIngamma(baselineset,1,i)),1)]);
                GrangerDownIngammanormal(:,i)=squeeze(GrangerDownIngamma(:,1,i))/mean([mean(squeeze(GrangerDownOutgamma(baselineset,1,i)),1) mean(squeeze(GrangerDownIngamma(baselineset,1,i)),1)]);
            end
            % sterrgrangerdownoutgammanormal=std(GrangerDownOutgammanormal,0,2)/sqrt(k);
            % sterrgrangerdowningammanormal=std(GrangerDownIngammanormal,0,2)/sqrt(k);
            
            
            GrangerDown=(mean(squeeze(GrangerDownIngamma),2)-mean(squeeze(GrangerDownOutgamma),2));
            
            GrangerDownnormal=(mean(GrangerDownIngammanormal,2)-mean(GrangerDownOutgammanormal,2));
            
            %     % subplot(2,1,1)
            %     plot(time,mean(squeeze(GrangerDownIngamma),2),'m')
            %     hold on
            %     dofill(time,squeeze(GrangerDownIngamma)','m',1,1)
            %     plot(time, mean(squeeze(GrangerDownOutgamma),2),'g')
            %     dofill(time,squeeze(GrangerDownOutgamma)','g',1,1)
            %     title(['GC ' num2str(lay1) '->' num2str(lay2) ', hi:red, lo:blue, unnormed'])
            %     xlabel('time (ms)'),ylabel('normalized GC')
            
            % axis([-300 500 0.015 .055])
            figure(106)
            %             pn = sub2ind([2 2],fi,2);
            subplot(1,2,sub2ind([2 1],1,1))
            
            plot(time2,mean(abs(GrangerDownIngammanormal),2),'r')
            hold on
            dofill(time2,GrangerDownIngammanormal','r',1,1)
            plot(time2, mean(GrangerDownOutgammanormal,2),'b')
            dofill(time2,GrangerDownOutgammanormal','b',1,1)
            dopptt(abs(GrangerDownOutgammanormal),abs(GrangerDownIngammanormal),time2,1,5,0)
            xlabel('time (ms)'),ylabel('normalized GC')
            title(['layer a -> layer b,' num2str(fsets(fi,1)) '-' num2str(fsets(fi,2)) ' Hz, n = ' num2str(size(GrangerDownIngammanormal,2))])
            %             legend('strong encoding', 'weak encoding','location','best')
            line([-200 650],[1 1],'color','k');hold on;            line([0 0],[-15 15],'color','k')
            %             axis([-100 600 0.8 2])
            %             subplot(3,2,sub2ind([2 3],fi,layi))
            m1 = max(mean(abs(GrangerDownIngammanormal),2));m2 = max(mean(abs(GrangerDownOutgammanormal),2));
            mval = max([m1 m2]);
            minval = min([mean(abs(GrangerDownIngammanormal),2)' mean(abs(GrangerDownOutgammanormal),2)']);
            axis([-100 tN .9 1.5])
            
            
            subplot(1,2,sub2ind([2 1],2,1))
            plot(time2,mean(GrangerUpIngammanormal,2),'r')
            hold on
            dofill(time2,GrangerUpIngammanormal','r',1,1)
            plot(time2, mean(GrangerUpOutgammanormal,2),'b')
            dofill(time2,GrangerUpOutgammanormal','b',1,1)
            dopptt(abs(GrangerUpOutgammanormal),abs(GrangerUpIngammanormal),time2,1,5,0)
            title(['layer b -> layer a,' num2str(fsets(fi,1)) '-' num2str(fsets(fi,2)) ' Hz, n = ' num2str(size(GrangerUpIngammanormal,2))])
            %             title(['GC ' num2str(lay2) '->' num2str(lay1) ', hi:red, lo:blue, theta, Normalized to baseline']);
            xlabel('time (ms)'),ylabel('normalized GC')
            %             legend('strong encoding', 'weak encoding','location','best')
            line([-200 650],[1 1],'color','k');hold on;            line([0 0],[-15 15],'color','k')
            axis([-100 tN 0.8 2])
            
            m1 = max(mean(abs(GrangerUpIngammanormal),2));m2 = max(mean(abs(GrangerUpOutgammanormal),2));
            mval = max([m1 m2]);
            minval = min([mean(abs(GrangerUpIngammanormal),2)' mean(abs(GrangerUpOutgammanormal),2)']);
            axis([-100 tN .9 1.5])
            
            
            
            %put layers on same plot
            figure(502);
            %             pn = sub2ind([2 2],fi,1);
            subplot(1,2,sub2ind([2 1],1,1))
            
            plot(time2,mean(abs(GrangerDownIngammanormal),2),'r')
            hold on
            dofill(time2,GrangerDownIngammanormal','r',1,1)
            
            
            
            plot(time2, mean(GrangerUpIngammanormal,2),'b')
            dofill(time2,GrangerUpIngammanormal','b',1,1)
            
            dopptt(abs(GrangerDownIngammanormal),abs(GrangerUpIngammanormal),time2,1,5,0)
            
            
            xlabel('time (ms)'),ylabel('normalized GC')
            title(['strong encoding' num2str(fsets(fi,1)) '-' num2str(fsets(fi,2)) ' Hz, n = ' num2str(size(GrangerDownIngammanormal,2))])
            %             legend('a to b', 'b to a','location','best')
            line([-200 650],[1 1],'color','k');hold on;            line([0 0],[-15 15],'color','k')
            %             axis([-100 600 0.8 2])
            %             subplot(3,2,sub2ind([2 3],fi,layi))
            m1 = max(mean(abs(GrangerDownIngammanormal),2));m2 = max(mean(abs(GrangerUpIngammanormal),2));
            mval = max([m1 m2]);
            minval = min([mean(abs(GrangerDownIngammanormal),2)' mean(abs(GrangerUpIngammanormal),2)']);
            axis([-100 tN minval-.1 mval+.1])
            
            
            subplot(1,2,sub2ind([2 1],2,1))
            plot(time2,mean(GrangerDownOutgammanormal,2),'r')
            hold on
            dofill(time2,GrangerDownOutgammanormal','r',1,1)
            
            plot(time2, mean(GrangerUpOutgammanormal,2),'b')
            dofill(time2,GrangerUpOutgammanormal','b',1,1)
            
            dopptt(abs(GrangerUpOutgammanormal),abs(GrangerDownOutgammanormal),time2,1,5,0)
            
            
            title(['weak encoding,' num2str(fsets(fi,1)) '-' num2str(fsets(fi,2)) ' Hz, n = ' num2str(size(GrangerUpIngammanormal,2))])
            %             title(['GC ' num2str(lay2) '->' num2str(lay1) ', hi:red, lo:blue, theta, Normalized to baseline']);
            xlabel('time (ms)'),ylabel('normalized GC')
            %             legend('a to b', 'b to a','location','best')
            line([-200 650],[1 1],'color','k');hold on;            line([0 0],[-15 15],'color','k')
            axis([-100 tN 0.8 2])
            
            
            m1 = max(mean(abs(GrangerDownOutgammanormal),2));m2 = max(mean(abs(GrangerUpOutgammanormal),2));
            mval = max([m1 m2]);
            minval = min([mean(abs(GrangerDownOutgammanormal),2)' mean(abs(GrangerUpOutgammanormal),2)']);
            axis([-100 tN minval-.1 mval+.1])
            
            %             figure(502)
            
            %             m1 = max(mean(abs(GrangerDownIngammanormal,2)));m2 = max(mean(abs(GrangerDownOutgammanormal,2)));
            %             mval = max([m1 m2]);
            %             minval = min([mean(abs(GrangerDownIngammanormal,2)) mean(abs(GrangerDownOutgammanormal,2))]);
            %             axis([t(1) t(end) minval-.1 mval+.1])
            % title(['GC ' num2str(lay1) '->' num2str(lay2) ', hi:red, lo:blue, Normalized to baseline');
            %             title(['b/r ' num2str(lay2) '->' num2str(lay1) ', m/g ' num2str(lay1) '->' num2str(lay2) ''])
            
            figure(602)
            %              GrangerUpOutgamma=mean(GrangerUpOutAll(tind,freqset,:),2);
            %             GrangerUpIngamma=mean(GrangerUpInAll(tind,freqset,:),2);
            GrangerUp = mean(cat(2,GrangerUpOutAll(tind,freqset,:),GrangerUpInAll(tind,freqset,:)),2);
            
            for i=1:size(GrangerUpOutgamma,3)
                %                 GrangerUpOutgammanormal(:,i)=squeeze(GrangerUpOutgamma(:,1,i))/mean([mean(squeeze(GrangerUpOutgamma(baselineset,1,i)),1) mean(squeeze(GrangerUpIngamma(baselineset,1,i)),1)]);
                %                 GrangerUpIngammanormal(:,i)=squeeze(GrangerUpIngamma(:,1,i))/mean([mean(squeeze(GrangerUpOutgamma(baselineset,1,i)),1) mean(squeeze(GrangerUpIngamma(baselineset,1,i)),1)]);
                GrangerUpnormal(:,i)=squeeze(GrangerUp(:,1,i))/([mean(squeeze(GrangerUp(baselineset,1,i)),1)]);
            end
            
            %                          GrangerDownOutgamma=mean(GrangerDownOutAll(tind,freqset,:),2);
            %             GrangerDownIngamma=mean(GrangerDownInAll(tind,freqset,:),2);
            GrangerDown = mean(cat(2,GrangerDownOutAll(tind,freqset,:),GrangerDownInAll(tind,freqset,:)),2);
            
            for i=1:size(GrangerDownOutgamma,3)
                %                 GrangerDownOutgammanormal(:,i)=squeeze(GrangerDownOutgamma(:,1,i))/mean([mean(squeeze(GrangerDownOutgamma(baselineset,1,i)),1) mean(squeeze(GrangerDownIngamma(baselineset,1,i)),1)]);
                %                 GrangerDownIngammanormal(:,i)=squeeze(GrangerDownIngamma(:,1,i))/mean([mean(squeeze(GrangerDownOutgamma(baselineset,1,i)),1) mean(squeeze(GrangerDownIngamma(baselineset,1,i)),1)]);
                GrangerDownnormal(:,i)=squeeze(GrangerDown(:,1,i))/([mean(squeeze(GrangerDown(baselineset,1,i)),1)]);
            end
            
            %             plot(time2,mean(GrangerDownnormal,2),'r')
            hold on
            dofill(time2,GrangerDownnormal','r',1,1)
            %             plot(time2, mean(GrangerUpnormal,2),'b')
            dofill(time2,GrangerUpnormal','b',1,1)
            dopptt(abs(GrangerUpnormal),abs(GrangerDownnormal),time2,1,5,0)
            
            xlabel('time (ms)'),ylabel('normalized GC')
            title(['all trials' num2str(fsets(fi,1)) '-' num2str(fsets(fi,2)) ' Hz, n = ' num2str(size(GrangerDownnormal,2))])
            %             legend('a to b', 'b to a','location','best')
            line([-200 650],[1 1],'color','k');hold on;            line([0 0],[-15 15],'color','k')
            %             axis([-100 600 0.8 2])
            %             subplot(3,2,sub2ind([2 3],fi,layi))
            m1 = max(mean(abs(GrangerDownnormal),2));m2 = max(mean(abs(GrangerUpnormal),2));
            mval = max([m1 m2]);
            minval = min([mean(abs(GrangerDownnormal),2)' mean(abs(GrangerUpnormal),2)']);
            axis off,box off
            axis([-100 tN .95 1.5])
            
            
            
        end
        %%
        
        combineconds = 0;
        norm2bl = 0;
        figure(701);clf(701);figure(701)
        freqset = [1:140];
        
        
        stimset = nearest(time,75):nearest(time,350);
        
        GrangerUpbl = mean(cat(1,GrangerUpOutAll(baselineset,freqset,:),GrangerUpInAll(baselineset,freqset,:)),1);
        GrangerDownbl = mean(cat(1,GrangerDownOutAll(baselineset,freqset,:),GrangerDownInAll(baselineset,freqset,:)),1);
        %         GrangerUpbl = mean(mean(cat(1,GrangerUpOutAll(baselineset,freqset,:),GrangerUpInAll(baselineset,freqset,:)),1),2);
        %         GrangerDownbl = mean(mean(cat(1,GrangerDownOutAll(baselineset,freqset,:),GrangerDownInAll(baselineset,freqset,:)),1),2);
        
        if combineconds
            GrangerUps = mean(cat(1,GrangerUpOutAll(stimset,freqset,:),GrangerUpInAll(stimset,freqset,:)),1);
            
            GrangerDowns = mean(cat(1,GrangerDownOutAll(stimset,freqset,:),GrangerDownInAll(stimset,freqset,:)),1);
            
            %
            % GrangerUpnormal = [];GrangerDownnormal = [];
            %         for i=1:size(GrangerUpOutgamma,3)
            % %             GrangerUpnormal(:,i) = squeeze(GrangerUps(1,:,i))./squeeze(GrangerUpbl(1,:,i));
            %             GrangerUpnormal(:,i)= squeeze(GrangerUps(1,:,i));
            %         end
            %         for i=1:size(GrangerDownOutgamma,3)
            % %             GrangerDownnormal(:,i)=squeeze(GrangerDowns(1,:,i))./squeeze(GrangerDownbl(1,:,i));
            %             GrangerDownnormal(:,i) = squeeze(GrangerDowns(1,:,i));
            %         end
            %
            %         plot(freqset,mean(GrangerDownnormal,2),'r')
            %         hold on
            %         dofill(freqset,GrangerDownnormal','r',1,1)
            %         plot(freqset, mean(GrangerUpnormal,2),'b')
            %         dofill(freqset,GrangerUpnormal','b',1,1)
            %         dopptt(abs(GrangerUpnormal),abs(GrangerDownnormal),freqset,1,1,0)
            %
            %         xlabel('freq (Hz)'),ylabel('raw GC')
            %         title(['all trials' num2str(fsets(fi,1)) '-' num2str(fsets(fi,2)) ' Hz, n = ' num2str(size(GrangerDownnormal,2))])
            %         %             legend('a to b', 'b to a','location','best')
            % %         line([-200 650],[1 1],'color','k');hold on;            line([0 0],[-15 15],'color','k')
            %         %             axis([-100 600 0.8 2])
            %         %             subplot(3,2,sub2ind([2 3],fi,layi))
            %         m1 = max(mean(abs(GrangerDownnormal),2));m2 = max(mean(abs(GrangerUpnormal),2));
            %         mval = max([m1 m2]);
            %         minval = min([mean(abs(GrangerDownnormal),2)' mean(abs(GrangerUpnormal),2)']);
            %
            %                 axis([1 150 minval-.1 mval+.1])
            %         subplot 122
            
            
            GrangerUpnormal = [];GrangerDownnormal = [];
            
            for i=1:size(GrangerUpOutgamma,3)
                if norm2bl
                    GrangerUpnormal(:,i) = squeeze(GrangerUps(1,:,i))./squeeze(GrangerUpbl(1,:,i));
                else
                    GrangerUpnormal(:,i)= squeeze(GrangerUps(1,:,i));
                end
            end
            for i=1:size(GrangerDownOutgamma,3)
                if norm2bl
                    GrangerDownnormal(:,i)=squeeze(GrangerDowns(1,:,i))./squeeze(GrangerDownbl(1,:,i));
                else
                    GrangerDownnormal(:,i) = squeeze(GrangerDowns(1,:,i));
                end
            end
            
            plot(freqset,nanmean(GrangerDownnormal,2),'r')
            hold on
            dofill(freqset,GrangerDownnormal','r',1,1)
            plot(freqset, nanmean(GrangerUpnormal,2),'b')
            dofill(freqset,GrangerUpnormal','b',1,1)
            %         dopptta(abs(GrangerUpnormal),abs(GrangerDownnormal),freqset,1,1,0,alpha)
            
            xlabel('freq (Hz)'),ylabel('normalized GC')
            title(['all trials, normalized' ' n = ' num2str(size(GrangerDownnormal,2))])
            %             legend('a to b', 'b to a','location','best')
            %         line([-200 650],[1 1],'color','k');hold on;            line([0 0],[-15 15],'color','k')
            %             axis([-100 600 0.8 2])
            %             subplot(3,2,sub2ind([2 3],fi,layi))
            m1 = max(mean(abs(GrangerDownnormal),2));m2 = max(mean(abs(GrangerUpnormal),2));
            mval = max([m1 m2]);
            minval = min([mean(abs(GrangerDownnormal),2)' mean(abs(GrangerUpnormal),2)']);
            %         axis([30 140 .95 1.5])
            %         box off
            %         axis off
            %         set(gca,'xtick',[30:10:140])
        else
            
            
            subplot 122
            GrangerUps = mean(cat(1,GrangerUpOutAll(stimset,freqset,:)),1);
            GrangerDowns = mean(cat(1,GrangerDownOutAll(stimset,freqset,:)),1);
            
            
            
            
            GrangerUpnormal = [];GrangerDownnormal = [];
            for i=1:size(GrangerUpOutgamma,3)
                if norm2bl
                    GrangerUpnormal(:,i) = squeeze(GrangerUps(1,:,i))./squeeze(GrangerUpbl(1,:,i));
                    ystr = 'normalized GC';
                else
                    GrangerUpnormal(:,i)= squeeze(GrangerUps(1,:,i));
                    ystr = 'raw GC';
                end
            end
            for i=1:size(GrangerDownOutgamma,3)
                if norm2bl
                    GrangerDownnormal(:,i)=squeeze(GrangerDowns(1,:,i))./squeeze(GrangerDownbl(1,:,i));
                else
                    GrangerDownnormal(:,i) = squeeze(GrangerDowns(1,:,i));
                end
            end
            
            plot(freqset,nanmean(GrangerDownnormal,2),'r')
            hold on
            dofill(freqset,GrangerDownnormal','r',1,1)
            plot(freqset, nanmean(GrangerUpnormal,2),'b')
            dofill(freqset,GrangerUpnormal','b',1,1)
            %         dopptta(abs(GrangerUpnormal),abs(GrangerDownnormal),freqset,1,1,0,alpha)
            
            xlabel('freq (Hz)'),ylabel(ystr)
            title(['all trials, normalized' ' n = ' num2str(size(GrangerDownnormal,2))])
            %             legend('a to b', 'b to a','location','best')
            %         line([-200 650],[1 1],'color','k');hold on;            line([0 0],[-15 15],'color','k')
            %             axis([-100 600 0.8 2])
            %             subplot(3,2,sub2ind([2 3],fi,layi))
            m1 = max(mean(abs(GrangerDownnormal),2));m2 = max(mean(abs(GrangerUpnormal),2));
            mval = max([m1 m2]);
            minval = min([mean(abs(GrangerDownnormal),2)' mean(abs(GrangerUpnormal),2)']);
            %         axis([30 140 .95 1.5])
            %         box off
            %         axis off
            %         set(gca,'xtick',[30:10:140])
            title('recognition')
            
            
            
            subplot 121
            GrangerUps = mean(cat(1,GrangerUpInAll(stimset,freqset,:)),1);
            GrangerDowns = mean(cat(1,GrangerDownInAll(stimset,freqset,:)),1);
            
            
            
            
            GrangerUpnormal = [];GrangerDownnormal = [];
            for i=1:size(GrangerUpOutgamma,3)
                if norm2bl
                    GrangerUpnormal(:,i) = squeeze(GrangerUps(1,:,i))./squeeze(GrangerUpbl(1,:,i));
                else
                    GrangerUpnormal(:,i)= squeeze(GrangerUps(1,:,i));
                end
            end
            for i=1:size(GrangerDownOutgamma,3)
                if norm2bl
                    GrangerDownnormal(:,i)=squeeze(GrangerDowns(1,:,i))./squeeze(GrangerDownbl(1,:,i));
                else
                    GrangerDownnormal(:,i) = squeeze(GrangerDowns(1,:,i));
                end
            end
            
            plot(freqset,nanmean(GrangerDownnormal,2),'r')
            hold on
            dofill(freqset,GrangerDownnormal','r',1,1)
            plot(freqset, nanmean(GrangerUpnormal,2),'b')
            dofill(freqset,GrangerUpnormal','b',1,1)
            %         dopptta(abs(GrangerUpnormal),abs(GrangerDownnormal),freqset,1,1,0,alpha)
            
            xlabel('freq (Hz)'),ylabel(ystr)
            title(['all trials, normalized' ' n = ' num2str(size(GrangerDownnormal,2))])
            %             legend('a to b', 'b to a','location','best')
            %         line([-200 650],[1 1],'color','k');hold on;            line([0 0],[-15 15],'color','k')
            %             axis([-100 600 0.8 2])
            %             subplot(3,2,sub2ind([2 3],fi,layi))
            m1 = max(mean(abs(GrangerDownnormal),2));m2 = max(mean(abs(GrangerUpnormal),2));
            mval = max([m1 m2]);
            minval = min([mean(abs(GrangerDownnormal),2)' mean(abs(GrangerUpnormal),2)']);
            %         axis([30 140 .95 1.5])
            %         box off
            %         axis off
            %         set(gca,'xtick',[30:10:140])
            title('encoding')
        end
        %%
        saveas(105,[figdir FT2 '_' num2str(105) '.png'])
        saveas(106,[figdir FT2 '_' num2str(106) '.png'])
        saveas(501,[figdir FT2 '_' num2str(501) '.png'])
        saveas(502,[figdir FT2 '_' num2str(502) '.png'])
        saveas(601,[figdir FT2 '_' num2str(601) '.png'])
        saveas(602,[figdir FT2 '_' num2str(602) '.png'])
        saveas(701,[figdir FT2 '_' num2str(701) '.png'])
        %         saveas(702,[figdir FT2 '_' num2str(702) '.png'])
        
    end
    
    
    
end
% saveas(101,[figdir FT '_' num2str(101) '.png'])
% saveas(102,[figdir FT '_' num2str(102) '.png'])
% saveas(103,[figdir FT '_' num2str(103) '.png'])
% saveas(104,[figdir FT '_' num2str(104) '.png'])
% saveas(105,[figdir FT '_' num2str(105) '.png'])